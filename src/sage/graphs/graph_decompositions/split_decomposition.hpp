/*
 *      Copyright (C) 2026 David Coudert <david.coudert@inria.fr>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *                 https://www.gnu.org/licenses/
 */

/**
 * @file split_decomposition.hpp
 * @brief Incremental split decomposition (split tree) of a connected graph,
 *        using graph-labelled trees.
 *
 * Follows the data structure and the algorithms of Section 5 of
 *   E. Gioan, C. Paul, M. Tedder, D. Corneil,
 *   "Practical and Efficient Split Decomposition via Graph-Labelled Trees",
 *   Algorithmica 69 (2014), doi:10.1007/s00453-013-9752-9.
 *
 * ## Overview
 *  The graph is processed along an LexBFS order x_1, ..., x_n (Algorithm 1).
 *  Vertices are inserted one at a time into the split tree ST(G) of the
 *  already inserted prefix; inserting x with neighbourhood S = N(x) costs time
 *  proportional to |S| (plus union-find overhead), because only the subtree
 *  T(S) spanning S is ever touched.
 *
 * ## Vocabulary (see the paper, Sections 2 and 4)
 *  - Graph-labelled tree (GLT): a tree whose internal nodes u carry a label
 *    graph G(u); the vertices of G(u) ("markers") are in bijection with the
 *    tree-edges incident to u. Leaves are the graph vertices.
 *  - Two leaves are adjacent in the graph iff the path joining them in the
 *    tree alternates tree-edges and label-edges ("accessibility graph").
 *  - Split tree = unique reduced GLT whose labels are prime, cliques or stars
 *    (reduced: no clique-join and no star-join possible, every node has >= 3
 *    markers).
 *  - State of a marker q with respect to S: it is a property of the part of
 *    the tree that lies beyond the tree-edge of q. Let L(q) be the leaves
 *    there and A(q) those accessible from q.
 *       empty   : S ∩ L(q) = ∅
 *       perfect : S ∩ L(q) = A(q)
 *       mixed   : otherwise
 *    A tree-edge (q, opposite q) is "perfect/empty/mixed" according to the
 *    states of its two extremities.
 *  - Theorem 4.14 splits the insertion of x into 7 cases, depending on where
 *    the non-mixed nodes/edges are. The code identifies the case in insert()
 *    and performs the update.
 *
 * ## Data structure (Section 5.1)
 *  - The tree is rooted at the first leaf root_record. Every record (node or
 *    leaf) has a "root marker" `root_marker`, the extremity of its edge to its
 *    parent.
 *  - Parent pointers are stored in a union-find structure over the children of
 *    prime nodes (the "children-set" of a prime node is one set whose 
 *    representative points to the node). This makes a node-join O(1) in the
 *    number of children moved, which is what yields the near-linear bound.
 *    Children of degenerate (clique/star) nodes are kept as singletons.
 *  - Label graphs: cliques and stars are implicit (marker count + centre);
 *    only prime nodes (and nodes freshly joined, Remark 5.2) store explicit
 *    adjacency sets.
 *
 * ## Node-joins (Algorithms 8 and 9)
 *  - Contraction of the fully-mixed subtree runs in the paper's three phases
 *    (see method insert()): (1) stars whose root marker is their centre absorb
 *    their children, (2) nodes whose root marker has degree one are absorbed
 *    by their parent, (3) everything else, bottom-up. This order is what
 *    bounds the number of label edges created (Lemma 6.16 of the paper).
 *  - nodeJoin() implements the shortcut of Algorithm 8: a star whose root
 *    marker has degree one is merged without creating any new label edge
 *    (method nodeJoinStarShortcut()).
 *
 * ## Deviations from the paper
 *  - Adjacency of prime nodes uses std::unordered_set (expected O(1)) instead
 *    of adjacency lists.
 *  - The input order MUST be an LexBFS order (the twin test in step 3 of
 *    insert() relies on it). The validity of the given order is not checked.
 */

#include <algorithm>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>
using namespace std;


/**
 * @brief Returns the number of '1' bits in a 32-bits integer.
 *
 * This is used in method verify which limits the number of vertices in a prime
 * node to 14 (otherwise the verification is too slow).
 * This is a standard bit-twiddling hack, see
 * https://graphics.stanford.edu/~seander/bithacks.html#CountBitsSetParallel
 * Also used in fast_digraph.pyx.
 */
static inline int popcount32(unsigned x) noexcept {
    x = x - ((x >> 1) & 0x55555555);
    x = (x & 0x33333333) + ((x >> 2) & 0x33333333);
    return ((x + (x >> 4) & 0x0F0F0F0F) * 0x01010101) >> 24;
}

/**
 * @brief The graph-labelled tree together with the vertex-insertion algorithm.
 *
 * All identifiers are plain ints indexing the vectors below:
 *  - "marker" ids index opposite, next_marker, previous_marker,
 *    leaf_vertex_of_marker, record_rooted_at_marker, perfect_stamp,
 *    scratch_stamp, label_adjacency;
 *  - "record" ids (nodes and leaves alike) index records and scratch;
 *  - "element" ids index the union-find arrays uf_parent, uf_rank,
 *    uf_parent_record.
 */
struct FastGLT {
    /**
     * @brief Kind of a record: a node or a leaf of the tree.
     * LEAF records are the graph vertices; the others label internal nodes.
     */
    enum Kind { PRIME, CLIQUE, STAR, LEAF };

    /**
     * @brief A node or a leaf of the tree.
     *
     * The markers of a record form an intrusive doubly-linked list
     * (first_marker/last_marker, next_marker/previous_marker) so that markers
     * can be moved between records in O(1) each and two lists can be
     * concatenated in O(1).
     */
    struct Record {
        Kind kind = PRIME;
        int star_center = -1;   // star centre marker
        int first_marker = -1;  // first marker in the list
        int last_marker = -1;   // last marker in the list
        int marker_count = 0;   // number of markers in the list
        int root_marker = -1;   // root marker: extremity of the tree-edge to the parent (-1 for root_record)
        int union_find_element = -1;  // union-find element standing for this record as a child
        int children_set_element = -1;  // prime nodes: an element of the set containing all its children
        int last_inserted_marker = -1;  // prime nodes: the marker inserted last (used to find the universal marker)
        int universal_marker = -1;  // prime nodes: the marker adjacent to all other markers (or -1 if none)
        bool is_alive = true;  // false once the record has been absorbed by a node-join
    };

    /**
     * @brief Per-insertion scratch data attached to a record.
     *
     * Fields compared with `insertion_number` are timestamps: they are valid
     * only if equal to the current insertion number, which avoids clearing
     * anything between insertions (only visited records are re-initialised).
     */
    struct InsertionScratch {
        // ---- stamps for membership tests and visited flags ----
        long visited_stamp = -1;
        long active_stamp = -1;
        long in_neighbourhood_stamp = -1;
        long peeled_stamp = -1;
        long root_marker_perfect_stamp = -1;
        long in_contraction_stamp = -1;
        // ---- parent in T(S); counters used by the peeling phases ----
        int parent_in_subtree = -1;  // parent in T(S)
        int pending_child_count = 0;
        int peeled_child_count = 0;
        int node_child_count = 0;
        vector<int> children_in_subtree;  // children in T(S)
        vector<int> perfect_markers;  // perfect markers of this node
    };

    // ---- marker arrays ----
    // opposite[m]: marker at the other extremity of m's tree-edge
    vector<int> opposite;
    // links of the per-record marker list
    vector<int> next_marker, previous_marker;
    // leaf_vertex_of_marker[m]: graph vertex if m is the marker of a leaf, else -1
    vector<int> leaf_vertex_of_marker;
    // record_rooted_at_marker[m] = record R if m == R.root_marker (child
    // reached through opposite[m]); -1 otherwise
    vector<int> record_rooted_at_marker;
    // perfect_stamp[m] == insertion_number  <=>  m is perfect for the current
    // insertion
    vector<long> perfect_stamp;
    // scratch stamps (membership tests, see `stamp_counter`)
    vector<long> scratch_stamp;
    // label adjacency of non-degenerate nodes (only filled for PRIME records)
    vector<unordered_set<int>> label_adjacency;
    // ---- records ----
    vector<Record> records;
    vector<InsertionScratch> scratch;
    // ---- union-find over tree positions ----
    //    uf_parent_record[representative] = parent record
    vector<int> uf_parent, uf_rank, uf_parent_record;
    // leaf_record_of_vertex: graph vertex -> its leaf record
    vector<int> leaf_record_of_vertex;
    // root_record: root of the tree: leaf of the first inserted vertex
    int root_record = -1;
    // insertion_number: insertion counter
    long insertion_number = 0;
    // stamp_counter: counter for scratch_stamp values
    long stamp_counter = 0;

    // ====================================================================================
    //  Low-level structure: union-find, markers, records, marker lists
    // ====================================================================================

    /**
     * @brief Creates a fresh union-find singleton (no parent yet).
     * @return the id of the new union-find element
     */
    int newUnionFindElement() {
        int element = uf_parent.size();
        uf_parent.push_back(element);
        uf_rank.push_back(0);
        uf_parent_record.push_back(-1);
        return element;
    }

    /**
     * @brief Union-find `findRoot` with path halving.
     * @param element  the element whose set representative is queried
     * @return the representative of the set containing `element`
     */
    int findRoot(int element) {
        while (uf_parent[element] != element) {
            uf_parent[element] = uf_parent[uf_parent[element]];
            element = uf_parent[element];
        }
        return element;
    }

    /**
     * @brief Union by rank of the sets of a and b.
     * 
     * NOTE: the caller must then set uf_parent_record[returned rep] to the
     * parent record, since `uf_parent_record` is only meaningful on
     * representatives.
     * 
     * @param set_a  union-find element of the first set
     * @param set_b  union-find element of the second set
     * @return the new representative of the union of the two sets
     */
    int uniteSets(int set_a, int set_b) {
        set_a = findRoot(set_a);
        set_b = findRoot(set_b);
        if (set_a == set_b)
            return set_a;
        if (uf_rank[set_a] < uf_rank[set_b])
            swap(set_a, set_b);
        uf_parent[set_b] = set_a;
        if (uf_rank[set_a] == uf_rank[set_b])
            ++uf_rank[set_a];
        return set_a;
    }

    /**
     * @brief Allocates a new marker (not linked to any record or tree-edge
     *        yet).
     * @return the id of the new marker
     */
    int newMarker() {
        opposite.push_back(-1);
        next_marker.push_back(-1);
        previous_marker.push_back(-1);
        leaf_vertex_of_marker.push_back(-1);
        record_rooted_at_marker.push_back(-1);
        perfect_stamp.push_back(-1);
        scratch_stamp.push_back(-1);
        label_adjacency.emplace_back();
        return (int)opposite.size() - 1;
    }

    /**
     * @brief Allocates a new record of kind record_kind, with its own
     *        union-find element (O(1)).
     * @param record_kind  the kind of the new record
     * @return the id of the new record
     */
    int newRecord(Kind record_kind) {
        int record_id = records.size();
        records.emplace_back();
        records.back().kind = record_kind;
        scratch.emplace_back();
        records[record_id].union_find_element = newUnionFindElement();
        return record_id;
    }

    /**
     * @brief Makes a and b the two extremities of one tree-edge (O(1)).
     * @param marker_a  the first marker
     * @param marker_b  the second marker
     */
    void link(int marker_a, int marker_b) {
        opposite[marker_a] = marker_b;
        opposite[marker_b] = marker_a;
    }

    /**
     * @brief Appends marker at the end of the marker list of record R (O(1)).
     * @param record_id  the id of the record
     * @param marker  the marker to append
     */
    void appendMarker(int record_id, int marker) {
        Record& record = records[record_id];
        previous_marker[marker] = record.last_marker;
        next_marker[marker] = -1;
        if (record.last_marker >= 0) {
            next_marker[record.last_marker] = marker;
        } else {
            record.first_marker = marker;
        }
        record.last_marker = marker;
        ++record.marker_count;
    }

    /**
     * @brief Removes marker from the marker list of record R (O(1)).
     * @param record_id  the id of the record
     * @param marker  the marker to remove
     */
    void removeMarker(int record_id, int marker) {
        Record& record = records[record_id];
        if (previous_marker[marker] >= 0) {
            next_marker[previous_marker[marker]] = next_marker[marker];
        } else {
            record.first_marker = next_marker[marker];
        }
        if (next_marker[marker] >= 0) {
            previous_marker[next_marker[marker]] = previous_marker[marker];
        } else {
            record.last_marker = previous_marker[marker];
        }
        previous_marker[marker] = next_marker[marker] = -1;
        --record.marker_count;
    }

    /**
     * @brief Returns the markers of record R as a vector (O(size); used by
     *        non-time-critical code only).
     * @param record_id  the id of the record
     * @return a vector containing the markers of the record
     */
    vector<int> markersOf(int record_id) const {
        vector<int> result;
        for (int marker = records[record_id].first_marker; marker >= 0; marker = next_marker[marker])
            result.push_back(marker);
        return result;
    }

    // ====================================================================================
    //  Tree navigation
    // ====================================================================================

    /**
     * @brief Returns the parent record of record_id (-1 for the root
     *        root_record). Amortised near-constant thanks to union-find.
     * @param record_id  the id of the record
     * @return the id of the parent record
     */
    int parentRecord(int record_id) {
        return record_id == root_record ? -1 : uf_parent_record[findRoot(records[record_id].union_find_element)];
    }

    /**
     * @brief Returns the record at the other end of the tree-edge of marker
     *        q, where q belongs to record u: the parent of u if q is u's root
     *        marker, otherwise the child hanging below q.
     * @param node  the id of the record u
     * @param marker  the marker q belonging to u
     * @return the id of the record at the other end of q's tree-edge
     */
    int recordAcrossEdge(int node, int marker) {
        return marker == records[node].root_marker ? parentRecord(node) : record_rooted_at_marker[opposite[marker]];
    }

    /**
     * @brief Declares the fresh record `child` (whose union-find element is
     *        still a singleton) a child of record `parent`. For a prime
     *        parent, the child joins the single children-set of the parent;
     *        for a degenerate parent it stays a singleton.
     * @param child  the id of the child record
     * @param parent  the id of the parent record
     */
    void placeChild(int child, int parent) {
        int child_element = records[child].union_find_element;
        Record& record_parent = records[parent];
        if (record_parent.kind == PRIME) {
            int representative = uniteSets(record_parent.children_set_element, child_element);
            uf_parent_record[representative] = parent;
            record_parent.children_set_element = representative;
        } else {
            uf_parent_record[findRoot(child_element)] = parent;
        }
    }

    /**
     * @brief Creates the leaf record of graph vertex `vertex` (one marker, the
     *        leaf's own marker, which is its root marker).
     * @param vertex  the id of the graph vertex
     * @return the id of the new leaf record
     */
    int newLeaf(int vertex) {
        int leaf = newRecord(LEAF);
        int marker = newMarker();
        leaf_vertex_of_marker[marker] = vertex;
        appendMarker(leaf, marker);
        records[leaf].root_marker = marker;
        if ((int)leaf_record_of_vertex.size() <= vertex)
            leaf_record_of_vertex.resize(vertex + 1, -1);
        leaf_record_of_vertex[vertex] = leaf;
        return leaf;
    }

    /**
     * @brief Creates the leaf of the new vertex x and hangs it as a child of
     *        record P through P's marker m.
     * @param parent  the id of the parent record
     * @param parent_marker  the marker of the parent record
     * @param new_vertex  the id of the new graph vertex
     */
    void hangLeaf(int parent, int parent_marker, int new_vertex) {
        int leaf = newLeaf(new_vertex);
        link(parent_marker, records[leaf].root_marker);
        record_rooted_at_marker[records[leaf].root_marker] = leaf;
        placeChild(leaf, parent);
    }

    // ====================================================================================
    //  Label graphs
    // ====================================================================================

    /**
     * @brief Are markers a and b of record R adjacent in the label graph G(R)?
     * 
     * Cliques: always; stars: iff one of them is the centre; prime (or
     * joined) nodes: adjacency set.
     * 
     * @param record_id  the record whose label graph is considered
     * @param marker_a   first marker
     * @param marker_b   second marker
     * @return true if a and b are adjacent in G(R), false otherwise
     */
    bool areAdjacent(int record_id, int marker_a, int marker_b) const {
        if (marker_a == marker_b) return false;
        switch (records[record_id].kind) {
            case CLIQUE: return true;
            case STAR:   return marker_a == records[record_id].star_center || marker_b == records[record_id].star_center;
            case PRIME:  return label_adjacency[marker_a].count(marker_b) > 0;
            default:     return false;
        }
    }

    /**
     * @brief Degree of marker m in the label graph of record record_id.
     * 
     * @param record_id  the record whose label graph is considered
     * @param marker     the marker whose degree is queried
     * @return the degree of marker m in the label graph of record record_id
     */
    int labelDegree(int record_id, int marker) const {
        switch (records[record_id].kind) {
            case CLIQUE: return records[record_id].marker_count - 1;
            case STAR:   return marker == records[record_id].star_center ? records[record_id].marker_count - 1 : 1;
            case PRIME:  return (int)label_adjacency[marker].size();
            default:     return 0;
        }
    }

    /**
     * @brief Converts a degenerate node to an explicit PRIME-kind node
     *        (adjacency sets filled in).
     * 
     * Needed before a node-join, whose result is generally neither a clique
     * nor a star (Remark 5.2). Only called by nodeJoin(); the label graph
     * is unchanged.
     */
    void materialize(int record_id) {
        if (records[record_id].kind == PRIME) return;
        auto all_markers = markersOf(record_id);
        for (int marker : all_markers)
            label_adjacency[marker].clear();
        if (records[record_id].kind == CLIQUE) {
            for (int marker_a : all_markers)
                for (int marker_b : all_markers)
                    if (marker_a != marker_b)
                        label_adjacency[marker_a].insert(marker_b);
        } else {
            int centre = records[record_id].star_center;
            for (int marker : all_markers)
                if (marker != centre) {
                    label_adjacency[centre].insert(marker);
                    label_adjacency[marker].insert(centre);
                }
        }
        records[record_id].kind = PRIME;
        records[record_id].star_center = -1;
    }

    /**
     * @brief For a degenerate node (children are singletons in the
     *        union-find), merges all its children into one children-set,
     *        as prime nodes require. Precondition: R has at least one child.
     * @param record_id  the id of the degenerate node
     */
    void formChildrenSet(int record_id) {
        int representative = -1;
        for (int marker = records[record_id].first_marker; marker >= 0; marker = next_marker[marker]) {
            if (marker == records[record_id].root_marker) continue;
            int child = record_rooted_at_marker[opposite[marker]];
            int child_element = findRoot(records[child].union_find_element);
            representative = representative < 0 ? child_element : uniteSets(representative, child_element);
        }
        uf_parent_record[representative] = record_id;
        records[record_id].children_set_element = representative;
    }

    // ====================================================================================
    //  Node-split (Algorithm 6) and node-join (Algorithm 8)
    // ====================================================================================

    /**
     * @brief Node-split of the degenerate node `node` along the bipartition
     *        (moved_markers, markersOf(node) \ moved_markers).
     *
     * Creates a new node new_node of the same kind (clique or star) holding
     * the moved_markers, plus a new tree-edge (node_marker, new_node_marker)
     * with node_marker in `node` and new_node_marker in new_node. The caller
     * guarantees
     * |moved_markers| >= 2 and |markersOf(node) \ moved_markers| >= 2
     * so that both halves have >= 3 markers.
     * For a star, the centre c stays where it is; the other half receives a
     * new centre marker (q or r) so that both halves remain stars and the
     * accessibility graph is unchanged.
     * Cost O(|moved_markers|) plus union-find updates (their children are
     *  re-parented to new_node).
     *
     * @param node           degenerate node to split
     * @param moved_markers markers moved to the new node
     * @param node_marker_out     [out] new marker in `node` (extremity of the new tree-edge)
     * @param new_node_marker_out [out] new marker in the new node (other extremity)
     * @return the new node
     */
    int nodeSplit(int node, const vector<int>& moved_markers,
                  int& node_marker_out, int& new_node_marker_out) {
        Kind node_kind = records[node].kind;
        int centre = records[node].star_center;
        long stamp = ++stamp_counter;
        for (int moved : moved_markers)
            scratch_stamp[moved] = stamp;      // mark membership in A
        bool root_marker_moved = records[node].root_marker >= 0 && scratch_stamp[records[node].root_marker] == stamp;  // is the edge to the parent on the A side?
        int old_parent = parentRecord(node);
        int new_node = newRecord(node_kind);
        for (int moved : moved_markers) {
            removeMarker(node, moved);
            appendMarker(new_node, moved);
        }
        int new_node_marker = newMarker();
        int node_marker = newMarker();
        appendMarker(new_node, new_node_marker);
        appendMarker(node, node_marker);
        link(node_marker, new_node_marker);
        if (node_kind == STAR) {
            // the half containing the old centre keeps it; the other half
            // gets its new marker as centre
            if (scratch_stamp[centre] == stamp) {
                records[new_node].star_center = centre;
                records[node].star_center = node_marker;
            } else {
                records[new_node].star_center = new_node_marker;
                records[node].star_center = centre;
            }
        }
        // children of the moved markers now hang below new_node
        for (int moved : moved_markers)
            if (moved != records[node].root_marker) {
                int child = record_rooted_at_marker[opposite[moved]];
                uf_parent_record[findRoot(records[child].union_find_element)] = new_node;
            }
        if (root_marker_moved) {
            // new_node takes node's place in the tree, node becomes a child
            // of new_node
            records[new_node].root_marker = records[node].root_marker;
            record_rooted_at_marker[records[new_node].root_marker] = new_node;
            records[new_node].union_find_element = newUnionFindElement();
            placeChild(new_node, old_parent);
            records[node].root_marker = node_marker;
            record_rooted_at_marker[node_marker] = node;
            records[node].union_find_element = newUnionFindElement();
            uf_parent_record[records[node].union_find_element] = new_node;
        } else {  // new_node becomes a child of node
            records[new_node].root_marker = new_node_marker;
            record_rooted_at_marker[new_node_marker] = new_node;
            records[new_node].union_find_element = newUnionFindElement();
            uf_parent_record[records[new_node].union_find_element] = node;
        }
        node_marker_out = node_marker;
        new_node_marker_out = new_node_marker;
        return new_node;
    }

    /**
     * @brief Node-join: contracts the tree-edge between node child_node and
     *        its parent parent_node. parent_node is kept.
     *
     * The new label is the join of the two label graphs along the edge: every
     * neighbour of parent_marker in the parent's label becomes adjacent to
     * every neighbour of child_root_marker in the child's label
     * (parent_marker and child_root_marker are deleted).
     * The result is stored as PRIME-kind with explicit adjacency, even if it
     * happens to be degenerate (Remark 5.2: reduction is performed lazily by
     * the insertion cases). The marker lists are concatenated and the
     * children-sets united, both in O(1)-ish.
     *
     * @param parent_node parent node (kept)
     * @param child_node  child node (absorbed, marked dead)
     */
    void nodeJoin(int parent_node, int child_node) {
        int child_root_marker = records[child_node].root_marker;
        int parent_marker = opposite[child_root_marker];
        if (records[child_node].kind == STAR && labelDegree(child_node, child_root_marker) == 1) {
            nodeJoinStarShortcut(parent_node, child_node);
            return;
        }
        if (records[parent_node].kind != PRIME) {
            formChildrenSet(parent_node);
            materialize(parent_node);
        }
        if (records[child_node].kind != PRIME) {
            formChildrenSet(child_node);
            materialize(child_node);
        }
        // rewire label edges: N(q) x N(r) become edges, then delete q and r
        vector<int> parent_marker_neighbours(label_adjacency[parent_marker].begin(),
                                             label_adjacency[parent_marker].end());
        vector<int> child_marker_neighbours(label_adjacency[child_root_marker].begin(),
                                            label_adjacency[child_root_marker].end());
        for (int neighbour_a : parent_marker_neighbours)
            label_adjacency[neighbour_a].erase(parent_marker);
        for (int neighbour_b : child_marker_neighbours)
            label_adjacency[neighbour_b].erase(child_root_marker);
        label_adjacency[parent_marker].clear();
        label_adjacency[child_root_marker].clear();
        for (int neighbour_a : parent_marker_neighbours)
            for (int neighbour_b : child_marker_neighbours) {
                label_adjacency[neighbour_a].insert(neighbour_b);
                label_adjacency[neighbour_b].insert(neighbour_a);
            }
        removeMarker(parent_node, parent_marker);
        removeMarker(child_node, child_root_marker);
        Record &parent_record = records[parent_node];
        Record &child_record = records[child_node];
        if (child_record.first_marker >= 0) {  // O(1) concatenation of the marker lists
            if (parent_record.last_marker >= 0)
                next_marker[parent_record.last_marker] = child_record.first_marker;
            else
                parent_record.first_marker = child_record.first_marker;
            previous_marker[child_record.first_marker] = parent_record.last_marker;
            parent_record.last_marker = child_record.last_marker;
            parent_record.marker_count += child_record.marker_count;
        }
        int representative = uniteSets(parent_record.children_set_element, child_record.children_set_element);
        uf_parent_record[representative] = parent_node;
        parent_record.children_set_element = representative;  // children of v become children of u
        record_rooted_at_marker[child_root_marker] = -1;
        child_record.is_alive = false;
        child_record.first_marker = child_record.last_marker = -1;
        child_record.marker_count = 0;
        parent_record.kind = PRIME;
        parent_record.star_center = -1;
    }

    /**
     * @brief Node-join shortcut (Algorithm 8, step 1): child_node is a star
     *        whose root marker r has degree one.
     *
     * Then r is adjacent only to the star's centre c, and joining along
     * (parent_marker, r) would make c adjacent to the neighbours of
     * parent_marker. Instead of creating those label edges, the marker
     * parent_marker itself takes the place of c: the tree-edge that was
     * attached to c is re-attached to parent_marker (which already has the
     * right neighbours in the parent's label), the other leaves of the star
     * are moved into the parent and made adjacent to parent_marker only.
     * No label edge is created except one per moved marker, and the cost is
     * O(#moved markers) union-find operations.
     *
     * @param parent_node parent node (kept; converted to explicit adjacency if degenerate)
     * @param child_node  star to absorb (marked dead)
     */
    void nodeJoinStarShortcut(int parent_node, int child_node) {
        int child_root_marker = records[child_node].root_marker;
        int parent_marker = opposite[child_root_marker];
        int centre = records[child_node].star_center;
        if (records[parent_node].kind != PRIME) {
            formChildrenSet(parent_node);
            materialize(parent_node);
        }
        // the tree-edge at the centre is re-attached to parent_marker
        int below_centre = opposite[centre];
        int centre_child = record_rooted_at_marker[below_centre];
        link(parent_marker, below_centre);
        opposite[child_root_marker] = -1;
        vector<int> moved_children{centre_child};
        // the remaining leaves of the star move to the parent,
        // adjacent to parent_marker only
        for (int marker : markersOf(child_node)) {
            if (marker == child_root_marker || marker == centre)
                continue;
            removeMarker(child_node, marker);
            appendMarker(parent_node, marker);
            label_adjacency[marker].clear();
            label_adjacency[marker].insert(parent_marker);
            label_adjacency[parent_marker].insert(marker);
            moved_children.push_back(record_rooted_at_marker[opposite[marker]]);
        }
        // all children of the star become children of the parent
        for (int child : moved_children) {
            int representative = uniteSets(records[parent_node].children_set_element, records[child].union_find_element);
            uf_parent_record[representative] = parent_node;
            records[parent_node].children_set_element = representative;
        }
        record_rooted_at_marker[child_root_marker] = -1;
        Record& child_record = records[child_node];
        child_record.is_alive = false;
        child_record.first_marker = child_record.last_marker = -1;
        child_record.marker_count = 0;
        records[parent_node].kind = PRIME;
        records[parent_node].star_center = -1;
    }

    // ====================================================================================
    //  Elementary insertions (Propositions 4.15 and 4.17)
    // ====================================================================================

    /**
     * @brief Adds the leaf of new_vertex to the degenerate node `node` (clique
     *        or star) (Prop. 4.15, cases 1-2).
     * 
     * For a clique the new marker is adjacent to all; for a star it is a new
     * degree-1 marker (adjacent to the centre).
     * 
     * @param node        the degenerate node (clique or star)
     * @param new_vertex  the vertex being inserted
     */
    void attachDegenerate(int node, int new_vertex) {
        int marker = newMarker();
        appendMarker(node, marker);
        hangLeaf(node, marker, new_vertex);
    }

    /**
     * @brief Adds the leaf of new_vertex to the prime (or freshly joined)
     *        `node` (Prop. 4.15, prime case, and the final step of case 7).
     *
     * A new marker is added to the label of `node`, adjacent exactly to
     * adjacent_markers (the perfect markers). It also maintains
     * records[node].last_inserted_marker (marker inserted last) and
     * records[node].universal_marker (universal marker, or -1), which are the
     * two candidates for the LexBFS twin test in insert() step 3.
     * 
     * @param node              the prime node (or freshly joined node)
     * @param adjacent_markers  the markers of `node` adjacent to the new marker
     * @param new_vertex        the vertex being inserted
     */
    void attachPrime(int node, const vector<int>& adjacent_markers, int new_vertex) {
        int old_marker_count = records[node].marker_count;
        int new_marker = newMarker();
        label_adjacency[new_marker].insert(adjacent_markers.begin(), adjacent_markers.end());
        for (int adjacent : adjacent_markers)
            label_adjacency[adjacent].insert(new_marker);
        long stamp = ++stamp_counter;
        for (int adjacent : adjacent_markers)
            scratch_stamp[adjacent] = stamp;
        int old_universal = records[node].universal_marker;
        int new_universal = -1;
        // the new marker is universal iff it is adjacent to every old marker;
        // an old universal marker stays universal iff it is adjacent to the new one
        if ((int)adjacent_markers.size() == old_marker_count)
            new_universal = new_marker;
        else if (old_universal >= 0 && scratch_stamp[old_universal] == stamp)
            new_universal = old_universal;
        appendMarker(node, new_marker);
        records[node].last_inserted_marker = new_marker;
        records[node].universal_marker = new_universal;
        hangLeaf(node, new_marker, new_vertex);
    }

    /**
     * @brief Subdivides the tree-edge (marker, opposite[marker]) by a new
     *        degree-3 node holding the leaf of new_vertex (Prop. 4.17).
     *
     * The new node is a clique if both extremities are perfect, otherwise a
     * star whose centre faces the perfect side (so the new vertex becomes
     * adjacent only to what is accessible through the perfect side).
     * Also used as the second step of case 4 (Prop. 4.16) after a node-split.
     *
     * @param node   record owning `marker`
     * @param marker marker of `node` whose tree-edge is subdivided
     * @param marker_perfect state of `marker` is perfect (otherwise empty)
     * @param opposite_perfect state of opposite[marker] is perfect (otherwise empty)
     * @param new_vertex the vertex being inserted
     */
    void subdivide(int node, int marker, bool marker_perfect, bool opposite_perfect, int new_vertex) {
        int opposite_marker = opposite[marker];
        int other_node = recordAcrossEdge(node, marker);
        // lower_node = the node whose marker is the lower endpoint of the edge (the child);
        // upper_node = the other node (the parent);
        // lower_marker = the marker at the lower endpoint
        // upper_marker = the other marker;
        // lower_perfect = the state of lower_marker
        // upper_perfect = the state of upper_marker;
        int lower_node, upper_node, lower_marker, upper_marker;
        bool lower_perfect, upper_perfect;
        if (marker == records[node].root_marker) {
            lower_node = node;
            lower_marker = marker;
            upper_node = other_node;
            upper_marker = opposite_marker;
            lower_perfect = marker_perfect;
            upper_perfect = opposite_perfect;
        }
        else {
            lower_node = other_node;
            lower_marker = opposite_marker;
            upper_node = node;
            upper_marker = marker;
            lower_perfect = opposite_perfect;
            upper_perfect = marker_perfect;
        }
        bool is_clique = lower_perfect && upper_perfect;
        int new_node = newRecord(is_clique ? CLIQUE : STAR);
        int new_up_marker = newMarker();     // towards parent
        int new_down_marker = newMarker();   // towards child
        int new_leaf_marker = newMarker();   // towards x
        appendMarker(new_node, new_up_marker);
        appendMarker(new_node, new_down_marker);
        appendMarker(new_node, new_leaf_marker);
        if (!is_clique)
            records[new_node].star_center = !upper_perfect ? new_up_marker : new_down_marker;       // centre faces the perfect side
        link(upper_marker, new_up_marker);
        link(lower_marker, new_down_marker);
        records[new_node].root_marker = new_up_marker;
        record_rooted_at_marker[new_up_marker] = new_node;
        placeChild(new_node, upper_node);   // new node becomes a child of the old parent
        records[lower_node].union_find_element = newUnionFindElement();
        uf_parent_record[records[lower_node].union_find_element] = new_node;   // lower_node becomes a child of the new node
        hangLeaf(new_node, new_leaf_marker, new_vertex);
    }

    /**
     * @brief Performs the insertion when a tree-edge with no mixed extremity
     *        has been identified (cases 1, 2, 5, 6 of Theorem 4.14).
     *
     * The edge is (marker, opposite[marker]) with marker in `node`. Both
     * extremities are perfect (PP) or one is perfect and the other empty (PE).
     * The new leaf is attached to an adjacent degenerate node when this keeps
     * the tree reduced (cases 1, 2), otherwise the edge is subdivided
     * (cases 5, 6).
     * 
     * @param node            record owning `marker`
     * @param marker          marker of `node` whose tree-edge is considered
     * @param marker_perfect  state of `marker` is perfect (otherwise empty)
     * @param opposite_perfect state of opposite[marker] is perfect (otherwise empty)
     * @param new_vertex      the vertex being inserted
     */
    void edgeCase(int node, int marker, bool marker_perfect, bool opposite_perfect, int new_vertex) {
        int opposite_marker = opposite[marker];
        int other_node = recordAcrossEdge(node, marker);
        if (marker_perfect && opposite_perfect) {                              // PP: case 1 or 5
            if (records[node].kind == CLIQUE) attachDegenerate(node, new_vertex);
            else if (records[other_node].kind == CLIQUE) attachDegenerate(other_node, new_vertex);
            else subdivide(node, marker, true, true, new_vertex);
            return;
        }
        // PE: one extremity is perfect, the other empty (case 2 or 6)
        int perfect_node, perfect_marker, empty_node, empty_marker;
        if (marker_perfect) {
            perfect_node = node;
            perfect_marker = marker;
            empty_node = other_node;
            empty_marker = opposite_marker;
        } else {
            perfect_node = other_node;
            perfect_marker = opposite_marker;
            empty_node = node;
            empty_marker = marker;
        }
        if (records[perfect_node].kind == STAR && records[perfect_node].star_center == perfect_marker)
            attachDegenerate(perfect_node, new_vertex);  // case 2, star centred on the perfect marker
        else if (records[empty_node].kind == STAR && records[empty_node].star_center != empty_marker)
            attachDegenerate(empty_node, new_vertex);    // case 2, star whose empty marker is a leaf of the star
        else
            subdivide(node, marker, marker_perfect, opposite_perfect, new_vertex);  // case 6
    }

    // ====================================================================================
    //  Initial vertices
    // ====================================================================================

    /**
     * @brief Adds the first vertex to the tree (Algorithm 3, step 1).
     * The tree is a single leaf, which becomes the root root_record.
     * @param vertex  the id of the first vertex
     */
    void addFirst(int vertex) {
        root_record = newLeaf(vertex);
        records[root_record].root_marker = -1;
    }

    /**
     * @brief Adds the second vertex to the tree (Algorithm 3, step 2).
     * The tree is two adjacent leaves (no internal node yet; nodes appear
     * from 3 vertices).
     * @param vertex  the id of the second vertex
     */
    void addSecond(int vertex) {
        int leaf = newLeaf(vertex);
        link(records[root_record].first_marker, records[leaf].root_marker);
        record_rooted_at_marker[records[leaf].root_marker] = leaf;
        uf_parent_record[records[leaf].union_find_element] = root_record;
    }

    // ====================================================================================
    //  Vertex insertion (Algorithm 3), driven by Theorem 4.14
    // ====================================================================================

    /**
     * @brief T' = T(S) after trimming: the records of the subtree spanned by
     *        the neighbourhood S
     * @param vertex  the id of the new vertex
     */
    struct Subtree {
        // all records of T' (stamp `visited_stamp` == insertion_number),
        // parents before children not guaranteed
        vector<int> records;
        // the record of T' closest to the root (the lowest common ancestor of S)
        int top = -1;
    };

    /**
     * @brief (Re)initialises the scratch data of `node` on its first visit
     *        during the current
     * @param node  the record being visited
     * @param visited_records  the list of records that have been visited
     */
    void markVisited(int node, vector<int>& visited_records) {
        InsertionScratch& scratch_data = scratch[node];
        scratch_data.visited_stamp = insertion_number;
        scratch_data.children_in_subtree.clear();
        scratch_data.perfect_markers.clear();
        scratch_data.pending_child_count = scratch_data.peeled_child_count = scratch_data.node_child_count = 0;
        visited_records.push_back(node);
    }

    /**
     * @brief One climbing step of Algorithm 4: visits `node`, and schedules
     *        its parent (once)
     * @param node  the record being visited
     * @param visited_records  the list of records that have been visited
     * @param next_active  the list of records that will be visited in the next iteration
     * @param root_reached  indicates whether the root has been reached
     */
    void climbOneLevel(int node, vector<int>& visited_records, vector<int>& next_active, bool& root_reached) {
        markVisited(node, visited_records);
        if (node == root_record) {
            root_reached = true;
            return;
        }
        int parent = parentRecord(node);
        scratch[node].parent_in_subtree = parent;
        if (scratch[parent].visited_stamp != insertion_number && scratch[parent].active_stamp != insertion_number) {
            scratch[parent].active_stamp = insertion_number;
            next_active.push_back(parent);
        }
    }

    /**
     * @brief Algorithm 4: computes T' = T(S) by synchronous climbing from the
     *        leaves of S.
     *
     * `active_records` holds marked records whose parent is not marked yet.
     * Each round, every active record is visited and marks its parent (stamp
     * `active_stamp` avoids duplicates). Climbing continues while at least two
     * active records remain (they must meet at the top) or, once the root is
     * reached, while at least one remains. Total work is O(|T(S)|).
     * Afterwards the children lists of T' are filled, its top is found and
     * trimmed (a top with a single child that is not in S is not part of T(S)).
     *
     * @param neighbour_leaves the leaf records of S (the neighbours of the new vertex)
     */
    Subtree buildSubtree(const vector<int>& neighbour_leaves) {
        vector<int> visited_records;
        vector<int> active_records = neighbour_leaves;
        vector<int> next_active;
        for (int leaf : neighbour_leaves) {
            scratch[leaf].active_stamp = insertion_number;
        }
        bool root_reached = false;
        while ((!root_reached && active_records.size() >= 2) || (root_reached && active_records.size() >= 1)) {
            next_active.clear();
            for (int node : active_records) {
                climbOneLevel(node, visited_records, next_active, root_reached);
            }
            active_records.swap(next_active);
        }
        // the remaining active record(s) (the top of T') are visited without
        // climbing further
        for (int node : active_records) {
            markVisited(node, visited_records);
            if (node == root_record) root_reached = true;
            else scratch[node].parent_in_subtree = parentRecord(node);
        }
        // children lists of T'
        for (int node : visited_records) {
            if (node != root_record && scratch[scratch[node].parent_in_subtree].visited_stamp == insertion_number) {
                scratch[scratch[node].parent_in_subtree].children_in_subtree.push_back(node);
            }
        }
        // top of T' = the visited record whose parent is not visited
        Subtree subtree;
        for (int node : visited_records) {
            if (node == root_record || scratch[scratch[node].parent_in_subtree].visited_stamp != insertion_number) {
                subtree.top = node;
            }
        }
        // trim: if the top has a single child in T' and is not itself in S,
        // T(S) starts lower
        while (scratch[subtree.top].in_neighbourhood_stamp != insertion_number && scratch[subtree.top].children_in_subtree.size() == 1) {
            int child = scratch[subtree.top].children_in_subtree[0];
            scratch[subtree.top].visited_stamp = -1;
            subtree.top = child;
        }
        for (int node : visited_records)
            if (scratch[node].visited_stamp == insertion_number)
                subtree.records.push_back(node);
        return subtree;
    }

    /**
     * @brief Initialisation of Algorithm 5: the leaves of S below a node give
     *        a perfect marker in that node (the marker facing the leaf);
     *        non-leaf children are counted so that the peeling phases know
     *        when a node has no pending child.
     * @param
     */
    void markPerfectLeafMarkers(const Subtree& subtree) {
        for (int node : subtree.records) {
            if (node == subtree.top) continue;
            int parent = scratch[node].parent_in_subtree;
            if (records[node].kind == LEAF) {
                int marker = opposite[records[node].root_marker];
                perfect_stamp[marker] = insertion_number;
                scratch[parent].perfect_markers.push_back(marker);
            }
            else ++scratch[parent].node_child_count;
        }
        for (int node : subtree.records) {
            scratch[node].pending_child_count = scratch[node].node_child_count;
        }
    }

    /**
     * @brief Algorithm 5, step 1: peels pendant perfect subtrees bottom-up.
     *
     * A node (not the root) is a perfect pendant subtree if all its
     * T'-children are perfect (leaves of S or already-peeled subtrees) and the
     * neighbourhood of its root marker in G(u) is exactly those perfect
     * markers. Then the marker facing the node in its parent is perfect and
     * the node is "removed" (peeled) from T'.
     *
     * @param subtree  the subtree T' = T(S) after trimming
     * @param new_vertex  the vertex being inserted
     * @return true if the whole of T' turned out to be one perfect subtree:
     *         the insertion is then complete (cases 2 or 6, handled by
     *         edgeCase()); false if the insertion must go on.
     */
    bool peelPendantPerfectSubtrees(const Subtree& subtree, int new_vertex) {
        vector<int> ready_nodes;
        for (int node : subtree.records)
            if (node != root_record && records[node].kind != LEAF && scratch[node].pending_child_count == 0)
                ready_nodes.push_back(node);
        while (!ready_nodes.empty()) {
            int node = ready_nodes.back();
            ready_nodes.pop_back();
            int node_root_marker = records[node].root_marker;
            bool is_perfect = scratch[node].perfect_markers.size() == scratch[node].children_in_subtree.size() && labelDegree(node, node_root_marker) == (int)scratch[node].perfect_markers.size();
            if (is_perfect)
                for (int perfect_marker : scratch[node].perfect_markers)
                    if (!areAdjacent(node, node_root_marker, perfect_marker)) {
                        is_perfect = false;
                        break;
                    }
            if (!is_perfect) continue;
            scratch[node].peeled_stamp = insertion_number;
            int marker_facing_node = opposite[node_root_marker];
            perfect_stamp[marker_facing_node] = insertion_number;
            // the whole T' is one perfect subtree: the edge above it is perfect/empty (cases 2 or 6)
            if (node == subtree.top) {
                edgeCase(parentRecord(node), marker_facing_node, true, false, new_vertex);
                return true;
            }
            int parent_node = scratch[node].parent_in_subtree;
            scratch[parent_node].perfect_markers.push_back(marker_facing_node);
            ++scratch[parent_node].peeled_child_count;
            if (--scratch[parent_node].pending_child_count == 0 && parent_node != root_record)
                ready_nodes.push_back(parent_node);
        }
        return false;
    }

    /**
     * @brief Algorithm 5, step 2: peels the perfect part containing the root,
     *        top-down.
     *
     * While exactly one non-peeled node child below_node remains below `node`,
     * checks whether the neighbourhood of the marker towards below_node is
     * exactly the perfect markers of `node` (plus its root marker if perfect).
     * If so, below_node's root marker is perfect (state seen from above) and
     * we descend into it; otherwise `node` is the top of a fully-mixed part.
     *
     * @param subtree_top     top of T'
     * @param fully_mixed_out [out] true if a fully-mixed part was found (case 7)
     * @return the node where the descent stopped
     */
    int peelRootPerfectPart(int subtree_top, bool& fully_mixed_out) {
        int node = subtree_top;
        fully_mixed_out = false;
        while (true) {
            int unpeeled_child_count = scratch[node].node_child_count - scratch[node].peeled_child_count;
            if (unpeeled_child_count != 1) break;
            int below_node = -1;
            for (int child : scratch[node].children_in_subtree)
                if (records[child].kind != LEAF && scratch[child].peeled_stamp != insertion_number) {
                    below_node = child;
                    break;
                }
            int marker_towards_child = opposite[records[below_node].root_marker];
            bool root_side_matches = true;
            if (node != root_record) {
                int perfect_count = scratch[node].perfect_markers.size() + (scratch[node].root_marker_perfect_stamp == insertion_number ? 1 : 0);
                root_side_matches = labelDegree(node, marker_towards_child) == perfect_count;
                if (root_side_matches)
                    for (int perfect_marker : scratch[node].perfect_markers)
                        if (!areAdjacent(node, marker_towards_child, perfect_marker)) {
                            root_side_matches = false;
                            break;
                        }
                if (root_side_matches && scratch[node].root_marker_perfect_stamp == insertion_number && !areAdjacent(node, marker_towards_child, records[node].root_marker))
                    root_side_matches = false;
            }
            if (!root_side_matches) {
                fully_mixed_out = true;
                break;
            }
            scratch[below_node].root_marker_perfect_stamp = insertion_number;
            perfect_stamp[records[below_node].root_marker] = insertion_number;
            node = below_node;
        }
        return node;
    }

    /**
     * @brief Cases 1-6 (Theorem 4.14): exactly one node remains with all its
     *        tree-edges non-mixed.
     *
     * Collects the perfect markers of `node` (those of its children plus its
     * root marker if perfect) and dispatches to the prime or the degenerate
     * case.
     * @param node        the record with all its tree-edges non-mixed
     * @param new_vertex  the vertex being inserted
     */
    void insertAtUnmixedNode(int node, int new_vertex) {
        if (node == root_record) { // PP edge at the root leaf
            edgeCase(root_record, records[root_record].first_marker, true, true, new_vertex);
            return;
        }
        vector<int> all_perfect_markers = scratch[node].perfect_markers;
        if (scratch[node].root_marker_perfect_stamp == insertion_number)
            all_perfect_markers.push_back(records[node].root_marker);
        if (records[node].kind == PRIME)
            insertAtPrimeNode(node, all_perfect_markers, new_vertex);
        else
            insertAtDegenerateNode(node, all_perfect_markers, new_vertex);
    }

    /**
     * @brief Prime node (cases 3, 5, 6): LexBFS twin test, otherwise x joins
     *        the node.
     *
     * If the new marker would have the same neighbourhood in G(u) as the last
     * or the universal marker `candidate` (Section 5.2), then x is a twin of
     * candidate and the tree-edge at candidate is non-mixed (cases 5/6):
     * subdivide it. Otherwise x joins the prime node (case 3).
     * 
     * @param node        the prime node
     * @param all_perfect_markers  the perfect markers of `node`
     * @param new_vertex  the vertex being inserted
     */
    void insertAtPrimeNode(int node, const vector<int>& all_perfect_markers, int new_vertex) {
        int perfect_total = all_perfect_markers.size();
        for (int candidate : {records[node].last_inserted_marker, records[node].universal_marker}) {
            if (candidate < 0) continue;
            bool candidate_is_perfect = perfect_stamp[candidate] == insertion_number;
            int needed_degree = perfect_total - (candidate_is_perfect ? 1 : 0);
            if (labelDegree(node, candidate) != needed_degree) continue;
            bool neighbourhood_matches = true;
            for (int perfect_marker : all_perfect_markers)
                if (perfect_marker != candidate && !areAdjacent(node, candidate, perfect_marker)) {
                    neighbourhood_matches = false;
                    break;
                }
            if (neighbourhood_matches) { // twin: perfect opposite
                edgeCase(node, candidate, candidate_is_perfect, true, new_vertex); 
                return;
            }
        }
        attachPrime(node, all_perfect_markers, new_vertex);  // case 3
    }

    /**
     * @brief Clique or star node (cases 1-2, 4-6, Lemma 5.6).
     *
     * Finds from |P| a marker `selected_marker` whose tree-edge is non-mixed
     * and subdivides it; if there is none, the node is a degenerate hybrid
     * node (case 4): node-split into (P*, E*), then subdivide the new edge
     * (Prop. 4.16).
     * @param node        the degenerate node (clique or star)
     * @param all_perfect_markers  the perfect markers of `node`
     * @param new_vertex  the vertex being inserted
     */
    void insertAtDegenerateNode(int node, const vector<int>& all_perfect_markers, int new_vertex) {
        int perfect_total = all_perfect_markers.size();
        int node_marker_count = records[node].marker_count;
        int selected_marker = -1;
        bool selected_is_perfect = false;
        if (records[node].kind == CLIQUE) {
            if (perfect_total == node_marker_count) {
                selected_marker = records[node].first_marker;
                selected_is_perfect = true;
            }
            else if (perfect_total == node_marker_count - 1) {
                for (int marker = records[node].first_marker; marker >= 0; marker = next_marker[marker])
                    if (perfect_stamp[marker] != insertion_number) {
                        selected_marker = marker;
                        break;
                    }
            }
        } else {
            int centre = records[node].star_center;
            bool centre_is_perfect = perfect_stamp[centre] == insertion_number;
            if (perfect_total == node_marker_count) {
                selected_marker = centre;
                selected_is_perfect = true;
            }
            else if (perfect_total == node_marker_count - 1 && !centre_is_perfect)
                selected_marker = centre;
            else if (perfect_total == 1 && centre_is_perfect) {
                selected_marker = records[node].first_marker == centre ? next_marker[centre] : records[node].first_marker;
            }
            else if (perfect_total == 2 && centre_is_perfect) {
                for (int perfect_marker : all_perfect_markers)
                    if (perfect_marker != centre)
                        selected_marker = perfect_marker;
                selected_is_perfect = true;
            }
        }
        if (selected_marker >= 0) {
            edgeCase(node, selected_marker, selected_is_perfect, true, new_vertex);
            return;
        }
        // case 4: degenerate hybrid node
        vector<int> perfect_non_centre;
        bool centre_is_perfect_in_P = false;
        for (int perfect_marker : all_perfect_markers) {
            if (records[node].kind == STAR && perfect_marker == records[node].star_center)
                centre_is_perfect_in_P = true;
            else perfect_non_centre.push_back(perfect_marker);
        }
        if (perfect_non_centre.size() < 2 || node_marker_count - (int)perfect_non_centre.size() < 2)
            throw logic_error("degenerate hybrid node without split");
        int node_marker, new_node_marker;
        nodeSplit(node, perfect_non_centre, node_marker, new_node_marker);
        subdivide(node, node_marker, true, centre_is_perfect_in_P, new_vertex);
    }

    /**
     * @brief Case 7 (Prop. 4.20): the fully-mixed subtree rooted at `top`.
     *
     * Cleaning of the degenerate nodes (node-splits), contraction of the
     * cleaned subtree (node-joins) into one node, then attachPrime of the new
     * leaf to the perfect markers.
     * @param top         the top of the fully-mixed subtree
     * @param new_vertex  the vertex being inserted
     */
    void insertFullyMixed(int top, int new_vertex) {
        // the nodes of the fully-mixed subtree, in preorder from its top (peeled subtrees excluded)
        vector<int> fully_mixed_nodes;
        vector<int> node_stack{top};
        while (!node_stack.empty()) {
            int order_node = node_stack.back();
            node_stack.pop_back();
            fully_mixed_nodes.push_back(order_node);
            for (int child : scratch[order_node].children_in_subtree)
                if (records[child].kind != LEAF && scratch[child].peeled_stamp != insertion_number)
                    node_stack.push_back(child);
        }
        unordered_map<int, int> representative_of;  // T' node -> its node in the cleaned tree cl(ST(G))
        vector<int> final_perfect_markers;  // perfect markers of the contracted node (neighbours of x)
        cleanDegenerateNodes(fully_mixed_nodes, top, representative_of, final_perfect_markers);
        int contracted_node = contractSubtree(fully_mixed_nodes, representative_of, top);
        attachPrime(contracted_node, final_perfect_markers, new_vertex);  // x adjacent to the perfect markers of the contracted node
    }

    /**
     * @brief Cleaning (Def. 4.18) of the fully-mixed subtree.
     *
     * For every degenerate node of T', splits off (i) the perfect markers P*
     * and (ii) the empty markers E*, so that the node left in T' only keeps
     * the markers that are tree-edges of T' (Ml) plus stand-ins for the peeled
     * parts. Prime nodes are left untouched.
     *
     * @param fully_mixed_nodes      nodes of the subtree in preorder
     * @param top                    top of the subtree
     * @param representative_of      [out] node of T' -> the node that now represents it in the cleaned subtree
     * @param final_perfect_markers  [out] the perfect markers of the cleaned subtree (neighbours of the new vertex)
     */
    void cleanDegenerateNodes(const vector<int>& fully_mixed_nodes, int top, unordered_map<int, int>& representative_of, vector<int>& final_perfect_markers) {
        for (int order_node : fully_mixed_nodes) {
            vector<int> perfect_markers_here = scratch[order_node].perfect_markers;
            if (scratch[order_node].root_marker_perfect_stamp == insertion_number)
                perfect_markers_here.push_back(records[order_node].root_marker);
            Kind node_kind = records[order_node].kind;
            if (node_kind != CLIQUE && node_kind != STAR) { // prime: untouched
                representative_of[order_node] = order_node;
                final_perfect_markers.insert(final_perfect_markers.end(), perfect_markers_here.begin(), perfect_markers_here.end());
                continue;
            }
            int centre = records[order_node].star_center;
            bool centre_is_perfect_in_P = false;
            vector<int> perfect_non_centre, connecting_markers;
            for (int perfect_marker : perfect_markers_here) {
                if (node_kind == STAR && perfect_marker == centre)
                    centre_is_perfect_in_P = true;
                else perfect_non_centre.push_back(perfect_marker);
            }
            for (int child_node : scratch[order_node].children_in_subtree)
                if (records[child_node].kind != LEAF && scratch[child_node].peeled_stamp != insertion_number)
                    connecting_markers.push_back(opposite[records[child_node].root_marker]);
            if (order_node != top)
                connecting_markers.push_back(records[order_node].root_marker);
            int node_marker_count = records[order_node].marker_count;
            int empty_count = node_marker_count - (int)perfect_non_centre.size() - (int)connecting_markers.size();   // number of empty markers
            vector<int> kept_perfect_markers;  // perfect markers staying in the node of T'
            vector<int> markers_to_move;       // markers moved to the node of T'
            int node_marker, new_node_marker;
            if (perfect_non_centre.size() > 1 && node_marker_count - (int)perfect_non_centre.size() > 1) {
                nodeSplit(order_node, perfect_non_centre, node_marker, new_node_marker);
                perfect_stamp[node_marker] = insertion_number;  // perfect, stands for the peeled P*
                kept_perfect_markers.push_back(node_marker);
                markers_to_move.push_back(node_marker);
            } else {
                kept_perfect_markers = perfect_non_centre;
                markers_to_move = perfect_non_centre;
            }
            if (centre_is_perfect_in_P)
                kept_perfect_markers.push_back(centre);
            markers_to_move.insert(markers_to_move.end(), connecting_markers.begin(), connecting_markers.end());
            int representative = order_node;
            if (empty_count > 1) {
                int new_node = nodeSplit(order_node, markers_to_move, node_marker, new_node_marker);               // new_node = node of T'; order_node keeps E* (peeled)
                representative = new_node;
                // a perfect centre stays with E*; its stand-in in the new node is perfect
                if (centre_is_perfect_in_P) {
                    kept_perfect_markers.erase(std::find(kept_perfect_markers.begin(), kept_perfect_markers.end(), centre));
                    kept_perfect_markers.push_back(new_node_marker);
                    perfect_stamp[new_node_marker] = insertion_number;
                }
            }
            representative_of[order_node] = representative;
            final_perfect_markers.insert(final_perfect_markers.end(), kept_perfect_markers.begin(), kept_perfect_markers.end());
        }
    }

    /**
     * @brief Contraction (Algorithm 9): joins all nodes of the cleaned
     *        fully-mixed subtree into one node.
     *
     * Done in three phases (the joins commute, Remark 2.14; the order only
     * matters for the cost, Lemma 6.16). Nodes are visited bottom-up (reverse
     * preorder), so a node's children are handled before it:
     *  1. a star whose root marker is its centre absorbs each of its children
     *     in the subtree;
     *  2. a node whose root marker has degree one is absorbed by its parent;
     *  3. everything left is joined into its parent.
     *
     * @param fully_mixed_nodes nodes of the subtree before cleaning, in preorder
     * @param representative_of node of T' -> its node in the cleaned subtree
     * @param top               top of the subtree before cleaning
     * @return the single node the subtree has been contracted into (the representative of `top`)
     */
    int contractSubtree(const vector<int>& fully_mixed_nodes, const unordered_map<int, int>& representative_of, int top) {
        vector<int> contraction_nodes;  // the nodes of the cleaned subtree, in preorder
        for (int order_node : fully_mixed_nodes) {
            int representative = representative_of.at(order_node);
            contraction_nodes.push_back(representative);
            scratch[representative].in_contraction_stamp = insertion_number;
        }
        const int contracted_node = representative_of.at(top);
        auto contractionParent = [&](int node) {  // parent of node if it is in the contracted subtree, else -1
            if (node == contracted_node) return -1;
            int parent = parentRecord(node);
            return scratch[parent].in_contraction_stamp == insertion_number ? parent : -1;
        };
        // phase 1: a star whose root marker is its centre absorbs each of its children in the subtree
        for (int index = (int)contraction_nodes.size() - 1; index >= 0; --index) {
            int node = contraction_nodes[index];
            if (!records[node].is_alive || records[node].kind != STAR || records[node].star_center != records[node].root_marker)
                continue;
            vector<int> subtree_children;
            for (int marker = records[node].first_marker; marker >= 0; marker = next_marker[marker]) {
                if (marker == records[node].root_marker) continue;
                int child = record_rooted_at_marker[opposite[marker]];
                if (records[child].kind != LEAF && scratch[child].in_contraction_stamp == insertion_number)
                    subtree_children.push_back(child);
            }
            for (int child : subtree_children)
                nodeJoin(node, child);
        }
        // phase 2: a node whose root marker has degree one is absorbed by its parent
        for (int index = (int)contraction_nodes.size() - 1; index >= 0; --index) {
            int node = contraction_nodes[index];
            if (!records[node].is_alive || labelDegree(node, records[node].root_marker) != 1)   
                continue;
            int parent = contractionParent(node);
            if (parent >= 0)
                nodeJoin(parent, node);
        }
        // phase 3: join everything that is left
        for (int index = (int)contraction_nodes.size() - 1; index >= 0; --index) {
            int node = contraction_nodes[index];
            if (!records[node].is_alive)   
                continue;
            int parent = contractionParent(node);
            if (parent >= 0)
                nodeJoin(parent, node);
        }
        return contracted_node;
    }

    /**
     * @brief Inserts vertex new_vertex, adjacent to the already-inserted
     *        vertices neighbour_vertices, into the split tree.
     *
     * Preconditions: at least two vertices are already inserted; the new
     * vertex has at least one neighbour in the current graph (guaranteed by
     * LexBFS order on a connected graph); the insertion order is an LexBFS
     * order.
     *
     * Steps:
     *  0. |S| = 1: the new vertex is pendant, handled directly.
     *  1. buildSubtree(): Algorithm 4, computes T' = T(S).
     *  2. markPerfectLeafMarkers(), peelPendantPerfectSubtrees(),
     *     peelRootPerfectPart():
     *     Algorithm 5, determines the states of the markers of T'. Afterwards
     *     either a unique node/edge identifies one of cases 1-6, or a
     *     fully-mixed subtree remains (case 7).
     *  3. Cases 1-6: insertAtUnmixedNode() -> edgeCase / attachPrime /
     *     node-split + subdivide.
     *  4. Case 7: insertFullyMixed() -> cleanDegenerateNodes(),
     *     contractSubtree(), attachPrime().
     * Only the markers of T' are touched.
     *
     * Throws std::logic_error if an internal invariant fails (should not
     * happen for LexBFS orders).
     */
    void insert(int new_vertex, const vector<int>& neighbour_vertices) {
        ++insertion_number;
        vector<int> neighbour_leaves;
        for (int vertex : neighbour_vertices)
            neighbour_leaves.push_back(leaf_record_of_vertex[vertex]);
        for (int leaf : neighbour_leaves)
            scratch[leaf].in_neighbourhood_stamp = insertion_number;

        // ---- step 0: single neighbour => pendant vertex (the edge to that leaf is perfect/empty) ----
        if (neighbour_leaves.size() == 1) {
            int leaf = neighbour_leaves[0];
            int leaf_marker = records[leaf].first_marker;
            int marker = opposite[leaf_marker];
            int node = (leaf == root_record) ? record_rooted_at_marker[marker] : parentRecord(leaf);
            edgeCase(node, marker, true, false, new_vertex);
            return;
        }
        // ---- steps 1-2: T(S) and the states of its markers ----
        Subtree subtree = buildSubtree(neighbour_leaves);
        markPerfectLeafMarkers(subtree);
        if (peelPendantPerfectSubtrees(subtree, new_vertex)) return;
        bool fully_mixed;
        int node = peelRootPerfectPart(subtree.top, fully_mixed);
        // ---- steps 3-4: cases 1-6, or case 7 ----
        if (!fully_mixed && scratch[node].node_child_count - scratch[node].peeled_child_count == 0) 
            insertAtUnmixedNode(node, new_vertex);
        else
            insertFullyMixed(node, new_vertex);
    }


    // ====================================================================================
    //  Verification (testing only; NOT part of the algorithm, and quadratic or worse)
    // ====================================================================================

    /**
     * @brief Result of verify(): ok == true, or a short reason in `why`.
     */
    struct Report {
        bool ok = true;
        string reason;
    };

    /**
     * @brief Checks the whole structure against `graph` (adjacency lists on
     *        vertices 0..n-1).
     *
     * Checks, in order: consistency of opposite / record_rooted_at_marker /
     * union-find parents; accessibility graph equals graph; every node has
     * >= 3 markers; prime-kind nodes have no split (brute force over subsets,
     * skipped when a prime node has more than 14 markers); the tree is reduced
     * (no clique-clique edge, no star-join).
     * @param graph  the adjacency lists of the graph being represented
     * @return a Report with ok == true if the structure is valid, or false
     *         and a short reason otherwise
     */
    Report verify(const vector<vector<int>>& graph) {
        Report report;
        auto fail = [&](string reason_text) {
            report.ok = false;
            report.reason = reason_text;
            return report;
        };
        // traverse from the root, structural consistency
        vector<int> record_of_marker(opposite.size(), -1);
        vector<int> records_to_visit{root_record};
        vector<int> reachable_records;
        while (!records_to_visit.empty()) {
            int record_id = records_to_visit.back();
            records_to_visit.pop_back();
            reachable_records.push_back(record_id);
            if (!records[record_id].is_alive)
                return fail("dead record reachable");
            for (int marker = records[record_id].first_marker; marker >= 0; marker = next_marker[marker]) {
                record_of_marker[marker] = record_id;
                if (opposite[opposite[marker]] != marker)
                    return fail("opposite not symmetric");
                if (marker == records[record_id].root_marker)
                    continue;
                int child = record_rooted_at_marker[opposite[marker]];
                if (child < 0 || records[child].root_marker != opposite[marker])
                    return fail("bad record_rooted_at_marker");
                if (parentRecord(child) != record_id)
                    return fail("union-find parent mismatch");
                records_to_visit.push_back(child);
            }
        }

        // accessibility graph = input graph
        set<pair<int,int>> accessible_pairs;
        set<pair<int,int>> graph_pairs;
        for (size_t source_vertex = 0; source_vertex < graph.size(); ++source_vertex)
            for (int neighbour_vertex : graph[source_vertex])
                if ((int)source_vertex < neighbour_vertex)
                    graph_pairs.insert({(int)source_vertex, neighbour_vertex});
        for (int record_id : reachable_records)
            if (records[record_id].kind == LEAF) {
                int leaf_vertex = leaf_vertex_of_marker[records[record_id].first_marker];
                vector<int> markers_to_follow{opposite[records[record_id].first_marker]};
                while (!markers_to_follow.empty()) {
                    int arrival_marker = markers_to_follow.back();
                    markers_to_follow.pop_back();
                    int arrival_record = record_of_marker[arrival_marker];
                    for (int label_neighbour = records[arrival_record].first_marker; label_neighbour >= 0; label_neighbour = next_marker[label_neighbour])
                        if (areAdjacent(arrival_record, arrival_marker, label_neighbour)) {
                            int reached_marker = opposite[label_neighbour];
                            if (records[record_of_marker[reached_marker]].kind == LEAF) {
                                int target_vertex = leaf_vertex_of_marker[reached_marker];
                                accessible_pairs.insert({min(leaf_vertex,target_vertex), max(leaf_vertex,target_vertex)});
                            }
                           else markers_to_follow.push_back(reached_marker);
                        }
                }
            }
        if (accessible_pairs != graph_pairs)
            return fail("accessibility graph differs");

        // node sizes and primality (a subset A of markers is a split iff
        // |A|,|rest| >= 2 and every marker of A with an outside neighbour is
        // adjacent to every marker of the rest with an inside neighbour)
        for (int record_id : reachable_records) {
            if (records[record_id].kind == LEAF)
                continue;
            int node_marker_count = records[record_id].marker_count;
            if (node_marker_count < 3)
                return fail("node with < 3 markers");
            if (records[record_id].kind == PRIME) {
                vector<int> node_markers = markersOf(record_id);
                if (node_marker_count > 14)
                    continue;
                for (int mask = 1; mask < (1 << node_marker_count); ++mask) {
                    int subset_size = popcount32(mask);
                    if (subset_size < 2 || node_marker_count - subset_size < 2)
                        continue;
                    vector<char> has_outside_neighbour(node_marker_count, 0);
                    for (int i = 0; i < node_marker_count; ++i)
                        for (int j = 0; j < node_marker_count; ++j)
                            if (((mask >> i) & 1) != ((mask >> j) & 1) && areAdjacent(record_id, node_markers[i], node_markers[j]))
                                has_outside_neighbour[i] = 1;
                    bool is_split = true;
                    for (int i = 0; i < node_marker_count && is_split; ++i)
                        for (int j = 0; j < node_marker_count && is_split; ++j)
                            if (((mask >> i) & 1) && !((mask >> j) & 1))
                                is_split = (areAdjacent(record_id, node_markers[i], node_markers[j]) == (has_outside_neighbour[i] && has_outside_neighbour[j]));
                    if (is_split)
                        return fail("prime node has a split");
                }
            }
        }
        // reducedness: no clique-join and no star-join possible on any internal tree-edge
        for (int record_id : reachable_records)
            for (int marker = records[record_id].first_marker; marker >= 0; marker = next_marker[marker]) {
                int opposite_marker = opposite[marker];
                if (marker > opposite_marker)
                    continue;
                int record_a = record_of_marker[marker];
                int record_b = record_of_marker[opposite_marker];
                if (records[record_a].kind == LEAF || records[record_b].kind == LEAF)
                    continue;
                if (records[record_a].kind == CLIQUE && records[record_b].kind == CLIQUE)
                    return fail("adjacent cliques");
                if (records[record_a].kind == STAR && records[record_b].kind == STAR && ((records[record_a].star_center == marker) != (records[record_b].star_center == opposite_marker)))
                    return fail("star-join possible");
            }
        return report;
    }
};

// ----- LexBFS (linear, partition refinement)

/**
 * @brief Linear-time lexicographic breadth-first search by partition
 *        refinement, array-based.
 *
 * Based on Habib, McConnell, Paul, Viennot (TCS 234, 2000) and
 * Tedder, Corneil, Habib, Paul (ICALP 2008, doi:10.1007/978-3-540-70575-8_52).
 *
 * The order under construction is an array `order`; the unvisited vertices
 * form consecutive "parts" of that array (vertices with the same label), each
 * part being a contiguous range of positions that starts at
 * `part_first_position[part]`. Visiting vertex v moves each unvisited
 * neighbour to the front of its part and cuts the moved prefix off as a new
 * part, so that neighbours of v come before non-neighbours inside every part.
 * Total time O(n + m).
 *
 * Part identifiers are never recycled and each (vertex, later neighbour) pair
 * creates at most one, hence at most 1 + #edges identifiers are used (the
 * arrays are sized accordingly).
 *
 * @param graph        simple undirected graph (no loops, no multi-edges) as adjacency lists on 0..n-1
 * @param start_vertex first vertex of the order
 * @return the vertices in LexBFS order
 */
vector<int> LexBFS(const vector<vector<int>>& graph, int start_vertex = 0) {
    const int vertex_count = (int)graph.size();
    vector<int> order(vertex_count);  // order[i] = vertex at position i
    vector<int> position(vertex_count);  // position = inverse
    iota(order.begin(), order.end(), 0);
    iota(position.begin(), position.end(), 0);
    if (start_vertex > 0) {                                           // put the start vertex first
        order[0] = start_vertex;
        order[start_vertex] = 0;
        position[0] = start_vertex;
        position[start_vertex] = 0;
    }
    size_t edge_count = 0;
    for (const auto& adjacent : graph)
        edge_count += adjacent.size();
    edge_count /= 2;
    const size_t max_parts = max<size_t>(vertex_count, edge_count) + 2;
    vector<int> part_of_position(vertex_count, 0);  // part containing the vertex at each position
    vector<int> part_first_position(max_parts, 0);  // first position of each part (not yet visited)
    vector<int> split_part(max_parts, 0);  // part created from this one during the current step (0 = none)
    int part_count = 1;  // initially one part holding every vertex

    for (int step = 0; step < vertex_count; ++step) {
        const int part_count_before_step = part_count;
        part_first_position[part_of_position[step]] += 1;   // the vertex at `step` leaves its part (it is visited now)
        const int visited_vertex = order[step];
        for (int neighbour : graph[visited_vertex]) {
            int neighbour_position = position[neighbour];
            if (neighbour_position <= step) continue;                 // already visited
            const int part = part_of_position[neighbour_position];
            const int first_position = part_first_position[part];
            // if the part has at least two vertices, swap the neighbour to the front of the part
            // and shrink the part from the front
            if (first_position < vertex_count - 1 && part_of_position[first_position + 1] == part) {
                if (first_position != neighbour_position) {
                    const int displaced = order[first_position];
                    position[displaced] = neighbour_position;
                    position[neighbour] = first_position;
                    order[neighbour_position] = displaced;
                    order[first_position] = neighbour;
                    neighbour_position = first_position;
                }
                part_first_position[part] += 1;
            }
            // first neighbour met in this part during this step: open the new part holding the neighbours
            if (split_part[part] < part_count_before_step) {
                split_part[part] = part_count;
                part_first_position[part_count] = neighbour_position;
                split_part[part_count] = 0;
                ++part_count;
            }
            part_of_position[neighbour_position] = split_part[part];
        }
    }
    return order;
}

// ---------------------------------------------------------------------- driver (Algorithm 1)

/**
 * @brief Builds the split tree of a connected graph by inserting its vertices
 *        in the given order.
 *
 * @param graph connected simple graph, adjacency lists on 0..n-1
 * @param order MUST be a LexBFS order of graph (e.g. from LexBFS(), or
 *              exported from Sage's slice decomposition); the insertion
 *              algorithm relies on it (twin test). Not checked here.
 * @return the split tree (use FastGLT::verify to check it against graph)
 */
FastGLT splitTreeFromOrder(const vector<vector<int>>& graph, const vector<int>& order) {
    FastGLT tree;
    vector<char> is_inserted(graph.size(), 0);  // in[v]: v already inserted
    for (size_t i = 0; i < order.size(); ++i) {
        int vertex = order[i];
        if (i == 0) tree.addFirst(vertex);
        else if (i == 1) tree.addSecond(vertex);
        else {
            vector<int> earlier_neighbours;
            for (int neighbour : graph[vertex])
                if (is_inserted[neighbour])
                    earlier_neighbours.push_back(neighbour);
            tree.insert(vertex, earlier_neighbours);
        }
        is_inserted[vertex] = 1;
    }
    return tree;
}

/**
 * @brief Builds the split tree of a connected graph by inserting its vertices
 *        in LexBFS order.
 * @param graph connected simple graph, adjacency lists on 0..n-1
 * @return the split tree (use FastGLT::verify to check it against graph)
 */
FastGLT splitTree(const vector<vector<int>>& graph) {
    return splitTreeFromOrder(graph, LexBFS(graph));
}



// =====================================================================================
// Exporting the split tree to a compact format suitable for external use
// =====================================================================================

/**
 * @brief Exporting the split tree to a compact format suitable for external use
 *        (e.g. Python, Sage, etc.)
 */
struct GltExport {
    int root_vertex = -1;  // the leaf the tree is rooted at (first vertex of the LexBFS order)
    vector<int> order; // the LexBFS order that was used
    // one entry per internal node j:
    vector<int> node_kind;  // 0 = PRIME, 1 = CLIQUE, 2 = STAR
    vector<int> node_parent;  // code of the tree-neighbour above the node (a leaf for the top node)
    vector<int> node_star_centre;  // STAR: code of the tree-neighbour facing the centre marker; otherwise 0 (unused)
    vector<int> node_marker_start;  // CSR offsets into node_marker_neighbour (size = #nodes + 1)
    vector<int> node_marker_neighbour;  // codes of the tree-neighbours of the node, i.e. its markers
    vector<int> node_edge_start;  // CSR offsets into node_edge_pairs, counted in pairs (size = #nodes + 1)
    vector<int> node_edge_pairs;  // PRIME nodes: label edges as consecutive pairs of codes (marker names)
};


/**
 * @brief Is the graph connected? (iterative DFS from vertex 0)
 * 
 * @param graph simple undirected graph, adjacency lists on 0..n-1
 * @return true if the graph is connected, false otherwise
 */
static bool isConnected(const vector<vector<int>>& graph) {
    vector<char> seen(graph.size(), 0);
    vector<int> stack{0};
    seen[0] = 1;
    int seen_count = 1;
    while (!stack.empty()) {
        int vertex = stack.back();
        stack.pop_back();
        for (int neighbour : graph[vertex])
            if (!seen[neighbour]) {
                seen[neighbour] = 1;
                ++seen_count;
                stack.push_back(neighbour);
            }
    }
    return seen_count == (int)graph.size();
}


/**
 * @brief Builds the split tree of a connected simple graph and exports it
 *        to GltExport.
 * 
 * @param graph connected simple graph, adjacency lists on 0..n-1, n >= 3
 * @param given_order empty => compute LexBFS from start_vertex; otherwise it
 *                    must be a LexBFS order (the function does not check it)
 * @param start_vertex first vertex of the order if given_order is empty
 * @param run_verify whether to verify the tree against the graph (slower,
 *                   but useful for testing). Note that the verification skips
 *                   the primality check for prime nodes with more than
 *                   14 markers.
 * @param out [out] the exported split tree
 * @param error [out] error message on failure
 * @return true on success, false on failure (error message filled)
 */
inline bool glt_build(const vector<vector<int>>& graph,
                      const vector<int>& given_order,
                      int start_vertex,
                      bool run_verify,
                      GltExport& out,
                      string& error) {
    try {
        int n = (int)graph.size();
        if (n < 3) {
            error = "graph must have at least 3 vertices";
            return false;
        }
        if (!isConnected(graph)) {
            error = "graph must be connected";
            return false;
        }
        if (!given_order.empty() && (int)given_order.size() != n) {
            error = "given order has the wrong length";
            return false;
        }
        // set the order: either the given one, or compute LexBFS from start_vertex
        out.order = given_order.empty() ? LexBFS(graph, start_vertex) : given_order;
        // build the split tree from the order
        FastGLT tree = splitTreeFromOrder(graph, out.order);

        if (run_verify) {
            auto report = tree.verify(graph);
            if (!report.ok) {
                error = "verification failed: " + report.reason;
                return false;
            }
        }
        out.root_vertex = out.order[0];

        // pass 1: number the internal nodes in preorder
        const int record_count = (int)tree.records.size();
        vector<int> node_number(record_count, -1);
        vector<int> internal_records;
        vector<int> stack_{tree.root_record};
        while (!stack_.empty()) {
            int record_id = stack_.back();
            stack_.pop_back();
            if (tree.records[record_id].kind != FastGLT::LEAF) {
                node_number[record_id] = (int)internal_records.size();
                internal_records.push_back(record_id);
            }
            vector<int> children;
            for (int m = tree.records[record_id].first_marker; m >= 0; m = tree.next_marker[m])
                if (m != tree.records[record_id].root_marker)
                    children.push_back(tree.record_rooted_at_marker[tree.opposite[m]]);
            for (int i = (int)children.size() - 1; i >= 0; --i)
                stack_.push_back(children[i]);
        }

        // helper: code of a tree-neighbour (leaf or internal node) of a record
        auto code_of = [&](int record_id) {
            const auto& r = tree.records[record_id];
            return r.kind == FastGLT::LEAF ? tree.leaf_vertex_of_marker[r.first_marker] : -(node_number[record_id] + 1);
        };

        // pass 2: export each internal node
        out.node_marker_start.push_back(0); out.node_edge_start.push_back(0);
        for (int record_id : internal_records) {
            const auto& r = tree.records[record_id];
            out.node_kind.push_back(r.kind == FastGLT::PRIME ? 0 : r.kind == FastGLT::CLIQUE ? 1 : 2);
            out.node_parent.push_back(code_of(tree.parentRecord(record_id)));
            vector<int> markers = tree.markersOf(record_id);
            // helper: code of the tree-neighbour of record_id facing marker m
            auto neighbour_code = [&](int m) {
                return code_of(m == r.root_marker ? tree.parentRecord(record_id) : tree.record_rooted_at_marker[tree.opposite[m]]);
            };
            out.node_star_centre.push_back(r.kind == FastGLT::STAR ? neighbour_code(r.star_center) : 0);
            for (int m : markers) {
                out.node_marker_neighbour.push_back(neighbour_code(m));
            }
            out.node_marker_start.push_back((int)out.node_marker_neighbour.size());
            if (r.kind == FastGLT::PRIME) {
                for (int a : markers) {
                    for (int b : tree.label_adjacency[a]) {
                        if (a < b) {
                            out.node_edge_pairs.push_back(neighbour_code(a));
                            out.node_edge_pairs.push_back(neighbour_code(b));
                        }
                    }
                }
            }
            out.node_edge_start.push_back((int)out.node_edge_pairs.size() / 2);
        }
        return true;
    } catch (const exception& e) {
        error = e.what();
        return false;
    }
}
