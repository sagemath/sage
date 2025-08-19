/*
 *      Copyright (C) 2024 Cyril Bouvier <cyril.bouvier@lirmm.fr>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *                 https://www.gnu.org/licenses/
 */

/*
 * This file contains inline implementations of utility structs, classes and
 * functions used in the implementation of the modular decomposition method
 *
 * AUTHORS:
 *
 *  - Cyril Bouvier (2024): code for second implementation of the linear time
 *  algorithm of D. Corneil, M. Habib, C. Paul and M. Tedder [TCHP2008]_
 */
#include <algorithm>
#include <iostream>
#include <deque>
#include <list>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

/* Labels attached to nodes of a partitive forest. For the meaning of the
 * different values, see algorithms 3 and 4 of [TCHP2008]_.
 */
enum class Label : uint8_t {
    EMPTY = 0b00,
    HOMOGENEOUS = 0b01,
    BROKEN = 0b10,
    DEAD = 0b11
};

/*
 * Flags attached to nodes of a partitive forest. For the meaning of the flag,
 * see algorithms 3 and 4 of [TCHP2008]_.
 */
enum class Flag : uint8_t {
    UNFLAGGED = 0b00,
    FLAGGED = 0b01
};

/*
 * A node of a modular decomposition tree is either a leaf or an internal node.
 * An internal node can be prime, series or parallel.
 */
enum class Type : uint8_t {
    PRIME = 0,
    SERIES = 1,
    PARALLEL = 2,
    LEAF = 3
};

/*
 * This struct is used to have a nice interface to the data describing a slice
 * decomposition. The slice decomposition is needed by the algorithm that
 * compute the modular decomposition.
 */
struct SDData {
    /*
     * Set the members of the struct; it is useful to initialized the struct
     * from the output of the extended_lex_BFS method.
     */
    void set_from_data(size_t lex_label_offset_arg,
                       const int* sigma_arg,
                       const size_t *xslice_len_arg,
                       const std::vector<int> *lex_label_arg) {
        lex_label_offset = lex_label_offset_arg;
        sigma = sigma_arg;
        xslice_len = xslice_len_arg;
        lex_label = lex_label_arg;
    }

    /*
     * Set the members of the struct to represent the subslice of sd starting at
     * the given offset.
     */
    void set_to_subslice(const SDData &sd, size_t offset) {
        lex_label_offset = sd.lex_label[offset].size();
        sigma = sd.sigma + offset;
        xslice_len = sd.xslice_len + offset;
        lex_label = sd.lex_label + offset;
    }

    /*
     * Return the number of lexicographic labels of the ith vertex of the slice
     * decomposition.
     */
    size_t lex_label_size(size_t i) const {
        if (lex_label[i].size() <= lex_label_offset) {
            return 0;
        } else {
            return lex_label[i].size() - lex_label_offset;
        }
    }

    /*
     * Return the pointer to the lexicographic labels of the ith vertex of the
     * slice decomposition.
     */
    const int * lex_label_ptr(size_t i) const {
        if (lex_label[i].size() <= lex_label_offset) {
            return nullptr;
        } else {
            return lex_label[i].data() + lex_label_offset;
        }
    }

    /* Return the index of the first slice (always 1). */
    size_t first_slice_index() const {
        return 1;
    }

    /* Return the index of the slice following the one starting at idx. */
    size_t next_slice_index(size_t idx) const {
        return idx + xslice_len[idx];
    }

    /* Return the number of vertices in the slice decomposition. */
    size_t size() const {
        return xslice_len[0];
    }

    /*
     * Check whether the pivot has any neighbor (i.e., the first slice has no
     * lexicographic labels).
     */
    bool is_pivot_isolated () const {
        return lex_label[1].size() <= lex_label_offset;
    }

    size_t lex_label_offset;
    const int *sigma;
    const size_t *xslice_len;
    const std::vector<int> *lex_label;
};

/*
 * This struct represents a node of a modular decomposition tree. It contains a
 * pointer to its parent (or nullptr for the root), a list of children (must be
 * empty for leaf) and a type (from the enum Type).
 * A leaf of a modular decomposition tree corresponds to a vertex of the graph.
 * The id of the vertex is stored in the vertex attribute of the struct. For
 * internal nodes, the vertex attribute is used to store the id of any vertex
 * belonging to the corresponding module.
 * The attributes label, flag, slice, cc_tag are used by the different parts of
 * the algorithm that computes the modular decomposition tree.
 */
struct md_tree_node {
    /* ctor for non-leaf node */
    md_tree_node(Type type, Label label, Flag flag)
        : parent(nullptr), vertex(INT_MAX), type(type),
          label(label), flag(flag),
          slice(SIZE_MAX), cc_tag(SIZE_MAX) {
    }

    /* ctor for leaf */
    md_tree_node(int vertex)
        : parent(nullptr), vertex(vertex), type(Type::LEAF),
          label(Label::EMPTY), flag(Flag::UNFLAGGED),
          slice(SIZE_MAX), cc_tag(SIZE_MAX) {
    }

    /* ctor for non-leaf node, with default label and flag. */
    md_tree_node(Type type)
        : md_tree_node(type, Label::EMPTY, Flag::UNFLAGGED) {
    }

    /* check whether the node is a leaf. */
    bool is_leaf() const {
        return type == Type::LEAF;
    }

    /* check whether the node is prime. */
    bool is_prime() const {
        return type == Type::PRIME;
    }

    /* check whether the node is series. */
    bool is_series() const {
        return type == Type::SERIES;
    }

    /* check whether the node is parallel. */
    bool is_parallel() const {
        return type == Type::PARALLEL;
    }

    /* check whether the node is degenerate (i.e., series or parallel). */
    bool is_degenerate() const {
        return type == Type::SERIES || type == Type::PARALLEL;
    }

    /*
     * The following methods are used to check the different possible status of
     * the label of the node
     */
    bool is_empty() const {
        return label == Label::EMPTY;
    }

    bool is_homogeneous() const {
        return label == Label::HOMOGENEOUS;
    }

    bool is_homogeneous_or_empty() const {
        return !(static_cast<uint8_t>(label) >> 1U);
    }

    bool is_broken() const {
        return label == Label::BROKEN;
    }

    bool is_dead() const {
        return label == Label::DEAD;
    }

    bool is_dead_or_broken() const {
        return static_cast<uint8_t>(label) >> 1U;
    }

    /* add the node c at the beginning of the list of children. */
    void prepend_new_child(md_tree_node *c) {
        c->parent = this;
        if (children.empty()) {
            vertex = c->vertex;
        }
        children.push_front(c);
    }

    /* add the node c at the end of the list of children. */
    void append_new_child(md_tree_node *c) {
        c->parent = this;
        if (children.empty()) {
            vertex = c->vertex;
        }
        children.push_back(c);
    }

    /*
     * "Stole" the children from the node n and add them at the end of the list
     * of children. ("stole" here means that at the end of this method, the list
     * of children of n will be empty).
     */
    void append_stolen_children_from(md_tree_node *n) {
        if (!n->children.empty()) {
            for (md_tree_node *c: n->children) {
                c->parent = this;
            }
            if (children.empty()) {
                vertex = n->children.front()->vertex;
            }
            children.splice(children.end(), n->children);
        }
    }

    /* Set the label and flag for all nodes of the tree. */
    void set_label_and_flag_recursively(Label l, Flag f) {
        label = l;
        flag = f;
        for (md_tree_node *c: children) {
            c->set_label_and_flag_recursively(l, f);
        }
    }

    md_tree_node *parent;
    std::list<md_tree_node *> children;
    int vertex;
    Type type;
    Label label;
    Flag flag;
    size_t slice;
    size_t cc_tag;
};

/* A forest is a list of tree */
using md_forest = std::list<md_tree_node *>;

/*
 * This struct is used to gather data structures needed by the different parts
 * of the algorithm computing the modular decomposition tree. It is created
 * at the beginning of the algorithm and pass along to any function that needs
 * it. The hope is that it will reduce the number of allocations/deallocations
 * because the number of creations and destructions of objects will be reduced.
 */
struct ScratchData {
    /*
     * Sub structure containing objects needed by the implementation of
     * algorithms 3 and 4 of [TCHP2008]_.
     */
    struct MDSequences {
        std::unordered_map<int, md_tree_node *> leaves;
        std::unordered_set<md_tree_node *> Marked;
        std::unordered_set<md_tree_node *> Full;
        std::deque<md_tree_node *> Explore;
    };

    /*
     * Sub structure containing objects needed by the implementation of
     * algorithms 5 and 6 of [TCHP2008]_.
     */
    struct Clusters {
        /* Assumes i < p < j where p is the index of the "cluster" {x}. */
        bool are_clusters_non_adjacent(size_t i, size_t j) const {
            size_t sj = clusters[j].front()->slice;
            auto e = std::make_pair(0, sj);
            for (const md_tree_node *mi: clusters[i]) {
                e.first = mi->vertex;
                auto it = module_slice_adjacency.find(e);
                if (it != module_slice_adjacency.end()) {
                    return false;
                }
            }
            return true;
        }

        struct pair_hash {
            inline size_t operator()(const std::pair<int, size_t> &v) const {
                return ((size_t) v.first)*31+v.second;
            }
        };

        std::vector<std::vector<md_tree_node *>> clusters;
        /* adjacency list between a module (represented using the corresponding
         * .vertex from the root node) and a slice.
         */
        std::unordered_set<std::pair<int, size_t>, pair_hash> module_slice_adjacency;
        std::vector<size_t> Left;
        std::vector<size_t> Right;
        std::unordered_map<int, size_t> cluster_of_v;
    };

    /* Allocate the node corresponding to the vertex, and add it to the map */
    md_tree_node *new_leaf(int vertex) {
        return mdseq.leaves[vertex] = new md_tree_node(vertex);
    }

    MDSequences mdseq;
    Clusters clusters;
};

/*
 * This function deallocate (using delete) the node n and all of its
 * descendants. It is needed to deallocate the modular decomposition tree
 * computed by the function corneil_habib_paul_tedder_inner.
 */
void dealloc_md_tree_nodes_recursively(md_tree_node *n) {
    for (md_tree_node *c: n->children) {
        dealloc_md_tree_nodes_recursively(c);
    }
    delete n;
}

/* Preprocess the trees in the forest: set the label to empty and the flag to
 * unflagged for all the nodes, and set the cc_tag needed to compute the cluster
 * later.
 */
void md_forest_preprocess(md_forest &MDi) {
    Type one_cc_type = Type::PARALLEL; /* only for first iteration */
    size_t s = 0;
    for (md_tree_node *md: MDi) {
        md->set_label_and_flag_recursively(Label::EMPTY, Flag::UNFLAGGED);
        md->slice = s;
        if (md->type == Type::PRIME || md->type == one_cc_type) {
            md->cc_tag = 0;
        } else {
            md->cc_tag = SIZE_MAX;
            size_t i = 0;
            for (md_tree_node *c: md->children) {
                c->cc_tag = i;
                i++;
            }
        }
        one_cc_type = Type::SERIES;
        s += 1;
    }
}

/*
 * This function set to BROKEN the ancestors of DEAD nodes and gather into one
 * node the HOMOGENEOUS and EMPTY childrend of a broken and degenerate node.
 * Corresponds to the end of algorithm 3 of [TCHP2008]_.
 */
void mark_partitive_forest_finish_inner_rec(md_tree_node *r) {
    size_t nb = 0; /* number of HOMOGENEOUS or EMPTY children */

    /* Do a postorder visit: so we first visit the children */
    for (md_tree_node *c: r->children) {
        mark_partitive_forest_finish_inner_rec(c);
        nb += c->is_homogeneous_or_empty();
    }

    if (r->is_dead_or_broken()) {
        if (r->parent != nullptr && !(r->parent->is_dead())) {
            /* if parent.label is not DEAD set it to BROKEN */
            r->parent->label = Label::BROKEN;
        }
        if (r->is_broken() && r->is_degenerate() && nb > 1) {
            md_tree_node *newnode = new md_tree_node(r->type, Label::EMPTY,
                                                              Flag::UNFLAGGED);

            /* Iterate over the children to gather HOMOGENEOUS and EMPTY
             * child under newnode
             * */
            auto it = r->children.begin();
            while (it != r->children.end()) {
                if ((*it)->is_homogeneous_or_empty()) {
                    newnode->append_new_child(*it);
                    it = r->children.erase(it); /* get iterator to next child */
                } else {
                    ++it;
                }

            }
            r->append_new_child(newnode);
        }
    }
}

/* This is an implementation of algorithm 3 of [TCHP2008]_. */
void md_forest_mark_partitive_forest(md_forest &MDi, const SDData &sd,
                                     ScratchData::MDSequences &scratch) {
    size_t i = sd.first_slice_index();
    i = sd.next_slice_index(i); /* skip first slice */
    scratch.Explore.clear();
    for (; i < sd.size(); i = sd.next_slice_index(i)) {
        scratch.Marked.clear();
        scratch.Full.clear();

        for (auto it = sd.lex_label[i].begin() + sd.lex_label_offset;
                                        it != sd.lex_label[i].end() ; ++it) {
            scratch.Explore.push_back(scratch.leaves.at(*it));
        }

        while (scratch.Explore.size() > 0) {
            md_tree_node *n = scratch.Explore.front();
            scratch.Explore.pop_front();
            md_tree_node *p = n->parent;
            scratch.Full.insert(n);
            if (n->is_empty()) {
                n->label = Label::HOMOGENEOUS;
            }
            if (p) {
                scratch.Marked.insert(p);
                /* if all children of p are Full, move p to Explore */
                bool b = true;
                for (md_tree_node *c: p->children) {
                    if (scratch.Full.find(c) == scratch.Full.end()) {
                        b = false;
                        break;
                    }
                }
                if (b) {
                    scratch.Marked.erase(p);
                    scratch.Explore.push_back(p);
                }
            }
        }

        for (md_tree_node *n: scratch.Marked) {
            /* If n is SERIES or PARALLEL => gather children of n in Full below
             * the same new node A (needed only if there are >= 2 such children)
             * and the children of n not in Full below the same new node B
             * (needed only if there is >= 2 such children).
             */
            if (n->is_degenerate() && n->children.size() > 2) {
                Type t = n->type;
                md_tree_node *newnodes[2] = {
                        new md_tree_node(t, Label::HOMOGENEOUS, Flag::FLAGGED),
                        new md_tree_node(t, Label::EMPTY, Flag::UNFLAGGED)
                    };
                auto end = n->children.end();
                for (auto it = n->children.begin(); it != end; ++it){
                    bool notFull = scratch.Full.find(*it) == scratch.Full.end();
                    newnodes[notFull]->append_new_child(*it);
                }
                n->children.clear();
                for (size_t i = 0; i < 2; i++) {
                    if (newnodes[i]->children.size() == 1) {
                        n->append_new_child(newnodes[i]->children.front());
                    } else {
                        n->append_new_child(newnodes[i]);
                    }
                }
            }

            if (n->label != Label::DEAD) {
                n->label = Label::DEAD;
                /* Set flag to * for children of n that are in Full */
                for (md_tree_node *c: n->children) {
                    if (scratch.Full.find(c) != scratch.Full.end()) {
                        c->flag = Flag::FLAGGED;
                    }
                }
            }
        }
    }

    for (md_tree_node *md: MDi) {
        mark_partitive_forest_finish_inner_rec(md);
    }
}

/* This function is needed by md_forest_extract_and_sort. */
void sort_broken_nodes_recursively(md_tree_node *n,
                                   bool dead_and_broken_first) {
    if (n->is_dead_or_broken()) {
        /* If the label is not DEAD or BROKEN, no need to go deeper: they
         * will not be any DEAD or BROKEN nodes.
         */
        for (md_tree_node *c: n->children) {
            sort_broken_nodes_recursively(c, dead_and_broken_first);
        }

        if (n->is_broken()) {
            /* if dead_and_broken_first is true
             *          => put DEAD and BROKEN children at the beginning
             * if dead_and_broken_first is false
             *          => put EMPTY and HOMOGENEOUS children at the beginning
             */
            auto b = n->children.begin();
            auto end = n->children.end();
            for (auto it = n->children.begin(); it != end; ++it) {
                if (dead_and_broken_first == (*it)->is_dead_or_broken()) {
                    std::iter_swap(it, b);
                    ++b;
                }
            }
        }
    }
}

/* This function is needed by md_forest_extract_and_sort. */
void sort_dead_nodes_recursively(md_tree_node *n, bool flagged_first) {
    if (n->is_dead_or_broken()) {
        /* If the label is not DEAD or BROKEN, no need to go deeper: they
         * will not be any DEAD or BROKEN nodes.
         */
        for (md_tree_node *c: n->children) {
            sort_dead_nodes_recursively(c, flagged_first);
        }

        if (n->is_dead()) {
        /* if flagged_first is true => put flagged children at the beginning
         * if flagged_first is talse => put unflagged children at the beginning
         */
            auto b = n->children.begin();
            auto end = n->children.end();
            for (auto it = n->children.begin(); it != end; ++it) {
                if (flagged_first == ((*it)->flag == Flag::FLAGGED)) {
                    std::iter_swap(it, b);
                    ++b;
                }
            }
        }
    }
}

/* This is an implementation of algorithm 4 of [TCHP2008]_. */
void md_forest_extract_and_sort(md_forest &MDi) {
    bool is_first_slice = true;
    auto it = MDi.begin();
    while (it != MDi.end()) {
        md_tree_node *md = *it;
        /* sort children of DEAD nodes */
        sort_dead_nodes_recursively(md, is_first_slice);
        /* sort children of BROKEN nodes */
        sort_broken_nodes_recursively(md, is_first_slice);

        /* remove DEAD and BROKEN nodes */
        auto it_delete = it;
        ++it;
        while (it_delete != it) {
            md_tree_node *r = *it_delete;
            if (r->is_dead_or_broken()) {
                for (md_tree_node *c: r->children) {
                    /* propagate cc tag and slice, if set */
                    c->cc_tag = r->cc_tag != SIZE_MAX ? r->cc_tag : c->cc_tag;
                    c->slice = r->slice != SIZE_MAX ? r->slice : c->slice;
                    c->parent = nullptr;
                }
                auto insert_it = MDi.erase(it_delete);
                it_delete = r->children.begin();
                MDi.splice(insert_it, r->children);
                delete r;
            } else {
                it_delete++;
            }
        }

        is_first_slice = false;
    }
}

/*
 * This function gathers the trees of the forest in clusters. A cluster is a set
 * of trees that belongs to the same slice and (co)connected components.
 * It also computes the Left and Right of the clusters (see section 5.2 of
 * [TCHP2008]_.
 */
void md_forest_clusters_computation(const md_forest &MDi, const SDData &sd,
                                    ScratchData::Clusters &scratch) {
    scratch.clusters.clear();
    scratch.cluster_of_v.clear();
    size_t prev_cc = SIZE_MAX, prev_slice = SIZE_MAX;
    for (md_tree_node *n: MDi) {
        size_t cc = n->cc_tag;
        size_t slice = n->slice;
        int v = n->vertex; /* a vertex belonging to the module */

        if (cc == SIZE_MAX) {  /* n is alone in the cluster */
            scratch.clusters.emplace_back(1, n);  /* new cluster */
        } else {
            if (cc != prev_cc || slice != prev_slice) {
                scratch.clusters.emplace_back();  /* start new cluster */
            }
            scratch.clusters.back().push_back(n);
        }

        prev_cc = cc;
        prev_slice = slice;
        scratch.cluster_of_v[v] = scratch.clusters.size()-1;
    }

    size_t p = scratch.cluster_of_v[sd.sigma[0]];
    size_t q = scratch.clusters.size();

    /* Left and Right computation.
     *  Left(i) == i if i <= p
     *  Left(i) == Left(j) for p < i,j if clusters Ki and Kj in same slice
     *  Right(i) = p for 0 <= i <= p
     * Also compute adjacency between modules and slices (except the first one).
     */
    scratch.Left.clear();
    scratch.Right.clear();
    scratch.module_slice_adjacency.clear();
    scratch.Left.reserve(q);
    scratch.Right.reserve(q);
    for (size_t i = 0; i <= p; i++) {
        scratch.Left.push_back(i);
    }
    scratch.Right.resize(p+1, p); /* Right[i] = p for 0 <= i <= p */
    for (size_t i = p+1; i < q; i++) {
        scratch.Right.push_back(i);
    }
    for (size_t i = sd.first_slice_index(), s = 0, j = 0; i < sd.size();
                                            i = sd.next_slice_index(i), s++) {
        size_t j0 = j;
        /* Compute j, the highest index of a cluster of the current slice s */
        for (; j+1 < q && scratch.clusters[j+1].front()->slice == s; j++);

        if (s == 0) { /* nothing to do for the first slice */
            j += 1; /* skip "cluster" {x} */
        } else {
            /* for Right and module_slice_adjacency: iterate over the
             * lexicographic labels of the slice.
             */
            for (auto it = sd.lex_label[i].begin() + sd.lex_label_offset;
                                        it != sd.lex_label[i].end() ; ++it) {
                auto c = scratch.cluster_of_v.find(*it);
                if (c != scratch.cluster_of_v.end()) {
                    scratch.module_slice_adjacency.emplace(*it, s);
                    /* cluster j is adjacent to cluster containing c so Right of
                     * the cluster c is >= j
                     */
                    scratch.Right[c->second] = j;
                }
            }

            /* for Left: find the cluster of the first non adjacent module */
            size_t lp; /* lp will be the Left for all clusters of the slice */
            for (lp = 0; lp < p; lp++) {
                bool adj = true;
                for (const md_tree_node *m: scratch.clusters[lp]) {
                    if (scratch.module_slice_adjacency.find(std::make_pair(m->vertex, s)) == scratch.module_slice_adjacency.end()) {
                        adj = false;
                        break;
                    }
                }
                if (!adj) {
                    break;
                }
            }
            /* Left is lp for all clusters of the slice s */
            std::fill_n(std::back_inserter(scratch.Left), j-j0, lp);
        }
    }

}

/* This is an implementation of algorithms 5 and 6 of [TCHP2008]_. */
md_tree_node *md_forest_parse_and_assemble(md_tree_node *root, size_t p,
                                        const ScratchData::Clusters &scratch) {
    size_t q = scratch.clusters.size();
    size_t l = p;
    size_t r = p;
    while (l > 0 || r+1 < q) {
        Type t;
        size_t i;
        size_t lp, old_l = l;
        size_t rp, old_r = r;

        if (r+1 == q || (l>0 && scratch.are_clusters_non_adjacent(l-1, r+1))) {
            lp = l-1;
            rp = r;
            t = Type::SERIES;
        } else {
            lp = l;
            rp = r+1;
            t = Type::PARALLEL;
        }

        while (lp < l || r < rp) {
            if (lp < l) {
                i = l = l-1;
            } else {
                i = r = r+1;
            }
            lp = std::min(lp, scratch.Left[i]);
            rp = std::max(rp, scratch.Right[i]);
        }

        t = (r-l)-(old_r-old_l) > 1 ? Type::PRIME : t;
        md_tree_node *old_root = root;
        root = new md_tree_node(t);

        for (size_t i = l; i <= r; i++) {
            if (i == old_l) { /* add the previous root */
                root->append_new_child(old_root);
                i = old_r;
            } else {
                for (md_tree_node *m: scratch.clusters[i]) {
                    if (t != Type::PRIME && m->type == t) {
                        root->append_stolen_children_from(m);
                        delete m;
                    } else {
                        root->append_new_child(m);
                    }
                }
            }
        }
    }
    return root;
}

/* This is the main function: it corresponds to algorithms 7 of [TCHP2008]_. */
md_tree_node *corneil_habib_paul_tedder_inner_rec(const SDData &sd,
                                                  ScratchData &scratch) {
    if (sd.size() == 0) { /* empty graph */
        return nullptr;
    }

    SDData sub_sd;
    std::list<md_tree_node *> MDi;
    int x = sd.sigma[0];

    /* First create a new leaf for x */
    md_tree_node *root = scratch.new_leaf(x);

    if (sd.size() == 1) { /* graph with one vertex */
        return root;
    } else if (sd.size() == 2) { /* graph with two vertices */
        int y = sd.sigma[1];
        /* root is SERIES if there is an edge between x and y, else PARALLEL */
        Type t = sd.is_pivot_isolated() ? Type::PARALLEL : Type::SERIES;
        root = new md_tree_node(t);
        root->append_new_child(scratch.mdseq.leaves[x]);
        root->append_new_child(scratch.new_leaf(y));
        return root;
    }

    /* Now it is known that the graph has more than two vertices */

    /* Recursive calls on all the slices */
    size_t first_of_last_slice = 1; /* set to 1 to remove warning */
    for (size_t i = sd.first_slice_index(); i < sd.size();
                                            i = sd.next_slice_index(i)) {
        first_of_last_slice = i;
        sub_sd.set_to_subslice(sd, i);
        md_tree_node *md = corneil_habib_paul_tedder_inner_rec(sub_sd, scratch);
        MDi.push_back(md);
    }

    if (sd.is_pivot_isolated()) { /* x is isolated (i.e., has no neighbor) */
        md_tree_node *md = MDi.front(); /* only one slice in this case */
        if (md->type == Type::PARALLEL) {
            root = md;
        } else {
            root = new md_tree_node(Type::PARALLEL);
            root->append_new_child(md);
        }
        root->prepend_new_child(scratch.mdseq.leaves[x]);
        return root;
    }

    /* Now, it is known that x has at least one neighbor, so the first slice
     * contains the neighborhood of x.
     * If G is not connected, the last slice is all the connected components
     * that do not contains x, it should be treated separately: the tree
     * corresponding to the last slice is removed from MDi and will be added
     * back before the clusters computation.
     */
    md_tree_node *last_md = nullptr;
    if (sd.lex_label_size(first_of_last_slice) == 0) { /* not connected */
        last_md = MDi.back();
        last_md->slice = MDi.size()-1; /* remember the slice */
        MDi.pop_back();
    }

    /* Preprocessing of the sub md trees:
     *  - set the slice attribute
     *  - set the connected components tag (for clusters computation later)
     *  - set label to EMPTY and flag to UNFLAGGED on all nodes
     */
    md_forest_preprocess(MDi);

    /* Add the pivot {x} in MDi after the first slice (= the neighbors of x) */
    MDi.insert(++MDi.begin(), root);

    /* Mark partitive forest */
    md_forest_mark_partitive_forest(MDi, sd, scratch.mdseq);

    /* Extract and sort */
    md_forest_extract_and_sort(MDi);

    /* Add back the last slice if graph is not connected */
    if (last_md != nullptr) {
        MDi.push_back(last_md);
    }

    /* Postprocessing for connected (co-)components of slices: compute the
     * factoring x-m-cluster sequence
     */
    md_forest_clusters_computation(MDi, sd, scratch.clusters);

    size_t p = scratch.clusters.cluster_of_v[x];

    /* Parse and assemble */
    root = md_forest_parse_and_assemble(root, p, scratch.clusters);

    return root;
}

/*
 * It is the function exported in the pxd file. It creates the ScratchData
 * struct before calling the algorithm.
 */
md_tree_node *corneil_habib_paul_tedder_inner(const SDData &sd) {
    ScratchData tmp;
    return corneil_habib_paul_tedder_inner_rec(sd, tmp);
}
