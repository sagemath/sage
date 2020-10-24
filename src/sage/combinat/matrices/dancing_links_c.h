//*****************************************************************************
//       Copyright (C) 2008 Carlo Hamalainen <carlo.hamalainen@gmail.com>,
//
//  Distributed under the terms of the GNU General Public License (GPL)
//
//    This code is distributed in the hope that it will be useful,
//    but WITHOUT ANY WARRANTY; without even the implied warranty of
//    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//    General Public License for more details.
//
//  The full text of the GPL is available at:
//
//                  http://www.gnu.org/licenses/
//*****************************************************************************

// This file contains a C++ port of Knuth's C implementation of the Dancing
// Links algorithm. Changes have been minor, including:
//
// * the main search loop doesn't use goto statements
// * there are no hard coded limits
// * the main search algorithm can be called multiple times, and a
//   solution is saved in the public variable 'solution' of the
//   dancing_links class.

#include <sstream>
#include <set>
#include <iostream>
#include <vector>
#include <string>
#include <cassert>


using namespace std;

template <class T>
string to_string(T x)
{
    std::ostringstream stream_out;
    stream_out << x;
    return stream_out.str();
}

typedef struct node_struct {
    int tag;
    struct node_struct *left, *right;
    struct node_struct *up, *down;
    struct column_struct *col;
} node;

typedef struct column_struct {
    struct node_struct head;
    int length;
    int index;

    struct column_struct *prev, *next;
} column;

#define SEARCH_FORWARD 1
#define SEARCH_ADVANCE 2
#define SEARCH_BACKUP 3
#define SEARCH_RECOVER 4
#define SEARCH_DONE 5

class dancing_links {
    int nr_columns;
    bool search_started;

    column *root;
    column * smallest_column()
    {
        int minLength = -1;
        column *minColumn = 0;

        set<column *> seenColumns;

        for(column *p = root->next; p != root; p = p->next) {
            if (minLength == -1 || p->length < minLength) {
                minLength = p->length;
                minColumn = p;
            }
        }

        return minColumn;
    }

    vector<node*> choice;
    vector<column*> col_array;
    vector<node*> node_array;

    // For use in search() only.
    int mode;
    node *current_node;
    column *best_col;

    void add_row(vector<int> row, int row_id)
    {
        node *rowStart = 0;

        for(vector<int>::iterator i = row.begin(); i != row.end(); i++) {
            link_entry(*i + 1, row_id, &rowStart);
        }
    }

    void cover(column *c)
    {
        // Unlink c from the column list.
        column *left = c->prev;
        column *right = c->next;
        left->next = right;
        right->prev = left;

        // For each row that this column points
        // to (i.e. has 1's), cover the row too.
        for(node *row = c->head.down; row != &(c->head); row = row->down) {
            for(node *n = row->right; n != row; n = n->right) {
                node *up = n->up;
                node *down = n->down;

                up->down = down;
                down->up = up;

                n->col->length--;
            }
        }
    }


    void cover_other_columns(node *c)
    {
        for(node *p = c->right; p != c; p = p->right)
            cover(p->col);
    }


    // Links a row, we call this for each column i that contains
    // a 1. We start things by calling with *rowStart == 0
    // and this function allocates a new row, otherwise it continues the
    // row given.
    void link_entry(int i, int row_id, node **rowStart)
    {
        assert(i <= nr_columns);

        // Link this in to column i.
        node *n = (node *) malloc(sizeof(struct node_struct));
        node_array.push_back(n);

        n->tag = row_id;
        n->col = col_array[i];

        if (*rowStart == 0) *rowStart = n;

        n->down = &(col_array[i]->head);

        if (col_array[i]->length == 0) {
            n->up = &(col_array[i]->head);

            col_array[i]->head.down = n;
        } else {
            n->up = col_array[i]->head.up;
            col_array[i]->head.up->down = n;
        }

        col_array[i]->head.up = n;

        if (*rowStart != n) {
            n->left = (*rowStart)->left;
            n->right = (*rowStart);

            (*rowStart)->left->right = n;
            (*rowStart)->left = n;
        } else {
            n->left = n;
            n->right = n;
        }

        col_array[i]->length++;
    }


    // Exact reverse of cover(c).
    void uncover(column *c)
    {
        for(node *row = c->head.up; row != &(c->head); row = row->up) {
            for(node *n = row->left; n != row; n = n->left) {
                node *up = n->up;
                node *down = n->down;

                up->down = down->up = n;

                n->col->length++;
            }
        }

        column *left = c->prev;
        column *right = c->next;
        left->next = right->prev = c;
    }

    void uncover_other_columns(node *c)
    {
        for(node *p = c->left; p != c; p = p->left)
            uncover(p->col);
    }


    void save_solution()
    {
        solution.clear();

        for(vector<node*>::iterator n = choice.begin(); n != choice.end(); n++) {
            solution.push_back((*n)->tag);
        }
    }


    void setup_columns()
    {
        assert(nr_columns > 0);

        // The root column, for bookkeeping.
        root = (column *) malloc(sizeof(struct column_struct));
        root->index = 0;
        root->head.tag = 0;
        col_array.push_back(root);

        for(int i = 0; i < nr_columns; i++) {
            column *current_col = (column *) malloc(sizeof(struct column_struct));

            current_col->length = 0;
            current_col->index = i+1;

            current_col->head.tag = -1;

            // The head initially points up/down to itself.
            current_col->head.up = &(current_col->head);
            current_col->head.down = &(current_col->head);

            root->prev = current_col;
            col_array.back()->next = current_col;

            current_col->prev = col_array.back();
            current_col->next = root;

            col_array.push_back(current_col);
        }
    }
public:
    vector<int> solution;

    dancing_links()
    {
        nr_columns = -1;
        current_node = NULL;
        best_col = NULL;
        search_started = false;
    }

    ~dancing_links()
    {
        for (vector<column*>::iterator i = col_array.begin(); i != col_array.end(); i++)
            free(*i);

        for (vector<node*>::iterator i = node_array.begin(); i != node_array.end(); i++)
            free(*i);
    }

    int number_of_columns() {
        return nr_columns;
    }

    void add_rows(vector<vector<int> > rows) {
        assert(nr_columns == -1);

        // Calculate the maximum value that appears
        for(vector<vector<int> >::iterator row = rows.begin(); row != rows.end(); row++) {
            for(vector<int>::iterator r = row->begin(); r != row->end(); r++) {
                if (nr_columns < *r) nr_columns = *r;
            }
        }

        nr_columns++;

        if (nr_columns == 0) {
            return;
        }

        setup_columns();

        // Add each row
        int row_id = 0;
        for(vector<vector<int> >::iterator row = rows.begin(); row != rows.end(); row++) {
            add_row(*row, row_id);
            row_id++;
        }
    }

    bool search_is_started()
    {
        return search_started;
    }

    bool search()
    {
        search_started = true;

        if (nr_columns <= 0) {
            return false;
        }

        if (mode == SEARCH_DONE) {
            return false;
        }

        // If current_node or best_col have changed from being NULL
        // then we must have already found a solution and we are
        // entering the search() function again, looking for the
        // next solution (if any).
        if (current_node != NULL || best_col != NULL) {
            mode = SEARCH_RECOVER;
        } else {
            mode = SEARCH_FORWARD;
        }

        while(true) {
            if (mode == SEARCH_FORWARD) {
                best_col = smallest_column();
                cover(best_col);

                current_node = best_col->head.down;
                choice.push_back(current_node);

                mode = SEARCH_ADVANCE;
            }

            if (mode == SEARCH_ADVANCE) {
                if (current_node == &(best_col->head)) {
                    mode = SEARCH_BACKUP; continue;
                }

                cover_other_columns(current_node);

                if (col_array[0]->next == col_array[0]) {
                    save_solution();
                    return true;
                }
                mode = SEARCH_FORWARD; continue;
            }

            if (mode == SEARCH_BACKUP) {
                uncover(best_col);
                if ((int) choice.size() == 1) {
                    mode = SEARCH_DONE; continue;
                }

                choice.pop_back();

                current_node = choice.back();
                best_col = current_node->col;

                mode = SEARCH_RECOVER;
            }

            if (mode == SEARCH_RECOVER) {
                uncover_other_columns(current_node);

                choice.pop_back();
                current_node = current_node->down;
                choice.push_back(current_node);

                mode = SEARCH_ADVANCE; continue;
            }

            if (mode == SEARCH_DONE) {
                return false;
            }
        }
    }
};
