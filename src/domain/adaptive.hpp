#ifndef ADAPTIVE_HPP_
#define ADAPTIVE_HPP_


#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"


namespace carpio {
    template<typename COO_VALUE, typename VALUE, int DIM>
    void _visit_current_info(Node <COO_VALUE, VALUE, DIM> *pn, utPointer utp) {
        ArrayListV<size_t> &arr = CAST_REF(ArrayListV<size_t>*, utp);
        size_t &min_l = arr[0];
        size_t &max_l = arr[1];
        size_t &num_n = arr[2];
        size_t &num_leaf = arr[3];
        num_n++;
        if (pn->get_level() > max_l) {
            max_l = pn->get_level();
        }
        if (pn->is_leaf()) {
            min_l = (pn->get_level() < min_l) ? pn->get_level() : min_l;
            num_leaf++;
        }
    }

    template<typename COO_VALUE, typename VALUE, int DIM>
    class Adaptive {
    public:
        static const size_t Dim = DIM;
        static const size_t NumFaces = DIM + DIM;
        static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
        static const size_t NumNeighbors = NumFaces;

        typedef COO_VALUE coo_value_t;
        typedef VALUE value_t;
        typedef Grid <COO_VALUE, VALUE, DIM> grid;
        typedef Grid <COO_VALUE, VALUE, DIM> *pgrid;
        typedef Cell <COO_VALUE, Dim> cell_t;
        typedef cell_t *pcell;
        typedef Data <VALUE, Dim> data_t;
        typedef data_t *pdata;
        typedef Node <COO_VALUE, VALUE, DIM> node;
        typedef Node <COO_VALUE, VALUE, DIM> *pnode;
        typedef typename SpaceT<pnode, Dim>::reference reference;
        typedef typename SpaceT<pnode, Dim>::const_reference const_reference;
        typedef typename SpaceT<pnode, Dim>::size_type size_type;

        typedef void (*pfunction)(pnode, utPointer);

        typedef void (*pfunction_conditional)(arrayList &, pnode, utPointer);

    protected:

        size_t _c_min_l;  //current min level
        size_t _c_max_l;  //current max level
        size_t _c_num_n;  //current num of node
        size_t _c_num_leaf; //current num of leaf

        size_t _min_l;
        size_t _max_l;
        //size_t _min_num_n;
        //size_t _max_num_n;

        pgrid _grid;

        /**
         *  protected function
         */


        void _get_current_info() {
            ArrayListV<size_t> arr_info(4);
            for (int i = 0; i < _grid->size(); i++) {
                pnode pn = _grid->nodes.at_1d(i);
                if (pn != nullptr) {
                    //
                    pn->traversal(_visit_current_info, &arr_info);
                }
            }
            _c_min_l = arr_info[0];
            _c_max_l = arr_info[1];
            _c_num_n = arr_info[2];
            _c_num_leaf = arr_info[3];
        }

    public:
        Adaptive(Grid <COO_VALUE, VALUE, DIM> *pg, size_t minl = 1, size_t maxl = 5) {
            size_t n_root = pg->get_num_root();
            _min_l = minl;
            _max_l = maxl;
            _grid = pg;
            ASSERT(minl <= maxl);
            _get_current_info();
        }

        void adapt() {
            std::function<void(pnode,size_t)> fun = [](pnode pn, size_t min_l){
                if(pn->is_leaf()&& pn->get_level()<min_l){
                    pn->new_full_child();
                }
            };
            for (int i = 0; i < _grid->size(); i++) {
                pnode pn = _grid->nodes.at_1d(i);
                if (pn != nullptr) {
                    pn->traversal(fun,_min_l);
                }
            }
        }


    };


}

#endif