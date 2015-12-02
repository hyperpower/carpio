#ifndef ADAPTIVE_HPP_
#define ADAPTIVE_HPP_


#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"
#include "grid.hpp"


namespace carpio {
    template<typename COO_VALUE, typename VALUE, int DIM>
    void _visit_current_info(Node_ <COO_VALUE, VALUE, DIM> *pn, utPointer utp) {
        ArrayListV<st> &arr = CAST_REF(ArrayListV<st>*, utp);
        st &min_l = arr[0];
        st &max_l = arr[1];
        st &num_n = arr[2];
        st &num_leaf = arr[3];
        num_n++;
        if (pn->get_level() > max_l) {
            max_l = pn->get_level();
        }
        if (pn->is_leaf()) {
            min_l = (pn->get_level() < min_l) ? pn->get_level() : min_l;
            num_leaf++;
        }
    }

    template<typename COO_VALUE, typename VALUE, st DIM>
    class Adaptive {
    public:
        static const st Dim = DIM;
        static const st NumFaces = DIM + DIM;
        static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
        static const st NumNeighbors = NumFaces;

        typedef COO_VALUE coo_value_t;
        typedef VALUE value_t;
        typedef Grid_ <COO_VALUE, VALUE, DIM> grid;
        typedef Grid_ <COO_VALUE, VALUE, DIM> *pgrid;
        typedef Cell_ <COO_VALUE, Dim> cell_t;
        typedef cell_t *pcell;
        typedef Data_ <VALUE, Dim> data_t;
        typedef data_t *pdata;
        typedef Node_ <COO_VALUE, VALUE, DIM> Node;
        typedef Node_ <COO_VALUE, VALUE, DIM> *pNode;
        typedef typename SpaceT<pNode, Dim>::reference reference;
        typedef typename SpaceT<pNode, Dim>::const_reference const_reference;
        typedef typename SpaceT<pNode, Dim>::size_type size_type;

        typedef void (*pfunction)(pNode, utPointer);

        typedef void (*pfunction_conditional)(arrayList &, pNode, utPointer);

    protected:

        st _c_min_l;  //current min level
        st _c_max_l;  //current max level
        st _c_num_n;  //current num of node
        st _c_num_leaf; //current num of leaf

        st _min_l;
        st _max_l;
        //st _min_num_n;
        //st _max_num_n;

        pgrid _grid;

        /**
         *  protected function
         */


        void _update_current_info() {
            ArrayListV<st> arr_info(4);
            for (int i = 0; i < _grid->size(); i++) {
                pNode pn = _grid->nodes.at_1d(i);
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
        Adaptive(Grid_ <COO_VALUE, VALUE, DIM> *pg,  //the pointer of grid
        		st minl = 1,  //min level
        		st maxl = 5   //max level
        		) {
            st n_root = pg->get_num_root();
            _min_l = minl;
            _max_l = maxl;
            _grid = pg;
            ASSERT(minl <= maxl);
            _update_current_info();
        }

        void adapt() {
            std::function<void(pNode,st)> fun = [](pNode pn, st min_l){
                if(pn->is_leaf()&& pn->get_level()<min_l){
                    pn->new_full_child();
                }
            };
            for (int i = 0; i < _grid->size(); i++) {
                pNode pn = _grid->nodes.at_1d(i);
                if (pn != nullptr) {
                    pn->traversal(fun,_min_l);
                }
            }
        }
    };


}

#endif
