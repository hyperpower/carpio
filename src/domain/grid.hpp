#ifndef GRID_H_
#define GRID_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

#include "../algebra/space.hpp"

namespace carpio {

    template<typename COO_VALUE, typename VALUE, int DIM>
    class Grid {
    public:
        static const size_t Dim = DIM;
        static const size_t NumFaces = DIM + DIM;
        static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
        static const size_t NumNeighbors = NumFaces;

        typedef COO_VALUE coo_value_t;
        typedef VALUE value_t;
        typedef Grid<COO_VALUE, VALUE, DIM> self;
        typedef Grid<COO_VALUE, VALUE, DIM> *pself;
        typedef Cell<COO_VALUE, Dim> cell_t;
        typedef cell_t *pcell;
        typedef Data<VALUE, Dim> data_t;
        typedef data_t *pdata;
        typedef Node<COO_VALUE, VALUE, DIM> node;
        typedef Node<COO_VALUE, VALUE, DIM> *pnode;
        typedef typename SpaceT<pnode, Dim>::reference reference;
        typedef typename SpaceT<pnode, Dim>::const_reference const_reference;
        typedef typename SpaceT<pnode, Dim>::size_type size_type;
        typedef void (*pfunction)(pnode, utPointer);

        typedef void (*pfunction_conditional)(arrayList &, pnode, utPointer);

        /*
         *  data
         */
        SpaceT<pnode, Dim> nodes;

        /*
         *  constructor
         */

        Grid(size_t ni, coo_value_t ox, coo_value_t dx, //
             size_t nj = 0, coo_value_t oy = 0, coo_value_t dy = 0, //
             size_t nk = 0, coo_value_t oz = 0, coo_value_t dz = 0) :
                nodes(ni, nj, nk) {
            size_t i_1d = 0;
            for (size_t i = 0; i < ni; i++) {
                for (size_t j = 0; j < ((Dim >= 2) ? nj : 1); j++) {
                    for (size_t k = 0; k < ((Dim == 3) ? nk : 1); k++) {
                        // new
                        nodes.at_1d(i_1d) =  //
                                new node(
                                        nullptr, // father
                                        0, // type
                                        0, //level
                                        i_1d, //root idx
                                        0, //path
                                        ox + (i + 0.5) * dx, 0.5 * dx, //x
                                        oy + (j + 0.5) * dy, 0.5 * dy, //y
                                        oz + (k + 0.5) * dz, 0.5 * dz); //z
                        i_1d++;
                    }
                }
            }
        }

    protected:
        void _delete() {
            for (int i = 0; i < nodes.size(); i++) {
                if (nodes.at_1d(i) != nullptr) {
                    delete nodes.at_1d(i);
                }
            }
        }

    public:
        ~Grid() {
            _delete();
        }

        reference operator()(size_type i, size_type j= 0, size_type k= 0){
            return nodes(i,j,k);
        }
        const_reference operator()(size_type i, size_type j= 0, size_type k= 0) const{
            return nodes(i,j,k);
        }

        /*
         *  size
         */
        inline size_t size_i() const {
            return nodes.iLen();
        }

        inline size_t size_j() const {
            return nodes.jLen();
        }

        inline size_t size_k() const {
            return (Dim < 3) ? 0 : nodes.kLen();
        }

        inline bool empty() const {
            if (nodes.size() <= 0) {
                return true;
            } else {
                return false;
            }
        }

        size_t get_dim() const {
            return Dim;
        }

        size_t size() const {
            return nodes.size();
        }

        void connect_root() {
            for (size_t i = 0; i < nodes.size_i(); i++) {
                for (size_t j = 0; j < (Dim >= 2 ? nodes.size_j() : 1); j++) {
                    for (size_t k = 0; k < (Dim == 3 ? nodes.size_k() : 1); k++) {
                        pnode xp = nullptr, xm = nullptr, //x
                                yp = nullptr, ym = nullptr, //y
                                zp = nullptr, zm = nullptr; //z
                        pnode cnode = nodes(i, j, k);
                        if (cnode != nullptr) {
                            // x m  and  p
                            xm = nodes.check_idx_ijk(i - 1, j, k) ? nodes(i - 1, j, k) : nullptr;
                            xp = nodes.check_idx_ijk(i + 1, j, k) ? nodes(i + 1, j, k) : nullptr;
                            ym = nodes.check_idx_ijk(i, j - 1, k) ? nodes(i, j - 1, k) : nullptr;
                            yp = nodes.check_idx_ijk(i, j + 1, k) ? nodes(i, j + 1, k) : nullptr;
                            zm = nodes.check_idx_ijk(i, j, k - 1) ? nodes(i, j, k - 1) : nullptr;
                            zp = nodes.check_idx_ijk(i, j, k + 1) ? nodes(i, j, k + 1) : nullptr;
                        }
                        cnode->set_neighbor(xm, xp, ym, yp, zm, zp);
                    }
                }
            }
        }

    };

}
#endif
