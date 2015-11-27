#ifndef GRID_H_
#define GRID_H_

#include "domain_define.hpp"
#include "node.hpp"
#include "cell.hpp"

#include "../algebra/space.hpp"

namespace carpio {

    template<typename COO_VALUE, typename VALUE, int DIM>
    class Grid_ {
    public:
        static const st Dim = DIM;
        static const st NumFaces = DIM + DIM;
        static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
        static const st NumNeighbors = NumFaces;

        typedef COO_VALUE cvt;
        typedef VALUE vt;
        typedef Grid_<COO_VALUE, VALUE, DIM> Self;
        typedef Grid_<COO_VALUE, VALUE, DIM> *pSelf;
        typedef Cell_<COO_VALUE, Dim> Cell;
        typedef Cell *pCell;
        typedef Data_<VALUE, Dim> Data;
        typedef Data *pData;
        typedef Node_<COO_VALUE, VALUE, DIM> Node;
        typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
        typedef typename SpaceT<pNode, Dim>::reference reference;
        typedef typename SpaceT<pNode, Dim>::const_reference const_reference;
        typedef typename SpaceT<pNode, Dim>::size_type size_type;

        typedef void (*pfunction)(pNode, utPointer);

        typedef void (*pfunction_conditional)(arrayList &, pNode, utPointer);

        typedef typename SpaceT<pNode, Dim>::iterator iterator;
        typedef typename SpaceT<pNode, Dim>::const_iterator const_iterator;

        /*
         *  data
         */
        SpaceT<pNode, Dim> nodes;

        /*
         *  constructor
         */

        Grid_(st ni, cvt ox, cvt dx, //
             st nj = 0, cvt oy = 0, cvt dy = 0, //
             st nk = 0, cvt oz = 0, cvt dz = 0) :
                nodes(ni, nj, nk) {
            st i_1d = 0;
            for (st i = 0; i < ni; i++) {
                for (st j = 0; j < ((Dim >= 2) ? nj : 1); j++) {
                    for (st k = 0; k < ((Dim == 3) ? nk : 1); k++) {
                        // new
                        nodes.at_1d(i_1d) =  //
                                new Node(
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
        ~Grid_() {
            _delete();
        }

        reference operator()(size_type i, size_type j = 0, size_type k = 0) {
            return nodes(i, j, k);
        }

        const_reference operator()(size_type i, size_type j = 0, size_type k = 0) const {
            return nodes(i, j, k);
        }

        /*
         *  size
         */
        inline st size_i() const {
            return nodes.iLen();
        }

        inline st size_j() const {
            return nodes.jLen();
        }

        inline st size_k() const {
            return (Dim < 3) ? 0 : nodes.kLen();
        }

        inline bool empty() const {
            if (nodes.size() <= 0) {
                return true;
            } else {
                return false;
            }
        }

        st get_dim() const {
            return Dim;
        }

        st size() const {
            return nodes.size();
        }

        st get_num_root() const {
            st num = 0;
            for (auto iter = nodes.begin(); iter != nodes.end(); ++iter) {
                if ((*iter) != nullptr) {
                    num++;
                }
            }
            return num;
        }

        /*
         *  iterator
         */
        iterator begin() {
            return nodes.begin();
        }

        const_iterator begin() const {
            return nodes.begin();
        }

        iterator end() {
            return nodes.end();
        }

        const_iterator end() const {
            return nodes.end();
        }

        void connect_root() {
            for (st i = 0; i < nodes.size_i(); i++) {
                for (st j = 0; j < (Dim >= 2 ? nodes.size_j() : 1); j++) {
                    for (st k = 0; k < (Dim == 3 ? nodes.size_k() : 1); k++) {
                        pNode xp = nullptr, xm = nullptr, //x
                                yp = nullptr, ym = nullptr, //y
                                zp = nullptr, zm = nullptr; //z
                        pNode cnode = nodes(i, j, k);
                        if (cnode != nullptr) {
                            // x m  and  p
                            xm = nodes.check_idx_ijk(i - 1, j, k) ? nodes(i - 1, j, k) : nullptr;
                            xp = nodes.check_idx_ijk(i + 1, j, k) ? nodes(i + 1, j, k) : nullptr;
                            ym = nodes.check_idx_ijk(i, j - 1, k) ? nodes(i, j - 1, k) : nullptr;
                            yp = nodes.check_idx_ijk(i, j + 1, k) ? nodes(i, j + 1, k) : nullptr;
                            zm = nodes.check_idx_ijk(i, j, k - 1) ? nodes(i, j, k - 1) : nullptr;
                            zp = nodes.check_idx_ijk(i, j, k + 1) ? nodes(i, j, k + 1) : nullptr;
                            cnode->set_neighbor(xm, xp, ym, yp, zm, zp);
                        }
                    }
                }
            }
        }


    };

    /*
     *  Function out of class =================================================
     */

    template<typename COO_VALUE, typename VALUE, int DIM>
    Node_<COO_VALUE, VALUE, DIM> *
    GetpNodeAt(Grid_<COO_VALUE, VALUE, DIM> *pg,
               const COO_VALUE &x,
               const COO_VALUE &y = 0,
               const COO_VALUE &z = 0) {
        typedef Node_<COO_VALUE, VALUE, DIM> *pnode;
        for (typename Grid_<COO_VALUE, VALUE, DIM>::iterator iter = pg->begin();
             iter != pg->end(); ++iter) {
            pnode pn = (*iter);
            if (pn != nullptr) {
                if (pn->cell->is_in_on(x, y, z)) {
                    return GetpNodeAt(pn, x, y, z);
                }
            }
        }
        return nullptr;

    }

    template<typename COO_VALUE, typename VALUE, int DIM>
    const Node_<COO_VALUE, VALUE, DIM> *
    GetpNodeAt(const Grid_<COO_VALUE, VALUE, DIM> *pg,
               const COO_VALUE &x,
               const COO_VALUE &y = 0,
               const COO_VALUE &z = 0) {


    }

}
#endif
