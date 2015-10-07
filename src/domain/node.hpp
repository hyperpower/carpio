#ifndef NODE_H_
#define NODE_H_

#include "../carpio_define.hpp"
#include "domain_define.hpp"
#include "cell.hpp"
#include "data.hpp"

#include <math.h>

namespace carpio {

    enum NodeIdx {
        //=========================
        //   y
        //   |
        //   ---------------
        //   |  PM  |  PP  |
        //   |  2   |  3   |
        //   |  WN  |  NE  |
        //   ---------------
        //   |  SW  |  SE  |
        //   |  0   |  1   |
        //   |  MM  |  MP  |
        //   ------------------->x
        //
        //   ---------------  ---------------
        //   |  MPM |  MPP |  |  PPM |  PPP |
        //   |  2   |  3   |  |  6   |  7   |
        //   |  WNB |  NEB |  |  WNF |  NEF |
        //   ---------------  ---------------
        //   |  SWB |  SEB |  |  SWP |  SEP |
        //   |  0   |  1   |  |  4   |  5   |
        //   |  MMM |  MMP |  |  PMM |  PMP |
        //   ---------------  ---------------
        //=========================

        //2D
                _MM_ = 0,
        _MP_ = 1,
        _PM_ = 2,
        _PP_ = 3,
        //3D
                _MMM_ = 0,
        _MMP_ = 1,
        _MPM_ = 2,
        _MPP_ = 3,
        _PMM_ = 4,
        _PMP_ = 5,
        _PPM_ = 6,
        _PPP_ = 7,
    };

    inline bool is_x_p(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 6) == 7;
    }

    inline bool is_x_m(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 6) == 6;
    }

    inline bool is_y_p(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 5) == 7;
    }

    inline bool is_y_m(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 5) == 5;
    }

    inline bool is_z_p(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 3) == 7;
    }

    inline bool is_z_m(size_t i) {
        ASSERT(i >= 0 && i < 8);
        return (i | 3) == 3;
    }
#define _TEMPLATE_COOV_V_DIM_ template<typename COO_VALUE, typename VALUE, int DIM>
#define _COOV_V_DIM_ COO_VALUE, VALUE, DIM

    _TEMPLATE_COOV_V_DIM_
    class Node {
    public:
        static const size_t Dim = DIM;
        static const size_t NumFaces = DIM + DIM;
        static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
        static const size_t NumNeighbors = NumFaces;
        static const size_t NumChildren = NumVertexes;

        typedef COO_VALUE coo_value_t;
        typedef VALUE value_t;
        typedef Node<_COOV_V_DIM_> Self;
        typedef Node<_COOV_V_DIM_> *pSelf;
        typedef Cell<COO_VALUE, Dim> Cell_;
        typedef Cell_ *pCell;
        typedef Data<VALUE, Dim> Data_;
        typedef Data_ *pData;
        typedef Self Node_;
        typedef Self *pNode;

        typedef void (*pFun)(pNode, utPointer);

        typedef void (*pFun_Conditional)(arrayList &, pNode, utPointer);

    protected:
        typedef size_t st;
        typedef COO_VALUE cvt;
        typedef VALUE vt;
        //
        int _node_type;
        st _level;
        st _root_idx;
        st _path;
    public:
        pNode father;
        pNode child[NumChildren];
        pNode neighbor[NumNeighbors];
        pCell cell;
        pData data;

    protected:
        int _height(const pNode Current) const {
            if (Current == nullptr) {
                return 0;
            }
            if (!Current->HasChild()) {
                return 0;
            } else {
                ArrayListV<st> arrh(NumChildren);
                for (st i = 0; i < this->NumChildren; ++i) {
                    arrh = _height(Current->child[i]);
                }
                return 1 + arrh.max();
            }
        }

        void _traversal_conditional(pNode pn, pFun_Conditional pfun_con,
                                    pFun pfun, utPointer utp) {
            if (pn == nullptr) {
                return;
            } else {
                (*pfun)(pn, utp);
                if (pn->HasChild()) {
                    arrayList avt(NumChildren);
                    pfun_con(avt, pn, utp);
                    for (int i = 0; i < NumChildren; i++) {
                        pNode c = pn->child[i];
                        if (c != nullptr && avt[i] == 1) {
                            _traversal_conditional(c, pfun_con, pfun, utp);
                        }
                    }
                }
            }
        }

        void _traversal(pNode pn, pFun pfun, utPointer utp) {
            if (pn == nullptr) {
                return;
            } else {
                (*pfun)(pn, utp);
                if (pn->has_child()) {
                    for (int i = 0; i < NumChildren; i++) {
                        pNode c = pn->child[i];
                        if (c != nullptr) {
                            _traversal(c, pfun, utp);
                        }
                    }
                }
            }
        }

        template <class Ret, class Args>
        void _traversal(pNode pn, std::function<Ret(pNode, Args)> fun, Args& args){
            if (pn == nullptr) {
                return;
            } else {
                fun(pn,args);
                if (pn->has_child()) {
                    for (int i = 0; i < NumChildren; i++) {
                        pNode c = pn->child[i];
                        if (c != nullptr) {
                            _traversal(c, fun, args);
                        }
                    }
                }
            }
        }

    public:
        /*
         *  constructor
         */
        Node(pNode f, int nt, st level, st root_idx, st path,    //
             const vt &x, const vt &dhx, //
             const vt &y = 0.0, const vt &dhy = 0.0, //
             const vt &z = 0.0, const vt &dhz = 0.0) {
            _node_type = nt;
            _level = level;
            cell = new Cell_(x, dhx, y, dhy, z, dhz);
            father = f;
            _root_idx = root_idx;
            _path = path;

            data = nullptr;
            for (int i = 0; i < this->NumChildren; i++) {
                child[i] = nullptr;
            }
            for (int i = 0; i < this->NumNeighbors; i++) {
                neighbor[i] = nullptr;
            }
        }
        /*
         *  delete
         */
    protected:
        /*
         *  before using this function, making sure that this node is a leaf
         */
        void _DeleteLeaf() {
            pNode f = this->father;
            if (f != nullptr) {
                f->child[get_idx()] = nullptr;
            }
            delete cell;
            if (data != nullptr) {
                delete data;
            }
        }

        void _Delete(pNode pn) {
            if (pn == nullptr) {
                return;
            }
            if (pn->has_child()) {
                for (int i = 0; i < NumChildren; i++) {
                    pNode ch = pn->child[i];
                    if (ch != nullptr) {
                        _Delete(ch);
                    }
                }
            } else { // is leaf
                pn->_DeleteLeaf();
            }
        }

    public:
        ~Node() {
            _Delete(this);
        }

        /*
         * type
         */
        inline int get_type() const {
            return _node_type;
        }

        inline void set_type(int type) {
            _node_type = type;
        }

        inline st get_level() const {
            return _level;
        }

        inline st get_idx() const {
            return (_path >> int(pow(Dim, _level))) & (NumVertexes - 1);
        }

        inline st get_path() const {
            return _path;
        }

        inline st get_root_idx() const {
            return _root_idx;
        }

        inline st height() const {
            return this->_height(this);
        }

        /*
         *  child
         */
        inline bool has_child() const {
            for (st i = 0; i < this->NumChildren; ++i) {
                if (this->child[i] != nullptr) {
                    return true;
                }
            }
            return false;
        }

        inline bool is_leaf() const {
            return !has_child();
        }

        inline bool has_child(st idx) const {
            return this->child[idx] != nullptr;
        }

        inline bool is_root() const {
            if (this->father == nullptr) {
                return true;
            } else {
                return false;
            }
        }

        inline bool is_full_child() const {
            bool res = this->child[0] != nullptr;
            for (st i = 1; i < this->NumChildren; ++i) {
                res = res && (this->child[i] != nullptr);
            }
            return res;
        }

        inline st count_children() const {
            st res = 0;
            for (st i = 0; i < this->NumChildren; ++i) {
                res += (this->child[i] != nullptr) ? 1 : 0;
            }
            return res;
        }

        /*
         *  new
         */
        void new_full_child() {
            if (!has_child()) {
                st ltmp = _level + 1;
                vt nhdx = this->cell->get_hd(_X_) * 0.5;
                vt nhdy = this->cell->get_hd(_Y_) * 0.5;
                vt nhdz = this->cell->get_hd(_Z_) * 0.5;
                vt cx = this->cell->get(_C_, _X_);
                vt cy = this->cell->get(_C_, _Y_);
                vt cz = this->cell->get(_C_, _Z_);
                for (st i = 0; i < this->NumChildren; ++i) {
                    pNode f = this;
                    int nt = 1;
                    st l = ltmp;
                    st ridx = _root_idx;
                    st npath = (i << int((pow(Dim, l)))) + _path;
                    this->child[i] = new Node_( //
                            f, nt, l, ridx, npath, //
                            cx + (is_x_p(i) ? nhdx : -nhdx), nhdx, //
                            cy + (is_y_p(i) ? nhdx : -nhdx), nhdy, //
                            cz + (is_z_p(i) ? nhdx : -nhdx), nhdz);
                }
            }
        }

        /*
         *  neighbor find
         */
        void set_neighbor(
                pNode xm, pNode xp, //x
                pNode ym, pNode yp, //y
                pNode zm, pNode zp) { //z
            //       yp 3
            //      ______
            //     |      |
            //xm 0 |      | xp 1
            //     |______|
            //       ym 2
            neighbor[0] = xm;
            neighbor[1] = xp;
            if (Dim >= 2) {
                neighbor[2] = ym;
                neighbor[3] = yp;
            }
            if (Dim == 3) {
                neighbor[4] = zm;
                neighbor[5] = zp;
            }
        }

        inline bool is_adjacent(const Direction &d) const {
            // Direction on x y or z
            unsigned short hi = d >> 3;
            return ((hi & get_idx()) ^ (hi & d)) == 0;
        }

        inline st reflect(const Direction &d) const {
            // Direction on x y or z
            return get_idx() ^ (d >> 3);
        }

        inline bool has_diagonal_sibling(const Direction &d) const {
            return (get_idx() ^ (d >> 3)) == (d & 7);
        }

        inline bool is_out_corner(const Direction &d) const {
            return get_idx() == (d & 7);
        }

        inline st out_common_direction(const Direction &d) const {

            // return direction on x y or z
            st hi = d >> 3;
            st low = d & 7;
            return (((low ^ get_idx()) ^ hi) << 3) + low;
        }

    protected:
        pNode _get_root_neighbor_step(const Direction &d1, const Direction &d2) const {
            pNode pstep1 = this->get_root_neighbor(d1);
            if (pstep1 == nullptr) {
                return nullptr;
            } else {
                return pstep1->get_root_neighbor(d2);
            }
        }

        pNode _get_root_neighbor_step(const Direction &d1,
                                      const Direction &d2,
                                      const Direction &d3) const {
            pNode pstep = this->get_root_neighbor(d1);
            if (pstep == nullptr) {
                return nullptr;
            } else {
                pstep = pstep->get_root_neighbor(d2);
                if (pstep == nullptr) {
                    return nullptr;
                } else {
                    return pstep->get_root_neighbor(d3);
                }
            }
        }

        pNode _get_root_neighbor_path(const Direction &d1,
                                      const Direction &d2) const {
            pNode pt = this->_get_root_neighbor_step(d1, d2);
            if (pt != nullptr) {
                return pt;
            } else {
                return this->_get_root_neighbor_step(d2, d1);
            }
        }

        pNode _get_root_neighbor_path(const Direction &d1,
                                      const Direction &d2,
                                      const Direction &d3) const {
            pNode pt = this->_get_root_neighbor_step(d1, d2, d3);
            if (pt != nullptr) {
                return pt;
            }
            pt = this->_get_root_neighbor_step(d1, d3, d2);
            if (pt != nullptr) {
                return pt;
            }
            pt = _get_root_neighbor_step(d2, d1, d3);
            if (pt != nullptr) {
                return pt;
            }
            pt = _get_root_neighbor_step(d2, d3, d1);
            if (pt != nullptr) {
                return pt;
            }
            pt = _get_root_neighbor_step(d3, d1, d2);
            if (pt != nullptr) {
                return pt;
            }
            pt = _get_root_neighbor_step(d3, d2, d1);
            if (pt != nullptr) {
                return pt;
            }
            return nullptr;
        }

    public:
        pNode get_root_neighbor(const Direction &d) const {
            unsigned short hi = HI(d);
            ASSERT(hi != 0);
            switch (HI(d)) {
                case 1:   //x
                    return this->neighbor[hi & LO(d)];
                case 2:   //y
                    return this->neighbor[((hi & LO(d)) >> 1) + hi];
                case 3:   //xy
                    return _get_root_neighbor_path(8 + LO(d), 16 + LO(d));
                case 4:   //z
                    return this->neighbor[((hi & LO(d)) >> 2) + hi];
                case 5:
                    return _get_root_neighbor_path(32 + LO(d), 8 + LO(d));
                case 6:
                    return _get_root_neighbor_path(16 + LO(d), 32 + LO(d));
                case 7:
                    return _get_root_neighbor_path(8 + LO(d), 16 + LO(d), 32 + LO(d));
                default:
                    return nullptr;
            }
        }
        /*
         *  Traverse
         */
        void traversal(pFun pfun, utPointer utp) {
            this->_traversal(this, pfun, utp);
        }
        template <class Args>
        void traversal(std::function<void(pNode, Args)> fun, Args& args){
            this->_traversal(this, fun, args);
        }

    };

    /*
     *  functions out of class
     */
    _TEMPLATE_COOV_V_DIM_
    int GetDataIdx(const Node<COO_VALUE, VALUE, DIM> *pn) {
        ASSERT(pn != nullptr);
        return pn->data->get_idx();
    }


}//

#endif /* NODE_H_ */
