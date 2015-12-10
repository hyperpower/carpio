#ifndef NODE_H_
#define NODE_H_

#include "../carpio_define.hpp"
#include "domain_define.hpp"
#include "cell.hpp"
#include "data.hpp"
#include <functional>

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

inline bool is_x_p(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 7;
}

inline bool is_x_m(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 6) == 6;
}

inline bool is_y_p(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 7;
}

inline bool is_y_m(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 5) == 5;
}

inline bool is_z_p(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 7;
}

inline bool is_z_m(st i) {
	ASSERT(i >= 0 && i < 8);
	return (i | 3) == 3;
}

static const int _E_ = 1;
static const int _F_C_ = 2;
static const int _C_F_ = 3;

#define _TEMPLATE_COOV_V_DIM_ template<typename COO_VALUE, typename VALUE, int DIM>
#define _COOV_V_DIM_ COO_VALUE, VALUE, DIM

_TEMPLATE_COOV_V_DIM_ class Node_ {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;
	static const st NumChildren = NumVertexes;

	typedef COO_VALUE cvt;
	typedef VALUE vt;
	typedef VALUE& ref_vt;
	typedef const VALUE& const_ref_vt;
	typedef Node_<_COOV_V_DIM_> Self;
	typedef Node_<_COOV_V_DIM_> *pSelf;
	typedef Cell_<COO_VALUE, Dim> Cell;
	typedef Cell_<COO_VALUE, Dim> *pCell;
	typedef Data_<VALUE, Dim> Data;
	typedef Data *pData;
	typedef Self Node;
	typedef Self *pNode;
	typedef const Self const_Node;
	typedef const Self* const_pNode;

	typedef void (*pFun)(pNode, utPointer);

	typedef void (*pFun_Conditional)(arrayList &, pNode, utPointer);

protected:
	//
		int _node_type;
		st _level;
		st _root_idx;
		st _idx;
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
			if (!Current->has_child()) {
				return 0;
			} else {
				ArrayListV<st> arrh(NumChildren);
				for (st i = 0; i < this->NumChildren; ++i) {
					arrh = _height(Current->child[i]);
				}
				return 1 + arrh.max();
			}
		}

		template<class Ret1, class Ret2, class Args1, class Args2>
		void _traversal_conditional(pNode pn,
				std::function<Ret1(bool[], pNode, Args1)> fun_con,
				Args1 argsc,
				std::function<Ret2(pNode, Args2)> fun,
				Args2 args) {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
				if (pn->has_child()) {
					bool ist[NumChildren];
					for (int i = 0; i < NumChildren; i++) {
						ist[i] = false;
					}
					fun_con(ist, pn, argsc);
					for (int i = 0; i < NumChildren; i++) {
						pNode c = pn->child[i];
						if (c != nullptr && ist[i]) {
							_traversal_conditional(c, fun_con, argsc, fun, args);
						}
					}
				}
			}
		}

		void _traversal_conditional(pNode pn, pFun_Conditional pfun_con,
				pFun pfun, utPointer utp) {
			if (pn == nullptr) {
				return;
			} else {
				(*pfun)(pn, utp);
				if (pn->has_child()) {
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

		template<class Ret, class Args>
		void _traversal(pNode pn, std::function<Ret(pNode, Args)> fun, Args &args) {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
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

		template<class Ret, class Args>
		void _traversal(const_pNode pn, std::function<Ret(const_pNode, Args)> fun, Args &args) const {
			if (pn == nullptr) {
				return;
			} else {
				fun(pn, args);
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
		Node_(pNode f, int nt, st level, st root_idx, st path,    //
				const vt &x, const vt &dhx,//
				const vt &y = 0.0, const vt &dhy = 0.0,//
				const vt &z = 0.0, const vt &dhz = 0.0) {
			_node_type = nt;
			_level = level;
			cell = new Cell(x, dhx, y, dhy, z, dhz);
			father = f;
			_root_idx = root_idx;
			_idx = path;

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
		~Node_() {
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
			//return (_path >> int(pow(Dim, _level))) & (NumVertexes - 1);
			return _idx;
		}

		inline st get_path() const {
			return _idx;
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

		inline bool is_alone() const {
			return this->is_leaf()&&this->is_root();
		}

		/*
		 *  count
		 */
		inline st count_children() const {
			st res = 0;
			for (st i = 0; i < this->NumChildren; ++i) {
				res += (this->child[i] != nullptr) ? 1 : 0;
			}
			return res;
		}

		inline st count_all() const {
			st res = 0;
			std::function<void(const_pNode, st)> fun = [&res](const_pNode pn, st dummy) {
				res++;
			};
			this->traversal(fun, res);
			return res;
		}

		inline st count_leaf() const {
			st res = 0;
			std::function<void(const_pNode, st)> fun = [&res](const_pNode pn, st dummy) {
				if(pn->is_leaf()) {
					res++;
				}
			};
			this->traversal(fun, res);
			return res;
		}

		inline st count_level(st le) const {
			st res = 0;
			std::function<void(const_pNode, st)> fun = [&res](const_pNode pn, st le) {
				if(pn->get_level()==le) {
					res++;
				}
			};
			this->traversal(fun, le);
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
					//st npath = (i << int((pow(Dim, l)))) + _path;
					st npath = i;
					this->child[i] = new Node_(//
							f, nt, l, ridx, npath,//
							cx + (is_x_p(i) ? nhdx : -nhdx), nhdx,//
							cy + (is_y_p(i) ? nhdx : -nhdx), nhdy,//
							cz + (is_z_p(i) ? nhdx : -nhdx), nhdz);
				}
			}
		}

		/*
		 *  make sure point is in this node
		 */
		st which_child(const COO_VALUE &x,
				const COO_VALUE &y = 0,
				const COO_VALUE &z = 0) {
			st idx = 0;
			if (Dim == 3) {
				idx = IsInRange(this->cell->get(_M_, _Z_), z, this->cell->get(_C_, _Z_), _cc_) ? 0 : 1;
			}
			idx = idx << 1;
			if (Dim >= 2) {
				idx = idx + (IsInRange(this->cell->get(_M_, _Y_), y, this->cell->get(_C_, _Y_), _cc_) ? 0 : 1);
			}
			idx = idx << 1;
			idx = idx + (IsInRange(this->cell->get(_M_, _X_), x, this->cell->get(_C_, _X_), _cc_) ? 0 : 1);
			return idx;
		}

		/*
		 *  neighbor find
		 */
		void set_neighbor(
				pNode xm, pNode xp, //x
				pNode ym, pNode yp,//y
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
			unsigned short hi = d >> 3;
			return ((get_idx() ^ hi) & hi) == (LO(d) & hi);
		}

		inline bool is_out_corner(const Direction &d) const {
			return (get_idx() & HI(d)) == (HI(d) & LO(d));
		}

		inline Direction out_common_direction(const Direction &d) const {
			// return direction on x y or z
			static const unsigned short COMMON_AXES[4][4] = { {0, 2, 1, 0},
				{	2, 0, 0, 1},
				{	1, 0, 0, 2},
				{	0, 1, 2, 0}};
			std::function<unsigned short(unsigned short)> to2bit = [](unsigned short num) {
				unsigned short z = (num & 4) >> 2;
				return ((num & 1) << 1) + z;
			};
			std::function<unsigned short(unsigned short)> to3bit = [](unsigned short num) {
				unsigned short z = (num & 1) << 2;
				return z + ((num & 2) >> 1);
			};
			unsigned short hi = d >> 3;
			unsigned short lo = LO(d) & hi;
			unsigned short id = get_idx() & hi;
			unsigned short com;
			switch (hi) {
				case 3: //xy
				com = COMMON_AXES[id][lo];
				break;
				case 6://yz
				com = COMMON_AXES[id >> 1][lo >> 1];
				com = com << 1;
				break;
				case 5: { //zx
					com = COMMON_AXES[to2bit(id)][to2bit(lo)];
					com = to3bit(com);
					break;
				}
				default:
				ASSERT(false);
			}
			return (com << 3) + (lo & com);
		}

		inline st diagonal_idx(Direction d) const {
			return get_idx() ^ HI(d);
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
				case 2://y
				return this->neighbor[((hi & LO(d)) >> 1) + hi];
				case 3://xy
				return _get_root_neighbor_path(8 + LO(d), 16 + LO(d));
				case 4://z
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

		pNode get_adj_neighbor(const pNode Current, Direction d) {
			// face direction
			std::cout << "--------" << "\n";
			std::cout << "d n     " << d << "\n";
			std::cout << "idx adj " << Current->get_idx() << "\n";
			pNode ca = nullptr;//common ancestor
			if (Current->father != nullptr
					&& Current->is_adjacent(d)) {
				ca = get_adj_neighbor(Current->father, d);
			} else {
				ca = Current->father;
			}
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->reflect(d)];
			} else if (ca == nullptr) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_cor_neighbor(const pNode Current, Direction d) {
			std::cout << "--------" << "\n";
			std::cout << "d n     " << d << "\n";
			std::cout << "idx cor " << Current->get_idx() << "\n";
			pNode ca = nullptr;    //common ancestor
			int flag = 0;
			if (Current->father != nullptr &&
					!Current->has_diagonal_sibling(d)) {
				//Find a common ancestor
				if (Current->is_out_corner(d)) {
					std::cout << "cor " << Current->get_idx() << "\n";
					ca = get_cor_neighbor(Current->father, d);
				} else {
					std::cout << "adj " << Current->get_idx() << "\n";
					ca = get_adj_neighbor(Current->father,
							Current->out_common_direction(d));
				}
			} else {
				std::cout << "dia " << Current->get_idx() << "\n";
				flag = 1;
				ca = Current->father;
			}
			//Follow opposite path to locate the neighbor
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->diagonal_idx(d)];
			} else if (ca == nullptr && flag == 1) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_neighbor_adj_cor(const pNode Current, //
				Direction d) {
			if (IsFaceDirection(d)) {
				return get_adj_neighbor(Current, d);
			}
			if (IsPlaneDirection(d)) {
				return get_cor_neighbor(Current, d);
			}
			return nullptr;
		}

		pNode get_cor_neighbor_xyz(const pNode Current, Direction d) {
			//std::cout << "xyz " << Current->get_idx() << "\n";
			pNode ca = nullptr;//common ancestor
			if (Current->father != nullptr &&
					Current->is_out_corner(d)) {
				ca = get_cor_neighbor_xyz(Current->father, d);
			} else {
				Direction idx = Current->get_idx();
				if (idx==((~LO(d))&7)) {
					ca = Current->father;
				} else {
					Direction nd = (((~(idx ^ (d & 7))) & 7) << 3) + idx;
					ca = get_neighbor_adj_cor(Current->father, nd);
				}
			}
			pNode pt = nullptr;
			if (ca != nullptr && ca->has_child()) {
				pt = ca->child[Current->diagonal_idx(d)];
			} else if (ca == nullptr) {
				pt = Current->get_root_neighbor(d);
			} else {
				pt = ca;
			}
			return pt;
		}

		pNode get_neighbor(const pNode Current, Direction d) {
			ASSERT(d > 7);
			Direction nd = d & 63;
			if ( IsXYZDirection(nd)) {
				return get_cor_neighbor_xyz(Current, nd);
			} else {
				return get_neighbor_adj_cor(Current, nd);
			}
		}

		pNode get_neighbor(Direction d) {
			return get_neighbor(this, d);
		}

		/*
		 *  Traverse
		 */
		void traversal(pFun pfun, utPointer utp) {
			this->_traversal(this, pfun, utp);
		}

		template<class Args>
		void traversal(std::function<void(pNode, Args)> fun, Args &args) {
			this->_traversal(this, fun, args);
		}

		template<class Args>
		void traversal(std::function<void(const_pNode, Args)> fun, Args &args) const {
			this->_traversal(this, fun, args);
		}

		template<class Ret1, class Ret2, class Args1, class Args2>
		void traversal_conditional(std::function<Ret1(bool[], pNode, Args1)> fun_con,
				Args1 argsc,
				std::function<Ret2(pNode, Args2)> fun,
				Args2 args) {
			this->_traversal_conditional(this, fun_con, argsc, fun, args);
		}

		/*
		 *  overload the function of cell
		 */
		bool is_in_on(const vt x,
				const vt y = 0,
				const vt z = 0) const {
			bool res = this->cell->is_in_on(x,y,z);
			if(this->is_leaf()||res == false) {
				return res;
			}
			std::function<void(const_pNode, int)> fun =
			[&res,&x,&y,&z](const_pNode pn, int dummy) {
				if (pn->is_leaf() && res==false) {
					res = pn->cell->is_in_on(x,y,z);
				}
			};
			int dummy = 1;
			this->_traversal(this, fun, dummy);
			return res;
		}
		cvt cp(Axes axes) const {  //center point
			return this->cell->get(_C_, axes);
		}
		cvt d(Axes axes) const {  //center point
			return this->cell->get_d(axes);
		}
		cvt p(Orientation ori, Axes axes) const {  //point
			return this->cell->get(ori, axes);
		}
		cvt p(Direction dir, Axes axes) const {
			if( Dim == 2 && axes == _Z_) {
				return 0.0;
			}
			return this->cell->get(ToOrientation(dir, axes), axes);
		}
		/*
		 *  overload the function of data
		 */
		ref_vt cd(st i) { //center data
			ASSERT(this->data != nullptr);
			return this->data->center(i);
		}
		const_ref_vt cd(st i) const { //center data
			ASSERT(this->data != nullptr);
			return this->data->center(i);
		}
		/*
		 *  show
		 */
		void show(const std::string& name = "" ) const {
			std::cout << "Node --- "<<name<<"\n";
			std::cout << "Dimension = " << Dim << "D\n";
			std::cout << "level      :" << this->_level << "\n";
			std::cout << "node type  :" << this->_node_type << "\n";
			std::cout << "idx        :" << this->_idx << "\n";
			std::cout << "CELL  show =========\n";
			std::cout << std::scientific;
			std::cout << "       min    " <<"     max    "<<"     d    "<< "     c    \n";
			std::cout << "x:" << p(_M_, _X_)<<" "<< p(_P_, _X_) <<" "<<d(_X_)<<" "<<cp(_X_)<< "\n";
			std::cout << "y:" << p(_M_, _Y_)<<" "<< p(_P_, _Y_) <<" "<<d(_Y_)<<" "<<cp(_Y_)<< "\n";
			if (Dim == 3) {
				std::cout << "z:" << p(_M_, _Z_)<<" "<< p(_P_, _Z_) <<" "<<d(_Z_)<<" "<<cp(_Z_)<< "\n";
			}
			std::cout << "DATA  show =========\n";
			if (nullptr == this->data) {
				std::cout << "No data \n";
			} else {
				//this->data->show_info();
			}
		}

	}
	;

	/*
	 *  functions out of class ================================================
	 */
	_TEMPLATE_COOV_V_DIM_ int GetDataIdx(const Node_<COO_VALUE, VALUE, DIM> *pn) {
		ASSERT(pn != nullptr);
		return pn->data->get_idx();
	}

	_TEMPLATE_COOV_V_DIM_
	Node_<COO_VALUE, VALUE
, DIM> *
GetpNodeAt(Node_<COO_VALUE, VALUE, DIM> *pn,
		const COO_VALUE &x,
		const COO_VALUE &y = 0,
		const COO_VALUE &z = 0) {
	typedef Node_<COO_VALUE, VALUE, DIM> *pNode;
	ASSERT(pn != nullptr);
	if (!pn->cell->is_in_on(x, y, z)) {
		return nullptr;
	}
	std::function<void(bool[], pNode, Float[])> condition =
	[](bool arr[], pNode pn, Float point[]) {
		st idx = pn->which_child(point[0],
				(pn->Dim >= 2) ? point[1] : 0,
				(pn->Dim == 3) ? point[2] : 0);
		arr[idx] = true;
	};
	pNode ret = nullptr;
	std::function<void(pNode, int)> fun =
	[&ret](pNode pn, int dummy) {
		if (pn->is_leaf()) {
			ret = pn;
		}
	};
	Float point[pn->Dim];
	point[0] = x;
	if (pn->Dim >= 2) {point[1] = y;}
	if (pn->Dim == 3) {point[2] = z;}
	int dummy = 0;
	pn->traversal_conditional(condition, point, fun, dummy);
	return ret;
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetpNodeSiblingPlus(const Node_<COO_VALUE, VALUE, DIM> *p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* f = p->father;
	if (f == nullptr) {
		return p;
	}
	for (int i = p->get_idx() + 1; i < Node::NumChildren; ++i) {
		Node* c = f->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return GetpNodeSiblingPlus(f);
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetpNodeSiblingPlus(Node_<COO_VALUE, VALUE, DIM> *p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* f = p->father;
	if (f == nullptr) {
		return p;
	}
	for (int i = p->get_idx() + 1; i < Node::NumChildren; ++i) {
		Node* c = f->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return GetpNodeSiblingPlus(f);
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetFirstChild(Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = nullptr;
	if(p==nullptr) {
		return c;
	}
	for (int i = 0; i < Node::NumChildren; ++i) {
		c = p->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return c;
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetFirstChild(const Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = nullptr;
	if(p==nullptr) {
		return c;
	}
	for (int i = 0; i < Node::NumChildren; ++i) {
		c = p->child[i];
		if (c != nullptr) {
			return c;
		}
	}
	return c;
}

_TEMPLATE_COOV_V_DIM_ Node_<COO_VALUE, VALUE, DIM> *
GetFirstLeaf(Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	Node* c = GetFirstChild(p);
	if (c == nullptr) {
		return p;
	} else {
		Node* resc = c;
		while (c != nullptr) {
			resc = c;
			c = GetFirstChild(c);
		}
		return resc;
	}
}

_TEMPLATE_COOV_V_DIM_ const Node_<COO_VALUE, VALUE, DIM> *
GetFirstLeaf(const Node_<COO_VALUE, VALUE, DIM> * p) {
	typedef Node_<COO_VALUE, VALUE, DIM> Node;
	const Node* c = GetFirstChild(p);
	if (c == nullptr) {
		return p;
	} else {
		const Node* resc = c;
		while (c != nullptr) {
			resc = c;
			c = GetFirstChild(c);
		}
		return resc;
	}
}

}
 //

#endif /* NODE_H_ */
