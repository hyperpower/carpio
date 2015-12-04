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
	typedef const Grid_<COO_VALUE, VALUE, DIM> * const_pSelf;
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
							new Node(nullptr, // father
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
	/*
	 *  index as 2 or 3 demension
	 */
	reference operator()(size_type i, size_type j = 0, size_type k = 0) {
		return nodes(i, j, k);
	}

	const_reference operator()(size_type i, size_type j = 0,
			size_type k = 0) const {
		return nodes(i, j, k);
	}
	/*
	 *  index as 1 demension
	 */
	reference at_1d(size_type i) {
		return nodes.at_1d(i);
	}

	const_reference at_1d(size_type i) const {
		return nodes.at_1d(i);
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

	const pNode get_last_root_pNode() const {
		for (st i = nodes.size() - 1; i >= 0; --i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}
	pNode get_last_root_pNode() {
		for (st i = nodes.size() - 1; i >= 0; --i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}
	const pNode get_first_root_pNode() const {
		for (st i = 0; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}
	pNode get_first_root_pNode() {
		for (st i = 0; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return nodes.at_1d(i);
			}
		}
		return nullptr;
	}

	st get_first_root_pNode_idx1d() const {
		st i = 0;
		for (; i < nodes.size(); ++i) {
			if (nodes.at_1d(i) != nullptr) {
				return i;
			}
		}
		return i;
	}

	st get_last_root_pNode_idx1d() const {
		size_type i = nodes.size() - 1;
		for (; i >= 0; --i) {
			if (nodes.at_1d(i) != nullptr) {
				return i;
			}
		}
		return i;
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
						xm = nodes.check_idx_ijk(i - 1, j, k) ?
								nodes(i - 1, j, k) : nullptr;
						xp = nodes.check_idx_ijk(i + 1, j, k) ?
								nodes(i + 1, j, k) : nullptr;
						ym = nodes.check_idx_ijk(i, j - 1, k) ?
								nodes(i, j - 1, k) : nullptr;
						yp = nodes.check_idx_ijk(i, j + 1, k) ?
								nodes(i, j + 1, k) : nullptr;
						zm = nodes.check_idx_ijk(i, j, k - 1) ?
								nodes(i, j, k - 1) : nullptr;
						zp = nodes.check_idx_ijk(i, j, k + 1) ?
								nodes(i, j, k + 1) : nullptr;
						cnode->set_neighbor(xm, xp, ym, yp, zm, zp);
					}
				}
			}
		}
	}
	/*
	 *  iterator leaf node
	 */
protected:
	template<typename COV, typename V, int D, class _Ref, class _Ptr>
	class iterator_leaf_ {
	public:
		typedef COV cvt;
		typedef VALUE vt;
		typedef Node_<COV, V, D> Node;
		typedef Node_<COV, V, D> *pNode;
		typedef const Node_<COV, V, D> const_Node;
		typedef const Node_<COV, V, D>* const_pNode;

		typedef iterator_leaf_<COV, V, D, Node&, pNode> iterator;
		typedef iterator_leaf_<COV, V, D, const_Node&, const_pNode> const_iterator;
		typedef iterator_leaf_<COV, V, D, _Ref, _Ptr> Self;

		typedef Grid_<COV, V, D> Grid;
		typedef Grid_<COV, V, D>* pGrid;
		typedef const Grid_<COV, V, D>* const_pGrid;

		typedef Node value_type;
		typedef _Ptr pointer;
		typedef _Ref reference;

		const_pGrid _f;
		st _idx;
		_Ptr _ptr;

		iterator_leaf_() {
			_ptr = nullptr;
			_f = nullptr;
			_idx = 0;
		}
		iterator_leaf_(const_pGrid f, st idx, _Ptr ptr) :
				_f(f), _idx(idx), _ptr(ptr) {
		}
		iterator_leaf_(const iterator& _x) :
			_f(_x._f), _idx(_x._idx), _ptr(_x._ptr){
		}
	protected:
		pNode _incr_root() {
			st count = 0;
			for (st ii = (_idx + 1) % _f->size(); count < _f->size();
					++ii, ++count) {
				pNode pt = _f->at_1d(ii);
				if (nullptr != pt) {
					_idx = ii;
					return pt;
				}
			}
			return nullptr;
		}
		void _incr() {
			pNode end = _f->get_last_root_pNode();
			if (_ptr == end) {
				return;  //will not increase
			}
			_Ptr s = GetpNodeSiblingPlus(_ptr);
			if (s->father != nullptr) {
				_ptr = GetFirstLeaf(s);
			} else {
				if (s == end) {
					_ptr = s;
				} else {
					pNode root_pt = _incr_root(); //this will change _idx
					_ptr = GetFirstLeaf(root_pt);
				}
			}
		}
	public:
		bool operator==(const iterator_leaf_& _x) const {
			return _ptr == _x._ptr;
		}
		bool operator!=(const iterator_leaf_& _x) const {
			return _ptr != _x._ptr;
		}

		reference operator*() const {
			return (*_ptr);
		}

		pointer operator->() const {
			return &(operator*());
		}

		Self & operator++() {
			this->_incr();
			return *this;
		}

		Self operator++(int) {
			Self __tmp = *this;
			this->_incr();
			return __tmp;
		}

		bool is_exist() {
			return _ptr != nullptr;
		}

		pointer get_pointer() {
			return _ptr;
		}

		const pointer get_pointer() const {
			return _ptr;
		}
	};

	/*
	 *  typedef of iterator leaf
	 */
public:
	typedef iterator_leaf_<cvt, vt, Dim, Node&, Node*> iterator_leaf;
	typedef iterator_leaf_<cvt, vt, Dim, const Node&, const Node*> const_iterator_leaf;

	iterator_leaf begin_leaf() {
		size_type idx = this->get_first_root_pNode_idx1d();
		pNode pt = this->get_first_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return iterator_leaf(this, idx, pn);
	}
	const_iterator_leaf begin_leaf() const {
		size_type idx = this->get_first_root_pNode_idx1d();
		const Node* pt = this->get_first_root_pNode();
		const Node* pn = GetFirstLeaf(pt);
		const_pSelf pg = this;
		return const_iterator_leaf(pg, idx, pn);
	}
	iterator_leaf last_leaf() {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_last_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return iterator_leaf(this, idx, pn);
	}
	const_iterator_leaf last_leaf() const {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_first_root_pNode();
		pNode pn = GetFirstLeaf(pt);
		return const_iterator_leaf(this, idx, pn);
	}
	iterator_leaf end_leaf() {
		size_type idx = this->get_last_root_pNode_idx1d();
		pNode pt = this->get_last_root_pNode();
		return iterator_leaf(this, idx, pt);
	}
	const_iterator_leaf end_leaf() const {
		size_type idx = this->get_last_root_pNode_idx1d();
		const pNode pt = this->get_last_root_pNode();
		return const_iterator_leaf(this, idx, pt);
	}

};

/*
 *  Function out of class =================================================
 */

template<typename COO_VALUE, typename VALUE, int DIM>
Node_<COO_VALUE, VALUE, DIM> *
GetpNodeAt(Grid_<COO_VALUE, VALUE, DIM> *pg, const COO_VALUE &x,
		const COO_VALUE &y = 0, const COO_VALUE &z = 0) {
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
GetpNodeAt(const Grid_<COO_VALUE, VALUE, DIM> *pg, const COO_VALUE &x,
		const COO_VALUE &y = 0, const COO_VALUE &z = 0) {

}

}
#endif
