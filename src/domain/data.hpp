#ifndef DATA_H_
#define DATA_H_

#include "domain_define.hpp"
#include "../algebra/array_list.hpp"

namespace carpio {

template<typename VALUE, size_t DIM>
class Data {
public:
	static const size_t Dim = DIM;
	static const size_t NumFaces = DIM + DIM;
	static const size_t NumVertexes = (DIM == 3) ? 8 : (DIM + DIM); //

	typedef VALUE value_t;
	typedef Data<VALUE, DIM> self;

	typedef void (*pfunction)(self*, utPointer);

protected:
	int _idx;
	ArrayListV<value_t> _center;
	ArrayListV<value_t> _face[NumFaces];
	ArrayListV<value_t> _vertex[NumVertexes];
	utPointer untype;
public:
	Data() {
		_idx = 0;
		untype = nullptr;
	}
	Data(const size_t& nc, const size_t& nf, const size_t& nv, utPointer utp) :
			_center(nc) {
		_idx = 0;
		for (int i = 0; i < NumFaces; ++i) {
			_face[i].reconstruct(nf);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			_face[i].reconstruct(nv);
		}
		untype = utp;
	}

	inline value_t& Center(size_t i){
		ASSERT(i<_center.size());
		return _center[i];
	}

	inline const value_t& Center(size_t i) const{
		ASSERT(i<_center.size());
		return _center[i];
	}

	inline value_t& Face(Direction d, size_t i){
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline const value_t& Face(Direction d, size_t i) const{
		ASSERT(i < this->NumFaces);
		return _face[i];
	}

	inline value_t& Vertex(Direction d, size_t i){
		ASSERT(i<this->NumVertexes);
		return _vertex[i];
	}

	inline const value_t& Vertex(Direction d, size_t i) const{
		ASSERT(i<this->NumVertexes);
		return _vertex[i];
	}


	bool Empty() const {
		bool res = true;
		res = res && (_center.size() == 0);
		for (int i = 0; i < NumFaces; ++i) {
			res = res && (_face[i].size() == 0);
		}
		for (int i = 0; i < NumVertexes; ++i) {
			res = res && (_vertex[i].size() == 0);
		}
		return res;
	}

	void ShowInfo() const {
		std::cout << "center data:" << this->_center.size() << "\n";
		std::cout << "face data   :" << "\n";
		for (int i = 0; i < NumFaces; ++i) {
			std::cout << "    " << i << "       :" << _face[i].size()
					<< "\n";
		}
		for (int i = 0; i < NumVertexes; ++i) {
			std::cout << "    " << i << "       :" << _vertex[i].size()
					<< "\n";
		}
	}

};
}


#endif /* DATA_H_ */
