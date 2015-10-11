/************************
 //  \file   ts_surface.h
 //  \brief
 // 
 //  \author czhou
 //  \date   15 juin 2015 
 ***********************/
#ifndef TS_SURFACE_H_
#define TS_SURFACE_H_

#include "ts_define.h"
#include "ts_point.h"
#include "ts_vertex.h"
#include "ts_segment.h"
#include "ts_edge.h"
#include "ts_face.h"
#include <fstream>
#include <sstream>

namespace LarusTS
{

template<class TYPE, st DIM>
class Surface
{
public:
	typedef Surface<TYPE, DIM> self_class;
	typedef st size_type;
	typedef Point<TYPE, DIM> Poi;
	typedef Point<TYPE, DIM>* pPoi;
	typedef Segment<TYPE, DIM> Seg;
	typedef Segment<TYPE, DIM>* pSeg;
	typedef Edge<TYPE, DIM> Edg;
	typedef Edge<TYPE, DIM>* pEdg;
	typedef Vertex<TYPE, DIM> Ver;
	typedef Vertex<TYPE, DIM>* pVer;
	typedef Triangle<TYPE, DIM> Tri;
	typedef Triangle<TYPE, DIM>* pTri;
	typedef Face<TYPE, DIM> Fac;
	typedef Face<TYPE, DIM>* pFac;
	typedef List<pSeg> list_pSeg;
	typedef List<pVer> list_pVer;
	typedef List<pTri> list_pTri;
	typedef List<pFac> list_pFac;
	static const st Dim;
public:
	Set<pFac> faces;
	Set<pVer> c_vertex;
	Set<pEdg> c_edge;
public:
	Surface()
	{
	}
	Surface(const String& filename);

	~Surface()
	{
		clear();
	}

	void clear(){
		for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
			delete (*iter);
		}
		c_vertex.clear();
		for (auto iter = c_edge.begin(); iter != c_edge.end(); ++iter) {
			delete (*iter);
		}
		c_edge.clear();
		for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
			delete (*iter);
		}
		faces.clear();

	}
	size_type size_edge() const
	{
		return c_edge.size();
	}
	size_type size_vertex() const
	{
		return c_vertex.size();
	}
	size_type size_face() const
	{
		return faces.size();
	}
	void transfer(TYPE dx, TYPE dy, TYPE dz = 0){
		for (auto iter = c_vertex.begin();iter != c_vertex.end(); ++iter) {
			Ver& ver = (*(*iter));
			ver[0]  = ver[0] + dx;
			ver[1]  = ver[1] + dy;
			if(DIM == 3){
				ver[2]  = ver[2] + dy;
			}
		}
	}
	// connection -------------------------------
	// check ------------------------------------

	// show -------------------------------------
	void show_vertex() const;
	void show_edge() const;
	void show_face() const;

	// output -----------------------------------
	void output_vtk(const String& fn) const;

};

template<class TYPE, st DIM>
Surface<TYPE, DIM>::Surface(const String& filename)
{
	std::ifstream fs;
	fs.open(filename.c_str(), std::fstream::in); //read
	if (!fs.is_open()) {
		std::cerr << "!> Can't find file. " << filename.c_str() << " \n";
		exit(-1);
	}
	uInt n = 0;
	uInt nv, ne, nf;
	Vector<pVer> v_vertex;
	Vector<Int> v_vc;
	Vector<pEdg> v_edge;
	Vector<Int> v_ec;
	Vector<pFac> v_face;
	while (!fs.eof()) {
		String sline;
		getline(fs, sline, '\n');
		if (sline.length() != 0) {
			if (n == 0) {
				std::istringstream istr(sline);
				istr >> nv >> ne >> nf;
			} else if (n <= nv) {
				std::istringstream istr(sline);
				TYPE x, y, z;
				istr >> x >> y >> z;
				pVer pver = new Ver(x, y, z);
				v_vertex.push_back(pver);
				v_vc.push_back(0);
			} else if (n <= nv + ne) {
				std::istringstream istr(sline);
				uInt i_v1, i_v2;
				istr >> i_v1 >> i_v2;
				pEdg pedg = new Edg(v_vertex[i_v1 - 1], v_vertex[i_v2 - 1]);
				v_vc[i_v1 - 1]++;
				v_vc[i_v2 - 1]++;
				v_edge.push_back(pedg);
				v_ec.push_back(0);
			} else if (n <= nv + ne + nf) {
				std::istringstream istr(sline);
				uInt i_e1, i_e2, i_e3;
				istr >> i_e1 >> i_e2 >> i_e3;
				pFac pfac = new Fac(v_edge[i_e1 - 1], v_edge[i_e2 - 1],
						v_edge[i_e3 - 1], this);
				v_ec[i_e1 - 1]++;
				v_ec[i_e2 - 1]++;
				v_ec[i_e3 - 1]++;
				v_face.push_back(pfac);
			}
		}
		n++;
	}
	for (uInt i = 0; i < v_vc.size(); ++i) {
		if (v_vc[i] == 0) {
			delete v_vertex[i];
		} else {
			c_vertex.insert(v_vertex[i]);
		}

	}
	for (uInt i = 0; i < v_ec.size(); ++i) {
		if (v_ec[i] == 0) {
			delete v_edge[i];
		} else {
			c_edge.insert(v_edge[i]);
		}
	}
	//insert pface to set
	for (auto iter = v_face.begin(); iter != v_face.end(); ++iter) {
		auto ptr = (*iter);
		faces.insert(ptr);
	}
	//for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
	//	auto ptr = (*iter);
	//	ptr->show();
	//}

}
template<class TYPE, st DIM>
void Surface<TYPE, DIM>::show_vertex() const
{
	std::cout << " size vertex = " << c_vertex.size() << "\n";
for (typename Set<pVer>::const_iterator iter = c_vertex.begin();
			iter != c_vertex.end(); ++iter) {
		(*iter)->show();
	}
}
template<class TYPE, st DIM>
void Surface<TYPE, DIM>::show_edge() const
{
	std::cout << " size edge = " << c_edge.size() << "\n";
	for (auto iter = c_edge.begin(); iter != c_edge.end(); ++iter) {
		(*iter)->show();
	}
}

template<class TYPE, st DIM>
void Surface<TYPE, DIM>::show_face() const
{
	std::cout << " size face = " << faces.size() << "\n";
	for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
		(*iter)->show();
	}
}

template<class TYPE, st DIM>
void Surface<TYPE, DIM>::output_vtk(const String& fn) const
{
	FILE* fptr = fopen(fn.c_str(), "w"); //write
	if (fptr == NULL) {
		std::cerr << "!> Open file error! " << fn << " \n";
		exit(-1);
	}
	fprintf(fptr, "# vtk DataFile Version 2.0\n"
			"Generated by LarusTS\n"
			"ASCII\n"
			"DATASET POLYDATA\n"
			"POINTS %lu float\n", c_vertex.size());
	Map<pVer, uInt> m_veridx;
	uInt count = 0;
	for (auto iter = c_vertex.begin(); iter != c_vertex.end(); ++iter) {
		auto pt = (*iter);
		fprintf(fptr, "%f %f %f \n", pt->x(), pt->y(), pt->z());
		m_veridx.insert(Pair<pVer, uInt>(pt, count));
		count++;
	}
	fprintf(fptr, "POLYGONS %lu %lu\n", faces.size(), faces.size() * 4);
	for (auto iter = faces.begin(); iter != faces.end(); ++iter) {
		auto pt = (*iter);
		uInt i1, i2, i3;
		i1 = m_veridx.find(pt->e1->v1)->second;
		i2 = m_veridx.find(pt->e1->v2)->second;
		if (pt->e2->v2 != pt->e1->v2 && pt->e2->v2 != pt->e1->v1) {
			i3 = m_veridx.find(pt->e2->v2)->second;
		} else {
			i3 = m_veridx.find(pt->e2->v1)->second;
		}
		fprintf(fptr, "3 %u %u %u\n", i1, i2, i3);
	}
	fclose(fptr);
}

}

#endif /* TS_SURFACE_H_ */
