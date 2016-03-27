/************************
 //  \file   test_AABBox.h
 //  \brief
 // 
 //  \author czhou
 //  \date   17 juin 2015 
 ***********************/
#ifndef TEST_AABBOX_H_
#define TEST_AABBOX_H_
#include <iostream>
#include "ts_AABBox.h"
#include "ts_BBTree.h"
#include "ts_intersect.h"
#include "ts_triangle.h"
#include "ts_io.h"
#include "ts_surface.h"
#include "ts_surface_constructor.h"

#include <cstdlib>
#include <iostream>
#include <ctime>
#include <assert.h>

using namespace std;
namespace LarusTS {
void test_aabbox1() {
	Surface<Float, 3> Sur("head.gts");
	Sur.output_vtk("out.vtk");
	cout << "num faces " << Sur.faces.size() << "\n";
	cout << "++++++++++--------------------------------\n";
	Set<AABBox<Float, 3> > set_box;
	int i = 0;
	for (auto iter = Sur.faces.begin(); iter != Sur.faces.end(); ++iter) {
		auto ptr = (*iter);
		AABBox<Float, 3> tbox(ptr);
		set_box.insert(tbox);
		i++;
	}
	cout << "num boxes " << set_box.size() << "\n";

	BBTree<AABBox<Float, 3> > bbtree(set_box);
	//bbtree.output_vtk_height("h.vtk", 1);
	//cout << "tree size " << bbtree.size() << "\n";
	//vtk_show(bbtree,bbtree.height()-1);
	//set_box.begin()->output_vtk("begin.vtk");
	//std::list<vtkSmartPointer<vtkProp> > actors;
	//actors.push_back(vtk_new_actor(bbtree.begin()->box));
	//actors.push_back(vtk_new_actor_axes(0, 0, 0));
	//vtk_show_actor(actors);
	std::cout<<"is Orientable  "<<Sur.is_orientable()<<endl;
	Surface<Float, 3> Sur_t("triangle.gts");
	//Sur_t.output_vtk("out_tri.vtk");
	auto iter = Sur_t.faces.begin();
	AABBox<Float, 3> bs((*iter));
	//bs.output_vtk("bs.vtk");

	cout << "inter " << do_intersect_box_box(&bbtree, &bs) << endl;
	List<AABBox<Float, 3>*> lres;
	do_intersect_box_obj(&bbtree, &bs, lres);
	std::list<vtkSmartPointer<vtkProp> > actors;
	actors.push_back(vtk_new_actor(lres));
	actors.push_back(vtk_new_actor(Sur));
	//actors.push_back(vtk_new_actor_normal(Sur));

	actors.push_back(vtk_new_actor_axes(0, 0, 0));
	vtk_show_actor(actors);

	//output_vtk("inter.vtk", lres);
	//cout << "end test ===============\n";

}


void test_boxtri() {
	// tri and its box -------------------------
	Surface<Float, 3> Sur_t("triangle.gts");
	Sur_t.output_vtk("out_tri.vtk");
	auto iter = Sur_t.faces.begin();
	AABBox<Float, 3> bt((*iter));
	bt.output_vtk("btri.vtk");

	// box set ---------------------------------
	AABBox<Float, 3> bs(EMPTY, nullptr, 0, 0, 0, 1, 2, 3);
	bs.output_vtk("bs.vtk");

	cout << "interset :" << bs.are_overlapping((*iter)) << endl;
	cout << "end test ===============\n";
}

void test_construct() {
	std::cout<<"test construct ====="<<endl;
	Surface<Float, 3> Sur;
	Float r = 1;
	uInt n = 10;
	ConstructCircle(Sur, n, r);

	std::list<vtkSmartPointer<vtkProp> > actors;
	actors.push_back(vtk_new_actor(Sur));
	actors.push_back(vtk_new_actor_normal(Sur));
	actors.push_back(vtk_new_actor_axes(0, 0, 0));
	vtk_show_actor(actors);
}

void test_triangle() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator
	for (int c = 0; c < 1; c++) {
		cout << c << " !!!--------------------------- \n";
		for (int i = 0; i < 3; i++) {
			v0[i] = std::rand() / float(RAND_MAX);
			v1[i] = std::rand() / float(RAND_MAX);
			v2[i] = std::rand() / float(RAND_MAX);
			u0[i] = std::rand() / float(RAND_MAX);
			u1[i] = std::rand() / float(RAND_MAX);
			u2[i] = std::rand() / float(RAND_MAX);
		}

		int res1, res2;
		res1 = TriTriIsect_Guigue(v0, v1, v2, u0, u1, u2);
		res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
		//return code change
		if (res1 == -1 && res2 == 0) {
			cout << "pass \n";
		} else if (res1 == 1 && res2 == 1) {
			cout << "pass \n";
		} else {
			output_vtk("t1.vtk", v0, v1, v2);
			output_vtk("t2.vtk", u0, u1, u2);
			cout << res1 << "  " << res2 << endl;
			exit(0);
		}
	}
}

void test_triangle_cal() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator
	for (int c = 0; c < 1; c++) {
		cout << c << " !!!--------------------------- \n";
		for (int i = 0; i < 3; i++) {
			v0[i] = std::rand() / float(RAND_MAX);
			v1[i] = std::rand() / float(RAND_MAX);
			v2[i] = std::rand() / float(RAND_MAX);
			u0[i] = std::rand() / float(RAND_MAX);
			u1[i] = std::rand() / float(RAND_MAX);
			u2[i] = std::rand() / float(RAND_MAX);
		}

		int res1, res2;
		Float* x = NULL;
		Float* y = NULL;
		Float* z = NULL;
		short len = 0;
		res1 = TriTriIsect_Guigue_calculation(v0, v1, v2, u0, u1, u2, x, y, z,
				len);
		cout << "Len " << len << endl;
		Float resv0[] = { x[0], y[0], z[0] };
		Float resv1[] = { x[1], y[1], z[1] };
		output_vtk("seg.vtk", resv0, resv1);

		res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
		//return code change
		if (res1 == -1 && res2 == 0) {
			cout << "pass \n";
		} else if (res1 == 1 && res2 == 1) {
			cout << "pass \n";
		} else {
			cout << res1 << "  " << res2 << endl;
			exit(0);
		}
		output_vtk("t1.vtk", v0, v1, v2);
		output_vtk("t2.vtk", u0, u1, u2);
	}
}

void test_triangle_sigle() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float u0[3];
	Float u1[3];
	Float u2[3];
	std::srand(std::time(0)); // use current time as seed for random generator

	cout << " !!!--------------------------- \n";
	v0[0] = 0;
	v0[1] = 0;
	v0[2] = 0;
	v1[0] = 1;
	v1[1] = 0;
	v1[2] = 0;
	v2[0] = 0;
	v2[1] = 1;
	v2[2] = 0;

	u0[0] = 0.25;
	u0[1] = 0.25;
	u0[2] = 1;
	u1[0] = 0.5;
	u1[1] = 0;
	u1[2] = -1;
	u2[0] = 0;
	u2[1] = 0.5;
	u2[2] = -1;

	int res1, res2;
	Float* x = NULL;
	Float* y = NULL;
	Float* z = NULL;
	short len = 0;
	res1 = TriTriIsect_Guigue_calculation(v0, v1, v2, u0, u1, u2, x, y, z, len);
	cout << "Len " << len << endl;
	Float resv0[] = { x[0], y[0], z[0] };
	Float resv1[] = { x[1], y[1], z[1] };
	output_vtk("seg.vtk", resv0, resv1);
	res2 = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2);
	//return code change
	if (res1 == -1 && res2 == 0) {
		cout << "pass \n";
	} else if (res1 == 1 && res2 == 1) {
		cout << "pass \n";
	} else {
		//	output_vtk("t1.vtk", v0, v1, v2);
		//	output_vtk("t2.vtk", u0, u1, u2);
		cout << res1 << "  " << res2 << endl;
		exit(0);
	}
	output_vtk("t1.vtk", v0, v1, v2);
	output_vtk("t2.vtk", u0, u1, u2);

}

void test_orientation() {
	Float v0[3];
	Float v1[3];
	Float v2[3];

	Float d[3];

	cout << " !!!--------------------------- \n";
	v0[0] = 0;
	v0[1] = 0;
	v0[2] = 0;
	v1[0] = 1;
	v1[1] = 0;
	v1[2] = 0;
	v2[0] = 0;
	v2[1] = 1;
	v2[2] = 0;

	d[0] = 0;
	d[1] = 0;
	d[2] = 1;

	Float* a = v2;
	Float* b = v0;
	Float* c = v1;
	cout << SIGN3(ORIENT3D(d, b, c, a)) << endl;

}

}

#endif /* TEST_AABBOX_H_ */
