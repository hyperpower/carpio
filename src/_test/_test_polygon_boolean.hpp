#ifndef _TEST_POLYGON_BOOLEAN_HPP_
#define _TEST_POLYGON_BOOLEAN_HPP_

#include "../geometry/geometry.hpp"
#include "../io/io_gnuplot_domain.h"
#include <iostream>
using namespace std;
namespace carpio {
void _case_1() {
	//Define two segments
	Segment_<Float, 2> sc(0.0, 1.0, 0.0, 1.0);
	Segment_<Float, 2> so(1.0, 0.5, 0.0, 0.5);
	//show segments
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Segment2D(ga, sc);
	lga.push_back(ga);
	GnuplotActor_Segment2D(ga, so);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//
	//int wt = IntersectWAType(sc, so);
	//ShowWAType(wt);
	//gp.set_label(1, ParseWAType(wt), 0.5, 0.1);
	gp.plot(lga);
	//
}
void _case_2() {
	//Define two segments
	Segment_<Float, 2> sc(0.0, 1.0, 0.0, 1.0);
	Segment_<Float, 2> so(0.5, 2.0, 0.5, 2.0);
	//show segments
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Segment2D(ga, sc);
	lga.push_back(ga);
	GnuplotActor_Segment2D(ga, so);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();
	//
	//int wt = IntersectWAType(sc, so);
	//ShowWAType(wt);
	//gp.set_label(1, ParseWAType(wt), 0.5, 0.1);
	gp.plot(lga);
	//
}
void test_wa_intersect() {
	_case_1();
	cout << "====fin test wa intersect=====\n";
}

void _polygon_intersect_case1() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon sub;
	x1 = 3;
	y1 = 3;
	x2 = 5;
	y2 = 5;
	CreatCube(clip, x1, y1, x2, y2);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	//
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}

void _polygon_intersect_case2() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(5);
	arrp[0].reconstruct(3, 0);
	arrp[1].reconstruct(5, 0);
	arrp[2].reconstruct(5, 5);
	arrp[3].reconstruct(4, 5);
	arrp[4].reconstruct(3, 4);
	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case3() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(6);
	arrp[0].reconstruct(3, 0);
	arrp[1].reconstruct(5, 0);
	arrp[2].reconstruct(5, 5);
	arrp[3].reconstruct(4, 5);
	arrp[4].reconstruct(3, 4);
	arrp[5].reconstruct(2, 1);
	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case4() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(6);
	arrp[0].reconstruct(3, 0);
	arrp[1].reconstruct(5, 0);
	arrp[2].reconstruct(5, 5);
	arrp[3].reconstruct(4.5, 5);
	arrp[4].reconstruct(4, 4);
	arrp[5].reconstruct(2, 1);
	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case5() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(3);
	arrp[0].reconstruct(4, 3);
	arrp[1].reconstruct(5, 0);
	arrp[2].reconstruct(5, 5);

	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case6() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(3);
	arrp[0].reconstruct(2, 2.5);
	arrp[1].reconstruct(4, 4);
	arrp[2].reconstruct(3, 2);

	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case7() {
	cout << "============test case 1 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 1, y1 = 1, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(3);
	arrp[0].reconstruct(4, 4);
	arrp[1].reconstruct(5, 0);
	arrp[2].reconstruct(5, 5);

	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);
	//Intersect(clip, sub, res);
}

void _polygon_intersect_case8() {
	cout << "============test case 8 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 2, y1 = 2, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(3);
	arrp[0].reconstruct(3, 5);
	arrp[1].reconstruct(5, 3);
	arrp[2].reconstruct(5, 5.1);

	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_title("case 8");
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case9() {
	cout << "============test case 9 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 2, y1 = 2, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(4);
	arrp[0].reconstruct(3, 1);
	arrp[1].reconstruct(5, 3);
	arrp[2].reconstruct(3, 5);
	arrp[3].reconstruct(1, 3);

	Polygon sub(arrp);
	//clip.reverse();
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon_vector(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon_vector(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_title("case 9");
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon_vector(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void _polygon_intersect_case10() {
	cout << "============test case 10 ===============\n";
	//new polygon
	Polygon clip;
	Float x1 = 2, y1 = 2, x2 = 4, y2 = 4;
	CreatCube(clip, x1, y1, x2, y2);
	Polygon::ArrP arrp(3);
	arrp[0].reconstruct(1, 3);
	arrp[1].reconstruct(3, 1);
	arrp[2].reconstruct(0, 0);

	Polygon sub(arrp);
	// show
	std::list<Gnuplot_actor> lga;
	Gnuplot_actor ga;
	GnuplotActor_Polygon(ga, clip);
	lga.push_back(ga);
	GnuplotActor_Polygon(ga, sub);
	lga.push_back(ga);
	Gnuplot gp;
	gp.set_title("case 10");
	gp.set_equal_ratio();

	//gp.set_xrange(-0.1, 1.1);
	//gp.set_yrange(-0.1, 1.1);
	std::list<Polygon> res;
	WAClipping_<Float> Clip(clip, sub);
	Clip.clipping(res);
	Clip.show();
	Clip.show(1);
	for (std::list<Polygon>::iterator iter = res.begin(); iter != res.end();
			++iter) {
		GnuplotActor_Polygon(ga, (*iter));
		lga.push_back(ga);
	}
	gp.plot(lga);

	//Intersect(clip, sub, res);
}
void test_list() {
	//test list
	std::list<int> mylist;
	std::list<int>::iterator it;

	// set some initial values:
	for (int i = 1; i <= 5; ++i)
		mylist.push_back(i); // 1 2 3 4 5

	it = mylist.begin();  //                                 1 2 3 4 5
	++it;                 // it points now to number 2         ^
	ASSERT(*it == 2);
	it = mylist.insert(it, 10);   // 1 10 2 3 4 5
	ASSERT(*it == 10);
	ASSERT(*(++it) == 2);

	it = mylist.end();
	++it;
	ASSERT(*it == 1);
	cout << "list test pass ==============\n";
}

//void test_WAType() {
//	cout << "WAType ===================\n";
//	// Define two segments
//	// colinear
//	Segment_<Float, 2> sc(0.0, 1.0, 0.0, 1.0);
//	Segment_<Float, 2> so(0.5, 2.0, 0.5, 2.0);
//
//	//
//	int wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ASSERT(wt == WA_ERR);
//
//	cout << "pass 1  =========================\n";
//	//
//	so.reconstruct(1.0, 1.0, 0.0, 1.0);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	int ass;
//	ass = WA_IN | WA_OUT_R | WA_HEAD_C | WA_HEAD_S | WA_CLIP | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 2  =========================\n";
//	so.reconstruct(1.0, 0.5, 0.0, 0.5);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ass = WA_IN | WA_OUT_R | WA_HEAD_S | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 3  =========================\n";
//	so.reconstruct(0.5, 0.0, 0.5, 1.0);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ass = WA_IN | WA_OUT_R | WA_TAIL_S | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 4  =========================\n";
//	so.reconstruct(0.0, 0.5, 1.0, 0.5);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ass = WA_OUT | WA_IN_R | WA_HEAD_S | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 5  =========================\n";
//	so.reconstruct(0.0, 1.0, 1.0, 1.0);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ass = WA_OUT | WA_IN_R | WA_HEAD_C | WA_HEAD_S | WA_CLIP | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 6  =========================\n";
//	so.reconstruct(0.0, 0.0, 1.0, 0.0);
//	wt = IntersectWAType(sc, so);
//	ShowWAType(wt);
//	ass = WA_OUT | WA_IN_R | WA_TAIL_C | WA_HEAD_S | WA_CLIP | WA_SUBJ;
//	ASSERT(ass == wt);
//	cout << "pass 7  =========================\n";
//	//show segments
//	std::list<Gnuplot_actor> lga;
//	Gnuplot_actor ga;
//	GnuplotActor_Segment2D(ga, sc);
//	lga.push_back(ga);
//	GnuplotActor_Segment2D(ga, so);
//	lga.push_back(ga);
//	Gnuplot gp;
//	gp.set_equal_ratio();
//	gp.set_xrange(-0.1, 1.1);
//	gp.set_yrange(-0.1, 1.1);
//	gp.set_label(1, ParseWAType(wt), 0.5, 0.1);
//	gp.plot(lga);
//}

void test_polygon_intersect() {
	//_polygon_intersect_case1();
	//test_WAType();
	_polygon_intersect_case9();
	cout << "====fin test polygon intersect=====\n";

}
}
#endif
