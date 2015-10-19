/************************
 //  \file   ts_tri_moller.h
 //  \brief
 // 
 //  \author czhou
 //  \date   19 juin 2015 
 ***********************/
#ifndef TS_TRI_MOLLER_H_
#define TS_TRI_MOLLER_H_

#include "ts_define.h"

namespace LarusTS
{
/********************************************************/
/* AABB-triangle overlap test code                      */
/* by Tomas Akenine-Möller                              */
/* Function: int triBoxOverlap(float boxcenter[3],      */
/*          float boxhalfsize[3],float triverts[3][3]); */
/* History:                                             */
/*   2001-03-05: released the code in its first version */
/*   2001-06-18: changed the order of the tests, faster */
/*                                                      */
/* Acknowledgement: Many thanks to Pierre Terdiman for  */
/* suggestions and discussions on how to optimize code. */
/* Thanks to David Hunt for finding a ">="-bug!         */
/********************************************************/
#include <math.h>
#include <stdio.h>

static const uInt X = 0;
static const uInt Y = 1;
static const uInt Z = 2;

#define CROSS(dest,v1,v2)                  \
          dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
          dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
          dest[2]=v1[0]*v2[1]-v1[1]*v2[0];

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2)        \
          dest[0]=v1[0]-v2[0]; \
          dest[1]=v1[1]-v2[1]; \
          dest[2]=v1[2]-v2[2];

#define FINDMINMAX(x0,x1,x2,min,max) \
          min = max = x0;            \
          if(x1<min) min=x1;         \
          if(x1>max) max=x1;         \
          if(x2<min) min=x2;         \
          if(x2>max) max=x2;

template<class VALUE>
int planeBoxOverlap(VALUE normal[3], VALUE vert[3], VALUE maxbox[3]) // -NJMP-
{
	int q;
	VALUE vmin[3], vmax[3], v;
	for (q = X; q <= Z; q++) {
		v = vert[q];					// -NJMP-
		if (normal[q] > 0.0) {
			vmin[q] = -maxbox[q] - v;	// -NJMP-
			vmax[q] = maxbox[q] - v;	// -NJMP-
		} else {
			vmin[q] = maxbox[q] - v;	// -NJMP-
			vmax[q] = -maxbox[q] - v;	// -NJMP-
		}
	}
	if (DOT(normal,vmin) > 0.0)
		return 0;	// -NJMP-
	if (DOT(normal,vmax) >= 0.0)
		return 1;	// -NJMP-

	return 0;
}

/*======================== X-tests ========================*/
#define AXISTEST_X01(a, b, fa, fb)			               \
	p0 = a*v0[Y] - b*v0[Z];			       	               \
	p2 = a*v2[Y] - b*v2[Z];			       	               \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];       \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_X2(a, b, fa, fb)			               \
	p0 = a*v0[Y] - b*v0[Z];			                       \
	p1 = a*v1[Y] - b*v1[Z];			       	               \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[Y] + fb * boxhalfsize[Z];       \
	if(min>rad || max<-rad) return 0;

/*======================== Y-tests ========================*/
#define AXISTEST_Y02(a, b, fa, fb)			               \
	p0 = -a*v0[X] + b*v0[Z];		      	               \
	p2 = -a*v2[X] + b*v2[Z];	       	       	           \
        if(p0<p2) {min=p0; max=p2;} else {min=p2; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];       \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_Y1(a, b, fa, fb)			               \
	p0 = -a*v0[X] + b*v0[Z];		      	               \
	p1 = -a*v1[X] + b*v1[Z];	     	       	           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Z];       \
	if(min>rad || max<-rad) return 0;

/*======================== Z-tests ========================*/

#define AXISTEST_Z12(a, b, fa, fb)			   \
	p1 = a*v1[X] - b*v1[Y];			           \
	p2 = a*v2[X] - b*v2[Y];			       	   \
        if(p2<p1) {min=p2; max=p1;} else {min=p1; max=p2;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;

#define AXISTEST_Z0(a, b, fa, fb)			   \
	p0 = a*v0[X] - b*v0[Y];				       \
	p1 = a*v1[X] - b*v1[Y];			           \
        if(p0<p1) {min=p0; max=p1;} else {min=p1; max=p0;} \
	rad = fa * boxhalfsize[X] + fb * boxhalfsize[Y];   \
	if(min>rad || max<-rad) return 0;

template<class VALUE>
int triBoxOverlap(VALUE boxcenter[3],    //
		VALUE boxhalfsize[3],  //
		VALUE triverts[3][3])  //
{

	/*    use separating axis theorem to test overlap between triangle and box */
	/*    need to test for overlap in these directions: */
	/*    1) the {x,y,z}-directions (actually, since we use the AABB of the triangle */
	/*       we do not even need to test these) */
	/*    2) normal of the triangle */
	/*    3) crossproduct(edge from tri, {x,y,z}-directin) */
	/*       this gives 3x3=9 more tests */
	VALUE v0[3], v1[3], v2[3];
//   double axis[3];
	VALUE min, max, p0, p1, p2, rad, fex, fey, fez; // -NJMP- "d" local variable removed
	VALUE normal[3], e0[3], e1[3], e2[3];

	/* This is the fastest branch on Sun */
	/* move everything so that the boxcenter is in (0,0,0) */
	SUB(v0, triverts[0], boxcenter);
	SUB(v1, triverts[1], boxcenter);
	SUB(v2, triverts[2], boxcenter);

	/* compute triangle edges */
	SUB(e0, v1, v0); /* tri edge 0 */
	SUB(e1, v2, v1); /* tri edge 1 */
	SUB(e2, v0, v2); /* tri edge 2 */

	/* Bullet 3:  */
	/*  test the 9 tests first (this was faster) */
	fex = ABS(e0[X]);
	fey = ABS(e0[Y]);
	fez = ABS(e0[Z]);
	AXISTEST_X01(e0[Z], e0[Y], fez, fey);
	AXISTEST_Y02(e0[Z], e0[X], fez, fex);
	AXISTEST_Z12(e0[Y], e0[X], fey, fex);

	fex = ABS(e1[X]);
	fey = ABS(e1[Y]);
	fez = ABS(e1[Z]);
	AXISTEST_X01(e1[Z], e1[Y], fez, fey);
	AXISTEST_Y02(e1[Z], e1[X], fez, fex);
	AXISTEST_Z0(e1[Y], e1[X], fey, fex);

	fex = ABS(e2[X]);
	fey = ABS(e2[Y]);
	fez = ABS(e2[Z]);
	AXISTEST_X2(e2[Z], e2[Y], fez, fey);
	AXISTEST_Y1(e2[Z], e2[X], fez, fex);
	AXISTEST_Z12(e2[Y], e2[X], fey, fex);

	/* Bullet 1: */
	/*  first test overlap in the {x,y,z}-directions */
	/*  find min, max of the triangle each direction, and test for overlap in */
	/*  that direction -- this is equivalent to testing a minimal AABB around */
	/*  the triangle against the AABB */

	/* test in X-direction */
	FINDMINMAX(v0[X], v1[X], v2[X], min, max);
	if (min > boxhalfsize[X] || max < -boxhalfsize[X])
		return 0;

	/* test in Y-direction */
	FINDMINMAX(v0[Y], v1[Y], v2[Y], min, max);
	if (min > boxhalfsize[Y] || max < -boxhalfsize[Y])
		return 0;

	/* test in Z-direction */
	FINDMINMAX(v0[Z], v1[Z], v2[Z], min, max);
	if (min > boxhalfsize[Z] || max < -boxhalfsize[Z])
		return 0;

	/* Bullet 2: */
	/*  test if the box intersects the plane of the triangle */
	/*  compute plane equation of triangle: normal*x+d=0 */
	CROSS(normal, e0, e1);
	// -NJMP- (line removed here)
	if (!planeBoxOverlap(normal, v0, boxhalfsize))
		return 0;	// -NJMP-

	return 1; /* box and triangle overlaps */
}

/********************************************************/
/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * Updated June 1999: removed the divisions -- a little faster now!
 * Updated October 1999: added {} to CROSS and SUB macros
 *
 * int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
 *                      float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 *********************************************************/

#define FABS(x) (float(fabs(x)))        /* implement as is fastest on your machine */

/* if USE_EPSILON_TEST is true then we do a check:
 if |dv|<EPSILON then dv=0.0;
 else no check is done (which is less robust)
 */
#define USE_EPSILON_TEST FALSE
#define EPSILON 0.000001

/* sort so that a<=b */
#define SORT(a,b)       \
             if(a>b)    \
             {          \
               VALUE c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }

/* this edge to edge test is based on Franlin Antonio's gem:
 "Faster Line Segment Intersection", in Graphics Gems III,
 pp. 199-202 */
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 1;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 1;                      \
    }                                                 \
  }

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  float Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  float a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 1;                 \
  }                                         \
}

template<class VALUE>
int coplanar_tri_tri(VALUE N[3], VALUE V0[3], VALUE V1[3], VALUE V2[3],
		VALUE U0[3], VALUE U1[3], VALUE U2[3])
{
	VALUE A[3];
	short i0, i1;
	/* first project onto an axis-aligned plane, that maximizes the area */
	/* of the triangles, compute indices: i0,i1. */
	A[0] = FABS(N[0]);
	A[1] = FABS(N[1]);
	A[2] = FABS(N[2]);
	if (A[0] > A[1]) {
		if (A[0] > A[2]) {
			i0 = 1; /* A[0] is greatest */
			i1 = 2;
		} else {
			i0 = 0; /* A[2] is greatest */
			i1 = 1;
		}
	} else /* A[0]<=A[1] */
	{
		if (A[2] > A[1]) {
			i0 = 0; /* A[2] is greatest */
			i1 = 1;
		} else {
			i0 = 0; /* A[1] is greatest */
			i1 = 2;
		}
	}

	/* test all edges of triangle 1 against the edges of triangle 2 */
	EDGE_AGAINST_TRI_EDGES(V0, V1, U0, U1, U2);
	EDGE_AGAINST_TRI_EDGES(V1, V2, U0, U1, U2);
	EDGE_AGAINST_TRI_EDGES(V2, V0, U0, U1, U2);

	/* finally, test if tri1 is totally contained in tri2 or vice versa */
	POINT_IN_TRI(V0, U0, U1, U2);
	POINT_IN_TRI(U0, V0, V1, V2);

	return 0;
}

#define NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1) \
{ \
        if(D0D1>0.0f) \
        { \
            /* here we know that D0D2<=0.0 */ \
            /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
            A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else if(D0D2>0.0f)\
        { \
            /* here we know that d0d1<=0.0 */ \
			/* that is D0, D2 are on the same side, D1 on the other or on the plane */ \
            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D1*D2>0.0f || D0!=0.0f) \
        { \
             /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
            A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
        } \
        else if(D1!=0.0f) \
        { \
            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
        } \
        else if(D2!=0.0f) \
        { \
            A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
        } \
        else \
        { \
                /* triangles are coplanar */ \
            return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
        } \
}

template<class VALUE>
int NoDivTriTriIsect( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //
		VALUE U0[3], VALUE U1[3], VALUE U2[3]) //
{
	VALUE E1[3], E2[3];
	VALUE N1[3], N2[3], d1, d2;
	VALUE du0, du1, du2, dv0, dv1, dv2;
	VALUE D[3];
	VALUE isect1[2], isect2[2];
	VALUE du0du1, du0du2, dv0dv1, dv0dv2;
	short index;
	VALUE vp0, vp1, vp2;
	VALUE up0, up1, up2;
	VALUE bb, cc, max;

	/* compute plane equation of triangle(V0,V1,V2) */
	SUB(E1, V1, V0);    // +0  -3  *0  /0 =3
	SUB(E2, V2, V0);    // +0  -3  *0  /0 =3
	CROSS(N1, E1, E2);  // +0  -3  *9  /0 =3
	d1 = -DOT(N1, V0);  // +2  -1  *3  /0 =1
	/* plane equation 1: N1.X+d1=0 */

	/* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
	du0 = DOT(N1,U0) + d1;  // +3  -0  *3  /0 =1
	du1 = DOT(N1,U1) + d1;  // +3  -0  *3  /0 =1
	du2 = DOT(N1,U2) + d1;  // +3  -0  *3  /0 =1

	/* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
	if (FABS(du0) < EPSILON)
		du0 = 0.0;
	if (FABS(du1) < EPSILON)
		du1 = 0.0;
	if (FABS(du2) < EPSILON)
		du2 = 0.0;
#endif
	du0du1 = du0 * du1;    // +0  -0  *1  /0 =1
	du0du2 = du0 * du2;    // +0  -0  *1  /0 =1

	if (du0du1 > 0.0f && du0du2 > 0.0f) /* same sign on all of them + not equal 0 ? */
		return 0; /* no intersection occurs */

	/* compute plane of triangle (U0,U1,U2) */
	SUB(E1, U1, U0);
	SUB(E2, U2, U0);
	CROSS(N2, E1, E2);
	d2 = -DOT(N2, U0);
	/* plane equation 2: N2.X+d2=0 */

	/* put V0,V1,V2 into plane equation 2 */
	dv0 = DOT(N2,V0) + d2;
	dv1 = DOT(N2,V1) + d2;
	dv2 = DOT(N2,V2) + d2;

#if USE_EPSILON_TEST==TRUE
	if (FABS(dv0) < EPSILON)
		dv0 = 0.0;
	if (FABS(dv1) < EPSILON)
		dv1 = 0.0;
	if (FABS(dv2) < EPSILON)
		dv2 = 0.0;
#endif

	dv0dv1 = dv0 * dv1;
	dv0dv2 = dv0 * dv2;

	if (dv0dv1 > 0.0f && dv0dv2 > 0.0f) /* same sign on all of them + not equal 0 ? */
		return 0; /* no intersection occurs */

	/* compute direction of intersection line */
	CROSS(D, N1, N2);   // +0  -3  *9  /0 =3

	/* compute and index to the largest component of D */
	max = (VALUE) FABS(D[0]);
	index = 0;
	bb = (VALUE) FABS(D[1]);
	cc = (VALUE) FABS(D[2]);
	if (bb > max)
		max = bb, index = 1;
	if (cc > max)
		max = cc, index = 2;

	/* this is the simplified projection onto L*/
	vp0 = V0[index];
	vp1 = V1[index];
	vp2 = V2[index];

	up0 = U0[index];
	up1 = U1[index];
	up2 = U2[index];

	/* compute interval for triangle 1 */
	VALUE a, b, c, x0, x1;
	NEWCOMPUTE_INTERVALS(vp0, vp1, vp2, dv0, dv1, dv2, dv0dv1, dv0dv2, a, b, c,
			x0, x1);  // +0  -4  *2  /0 =5

	/* compute interval for triangle 2 */
	VALUE d, e, f, y0, y1;
	NEWCOMPUTE_INTERVALS(up0, up1, up2, du0, du1, du2, du0du1, du0du2, d, e, f,
			y0, y1);

	VALUE xx, yy, xxyy, tmp;
	xx = x0 * x1;
	yy = y0 * y1;
	xxyy = xx * yy;

	tmp = a * xxyy;
	isect1[0] = tmp + b * x1 * yy;
	isect1[1] = tmp + c * x0 * yy;

	tmp = d * xxyy;
	isect2[0] = tmp + e * xx * y1;
	isect2[1] = tmp + f * xx * y0;

	SORT(isect1[0], isect1[1]);
	SORT(isect2[0], isect2[1]);

	if (isect1[1] < isect2[0] || isect2[1] < isect1[0])
		return 0;
	return 1;
}

// IMPLIMENT OF
// --------------------------------------------------------
// [a, b, c, d]
//   | ax ay az 1 |
// = | bx by bz 1 |
//   | cx cy cz 1 |
//
//   | ax-dx ay-dy az-dz |
// = | bx-dx by-dy bz-dz |
//   | cx-dx cy-dy cz-dz |
template<class VALUE>
VALUE ORIENT3D(VALUE *pa, VALUE *pb, VALUE *pc, VALUE *pd)
{
	VALUE adx, bdx, cdx;
	VALUE ady, bdy, cdy;
	VALUE adz, bdz, cdz;

	adx = pa[0] - pd[0];
	bdx = pb[0] - pd[0];
	cdx = pc[0] - pd[0];
	ady = pa[1] - pd[1];
	bdy = pb[1] - pd[1];
	cdy = pc[1] - pd[1];
	adz = pa[2] - pd[2];
	bdz = pb[2] - pd[2];
	cdz = pc[2] - pd[2];

	return adx * (bdy * cdz - bdz * cdy) + bdx * (cdy * adz - cdz * ady)
			+ cdx * (ady * bdz - adz * bdy);
}
// the method 2 is based on the Devillers O, Guigue P. Faster triangle-triangle intersection tests[J]. 2002.
//
//template specialization ---------------------------------
#include "ts_predicates.h"
//template <>
//double ORIENT3D<double>(double *pa,
//		double *pb,
//		double *pc,
//		double *pd)
//{
//	return orient3d(pa, pb, pc, pd);
//}


inline short COUNT( //
		const short &v, //
		const short& a, const short& b, const short& c)
{
	short res = 0;
	if (v == a) {
		res++;
	}
	if (v == b) {
		res++;
	}
	if (v == c) {
		res++;
	}
	return res;
}
template<class VALUE>
inline short SIGN3(const VALUE& v)
{
	if (v < 0) {
		return -1;
	} else if (v == 0) {
		return 0;
	} else {
		return 1;
	}
}
template<class VALUE>
inline void _reverse(VALUE*& V0, VALUE*& V1, VALUE*& V2)
{
	VALUE* tmp;
	tmp = V1;
	V1 = V2;
	V2 = tmp;
}
template<class VALUE>
inline void _permutation(VALUE* V0, VALUE* V1, VALUE* V2, //Triangle V
		short arr[3], short flag, //array sign
		VALUE*& p, VALUE*& q, VALUE*& r)
{
	if (arr[0] == flag) {
		p = V0;
		q = V1;
		r = V2;
		return;
	}
	if (arr[1] == flag) {
		p = V1;
		q = V2;
		r = V0;
		return;
	}
	if (arr[2] == flag) {
		p = V2;
		q = V0;
		r = V1;
		return;
	}
}
template<class VALUE>
inline void show(String name, VALUE v[3])
{
	std::cout << name << ' ' << v[0] << ' ' << v[1] << ' ' << v[2] << "\n";
}

template<class VALUE>
inline int general_case( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //Triangle V
		VALUE U0[3], VALUE U1[3], VALUE U2[3], //Triangle U
		short arr1[3], short arr2[3])
{
	short sum1 = arr1[0] + arr1[1] + arr1[2];
	VALUE *UT0 = U0, *UT1 = U1, *UT2 = U2;
	VALUE *VT0 = V0, *VT1 = V1, *VT2 = V2;

	//std::cout << "arr1 " << arr1[0] << ' ' << arr1[1] << ' ' << arr1[2]
	//		<< std::endl;
	//std::cout << "arr2 " << arr2[0] << ' ' << arr2[1] << ' ' << arr2[2]
	//		<< std::endl;
	if (sum1 == 1) {
		//reverse
		//std::cout << "reverse 1\n";
		_permutation(U0, U1, U2, arr1, -1, UT0, UT1, UT2);
		//_reverse(V0, V1, V2, arr1, p1, q1, r1);
	} else {
		//premutation
		//std::cout << "permuta 1\n";
		_permutation(U0, U1, U2, arr1, 1, UT0, UT1, UT2);
	}

	short sum2 = arr2[0] + arr2[1] + arr2[2];

	if (sum2 == 1) {
		//reverse
		//std::cout << "reverse 2\n";
		_permutation(V0, V1, V2, arr2, -1, VT0, VT1, VT2);

		//_reverse(U0, U1, U2, arr2, p2, q2, r2);
	} else {
		//premutation
		//std::cout << "permuta 2\n";
		_permutation(V0, V1, V2, arr2, 1, VT0, VT1, VT2);
	}
	if (sum1 == 1) {
		_reverse(VT0, VT1, VT2);
	}
	if (sum2 == 1) {
		_reverse(UT0, UT1, UT2);
	}
//decision tree
//i k

	//std::cout << "U->V " << SIGN3(ORIENT3D(UT0, VT1, VT2, VT0)) << std::endl;
	//std::cout << "     " << SIGN3(ORIENT3D(UT1, VT1, VT2, VT0)) << std::endl;
	//std::cout << "     " << SIGN3(ORIENT3D(UT2, VT1, VT2, VT0)) << std::endl;

	//std::cout << "V->U " << SIGN3(ORIENT3D(VT0, UT1, UT2, UT0)) << std::endl;
	//std::cout << "     " << SIGN3(ORIENT3D(VT1, UT1, UT2, UT0)) << std::endl;
	//std::cout << "     " << SIGN3(ORIENT3D(VT2, UT1, UT2, UT0)) << std::endl;

	VALUE cik = ORIENT3D(VT1, UT1, VT0, UT0);
	//show("UT0 ", VT0);
	//show("UT1 ", VT1);
	//show("UT2 ", VT2);

	//std::cout << "step 1 " << SIGN3(cik) << std::endl;
	if (cik > 0) {
		return -1;
	} else { //cik<0 i>k
		VALUE cjl = ORIENT3D(VT2, UT2, VT0, UT0);
		//std::cout << "step 2 " << SIGN3(cjl) << std::endl;
		if (cjl < 0) {
			return -1;
		} else { //cjl>0 j<l
			return 1;  //intersect
		}
	}
}

template<class VALUE>
inline void _clear_and_rebuild_pointer(VALUE*& X, VALUE*& Y, VALUE*& Z,
		short& len)
{
	if (X != NULL) {
		delete X;
	}
	if (Y != NULL) {
		delete Y;
	}
	if (Z != NULL) {
		delete Z;
	}
	X = new VALUE[len];
	Y = new VALUE[len];
	Z = new VALUE[len];
}

template<class VALUE>
inline void _plane_equation(VALUE* V0, VALUE* V1, VALUE* V2, VALUE* VN,
		VALUE& D)
{
	VALUE v01[3], v02[3];
	SUB(v01, V1, V0);
	SUB(v02, V2, V0);
	CROSS(VN, v01, v02);
	D = DOT(VN, V0);
}

template<class VALUE>
inline void _intersect_line_plane(VALUE* LV, VALUE* LD, VALUE* PN,
		VALUE& PD, VALUE* RES)
{
	const VALUE& A = PN[0];
	const VALUE& B = PN[1];
	const VALUE& C = PN[2];
	const VALUE& D = PD;
	VALUE t;
	VALUE s = (D - A * LD[0] - B * LD[1] - C * LD[2]);
	VALUE d = (A * LV[0] + B * LV[1] + C * LV[2]);
	if (d == 0) {
		t = s / 1e-10;
	} else {
		t = s / d;
	}
	RES[0] = LD[0] + LV[0] * t;
	RES[1] = LD[1] + LV[1] * t;
	RES[2] = LD[2] + LV[2] * t;
}

template<class VALUE>
inline int general_case_calculation( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //Triangle V
		VALUE U0[3], VALUE U1[3], VALUE U2[3], //Triangle U
		short arr1[3], short arr2[3],          //orientation
		VALUE*& X, VALUE*& Y, VALUE*& Z, short& len)
{
	short sum1 = arr1[0] + arr1[1] + arr1[2];
	short sum2 = arr2[0] + arr2[1] + arr2[2];
	VALUE *UT0 = U0, *UT1 = U1, *UT2 = U2;
	VALUE *VT0 = V0, *VT1 = V1, *VT2 = V2;

	if (sum1 == 1) {
		_permutation(U0, U1, U2, arr1, -1, UT0, UT1, UT2);
	} else {
		_permutation(U0, U1, U2, arr1, 1, UT0, UT1, UT2);
	}
	if (sum2 == 1) {
		_permutation(V0, V1, V2, arr2, -1, VT0, VT1, VT2);
	} else {
		_permutation(V0, V1, V2, arr2, 1, VT0, VT1, VT2);
	}
	if (sum1 == 1) {
		_reverse(VT0, VT1, VT2);
	}
	if (sum2 == 1) {
		_reverse(UT0, UT1, UT2);
	}
	//decision tree
	//i k
	VALUE cik = ORIENT3D(VT1, UT1, VT0, UT0);
	if (cik > 0) {
		return -1;
	} else {          //cik<0 i>k
		VALUE cjl = ORIENT3D(VT2, UT2, VT0, UT0);
		if (cjl < 0) {
			return -1;
		} else { //cjl>0 j<l

			VALUE cjk = ORIENT3D(VT1, UT2, VT0, UT0);
			VALUE cil = ORIENT3D(VT2, UT1, VT0, UT0);
			VALUE LV[3];
			VALUE* LD;
			VALUE PV[3], PD;
			if (cjk >= 0) { // get k
				SUB(LV, VT1, VT0);
				LD = VT0;
				_plane_equation(UT0, UT1, UT2, PV, PD);
			} else {       //get j
				SUB(LV, UT2, UT0);
				LD = UT0;
				_plane_equation(VT0, VT1, VT2, PV, PD);
			}
			VALUE PS[3];
			_intersect_line_plane(LV, LD, PV, PD, PS);
			if (cil >= 0) { // get i
				SUB(LV, UT1, UT0);
				LD = UT0;
				_plane_equation(VT0, VT1, VT2, PV, PD);
			} else {     //get l
				SUB(LV, VT2, VT0);
				LD = VT0;
				_plane_equation(UT0, UT1, UT2, PV, PD);
			}
			VALUE PE[3];
			_intersect_line_plane(LV, LD, PV, PD, PE);
			len = 2;
			_clear_and_rebuild_pointer(X, Y, Z, len);
			X[0] = PS[0];
			X[1] = PE[0];
			Y[0] = PS[1];
			Y[1] = PE[1];
			Z[0] = PS[2];
			Z[1] = PE[2];

			return 1;  //intersect
		}
	}
}

//step 2 analyse cases
// return code
// 1            same sign
// 2 2D case    co-plane
// 3 2D case    two points on plane
// 4 2D case    one point  on plane other tow at same side
// 5 3D case    one point  on plane other tow at different side
// 6 3D case    general case
//                 *  +
// -1 -1 -1  3     - -3 -> same sign
// -1 -1  0  2 1     -2 =
// -1 -1  1  2   1 + -1 ----->
// -1  0 -1  2 1     -2 =
// -1  0  0  1 2     -1 ==
// -1  0  1  1 1 1      ----->
// -1  1 -1  2   1 + -1 ----->
// -1  1  0  1 1 1      ----->
// -1  1  1  1   2 -  1 ----->
//  0 -1 -1  2 1     -2 =
//  0 -1  0  1 2     -1 ==
//  0 -1  1  1 1 1      ----->
//  0  0 -1  1 2     -1 ==
//  0  0  0    3        -------------------------
//  0  0  1    2 1    1
//  0  1 -1  1 1 1
//  0  1  0    2 1    1
//  0  1  1    1 2    2
//  1 -1 -1  2   1 +  1
//  1 -1  0  1 1 1
//  1 -1  1  1   2 -  1
//  1  0 -1  1 1 1
//  1  0  0    2 1    1
//  1  0  1    1 2    2
//  1  1 -1  1   2 -  1
//  1  1  0    1 2    2
//  1  1  1      3 +  3 -> same sign

inline short analyse_cases(short arr[3])
{
	short sum = arr[0] + arr[1] + arr[2];
//short mtp = arr[0] * arr[1] * arr[2];
	short abssum = ABS(sum);
//short absmtp = ABS(mtp);

//case 1 same sign
	if (abssum == 3) {
		return 1;
	}
//case 2 all equal to zero, co-plane
	short num_zero = COUNT(0, arr[0], arr[1], arr[2]);
	if (num_zero == 3) {
		return 2;  //==>>> 2D
	}
//case 3 diffrent signs
//3-2 2 zero
	if (num_zero == 2) {  //two points on PI
		return 3;  //==>>> 2D
	}
//3-1 1 zero
	if (num_zero == 1) { //only one point on PI
		if (abssum == 2) {
			return 4; //other two are at same side  //==>>> 2D
		}
		return 5;
		//other two are at different side
	}
	return 6;  //general case 3D
}

template<class VALUE>
int TriTriIsect_Guigue( //
		VALUE V0[3], VALUE V1[3], VALUE V2[3], //Triangle V
		VALUE U0[3], VALUE U1[3], VALUE U2[3]) //Triangle U
{
//T1 intersect PI2
	VALUE oVU[3]; //
	oVU[0] = ORIENT3D(U0, V1, V2, V0);
	oVU[1] = ORIENT3D(U1, V1, V2, V0);
	oVU[2] = ORIENT3D(U2, V1, V2, V0);
	short sVU[3];
	sVU[0] = SIGN3(oVU[0]);
	sVU[1] = SIGN3(oVU[1]);
	sVU[2] = SIGN3(oVU[2]);

	short caseVU = analyse_cases(sVU);

//T1 intersect PI2
	if (caseVU == 1) {
		return -1;  // no intersect---------------------
	} else if (2 <= caseVU && caseVU <= 4) { // 2D condition;
	//==>>> 2D
		std::cout << "2D case\n";
		return -2;
	} else {  //Additional tests
		VALUE oUV[3]; //
		oUV[0] = ORIENT3D(V0, U1, U2, U0);
		oUV[1] = ORIENT3D(V1, U1, U2, U0);
		oUV[2] = ORIENT3D(V2, U1, U2, U0);
		short sUV[3];
		sUV[0] = SIGN3(oUV[0]);
		sUV[1] = SIGN3(oUV[1]);
		sUV[2] = SIGN3(oUV[2]);
		short caseUV = analyse_cases(sUV);
		if (caseUV == 1) {
			return -1;  // no intersect ------not reach--
		} else if (2 <= caseUV && caseUV <= 4) {  // 2D condition;
		//==>>> 2D
			std::cout << "2D case\n";
			return -2;
		} else {  // 3D condition
			return general_case(V0, V1, V2, U0, U1, U2, sVU, sUV);
		}
	}
	return -1; // no intersect----------------------------
}

template<class VALUE>
int TriTriIsect_Guigue_calculation(
		//
		VALUE V0[3], VALUE V1[3],
		VALUE V2[3], //Triangle V
		VALUE U0[3], VALUE U1[3], VALUE U2[3], VALUE*& X, VALUE*& Y, VALUE*& Z,
		short& len) //Triangle U
{
//T1 intersect PI2
	VALUE oVU[3]; //
	oVU[0] = ORIENT3D(U0, V1, V2, V0);
	oVU[1] = ORIENT3D(U1, V1, V2, V0);
	oVU[2] = ORIENT3D(U2, V1, V2, V0);
	short sVU[3];
	sVU[0] = SIGN3(oVU[0]);
	sVU[1] = SIGN3(oVU[1]);
	sVU[2] = SIGN3(oVU[2]);

	short caseVU = analyse_cases(sVU);

//T1 intersect PI2
	if (caseVU == 1) {
		return -1;  // no intersect---------------------
	} else if (2 <= caseVU && caseVU <= 4) { // 2D condition;
	//==>>> 2D
		std::cout << "2D case\n";
		return -2;
	} else {  //Additional tests
		VALUE oUV[3]; //
		oUV[0] = ORIENT3D(V0, U1, U2, U0);
		oUV[1] = ORIENT3D(V1, U1, U2, U0);
		oUV[2] = ORIENT3D(V2, U1, U2, U0);
		short sUV[3];
		sUV[0] = SIGN3(oUV[0]);
		sUV[1] = SIGN3(oUV[1]);
		sUV[2] = SIGN3(oUV[2]);
		short caseUV = analyse_cases(sUV);
		if (caseUV == 1) {
			return -1;  // no intersect ------not reach--
		} else if (2 <= caseUV && caseUV <= 4) {  // 2D condition;
		//==>>> 2D
			std::cout << "2D case\n";
			return -2;
		} else {  // 3D condition
			return general_case_calculation(V0, V1, V2, U0, U1, U2, sVU, sUV, X,
					Y, Z, len);
		}
	}
	return -1; // no intersect----------------------------
}

}

#endif /* TS_TRI_MOLLER_H_ */
