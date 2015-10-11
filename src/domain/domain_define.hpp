#ifndef GRIDDEF_H_
#define GRIDDEF_H_

#include "../carpio_define.hpp"
#include "../algebra/arithmetic.hpp"

namespace carpio {

    typedef short code;


    /*enum Direction {
        //D_X==================
                D_MMX_M = 8 + 0, //001 000
        D_MPX_M = 8 + 2, //001 010
        D_PMX_M = 8 + 4, //001 100
        D_PPX_M = 8 + 6, //001 110
        //
                D_MMX_P = 8 + 1, //001 001
        D_MPX_P = 8 + 3, //001 011
        D_PMX_P = 8 + 5, //001 101
        D_PPX_P = 8 + 7, //001 111

        //D_Y==================
                D_MYM_M = 16 + 0, //010 000
        D_MYP_M = 16 + 1, //010 001
        D_PYM_M = 16 + 4, //010 100
        D_PYP_M = 16 + 5, //010 101
        //
                D_MYM_P = 16 + 2, //010 010
        D_MYP_P = 16 + 3, //010 011
        D_PYM_P = 16 + 6, //010 110
        D_PYP_P = 16 + 7, //010 111

        //D_Z==================
        D_ZMM_M = 32 + 0, //100 000
        D_ZMP_M = 32 + 1, //100 001
        D_ZPM_M = 32 + 2, //100 010
        D_ZPP_M = 32 + 3, //100 011
        //
        D_ZMM_P = 32 + 4, //100 100
        D_ZMP_P = 32 + 5, //100 101
        D_ZPM_P = 32 + 6, //100 110
        D_ZPP_P = 32 + 7, //100 111

        //D_ZX==================
        D_ZMX_MM = 40 + 0, //101 000
        D_ZMX_MP = 40 + 1, //101 001
        D_ZMX_PM = 40 + 4, //101 100
        D_ZMX_PP = 40 + 5, //101 101
        //
        D_ZPX_MM = 40 + 2, //101 010
        D_ZPX_MP = 40 + 3, //101 011
        D_ZPX_PM = 40 + 6, //101 110
        D_ZPX_PP = 40 + 7, //101 111

        //D_ZY==================
        D_ZYM_MM = 48 + 0, //110 000
        D_ZYM_MP = 48 + 2, //110 010
        D_ZYM_PM = 48 + 4, //110 100
        D_ZYM_PP = 48 + 6, //110 110
        //
        D_ZYP_MM = 48 + 1, //110 001
        D_ZYP_MP = 48 + 3, //110 011
        D_ZYP_PM = 48 + 5, //110 101
        D_ZYP_PP = 48 + 7, //110 111

        //D_YX==================
        D_MYX_MM = 24 + 0, //011 000
        D_MYX_MP = 24 + 1, //011 001
        D_MYX_PM = 24 + 2, //011 010
        D_MYX_PP = 24 + 3, //011 011
        //
                D_PYX_MM = 24 + 4, //011 100
        D_PYX_MP = 24 + 5, //011 101
        D_PYX_PM = 24 + 6, //011 110
        D_PYX_PP = 24 + 7, //011 111

        //D_ZYX==================
        D_ZYX_MMM = 56 + 0, //111 000
        D_ZYX_MMP = 56 + 1, //111 001
        D_ZYX_MPM = 56 + 2, //111 010
        D_ZYX_MPP = 56 + 3, //111 011
        D_ZYX_PMM = 56 + 4, //111 100
        D_ZYX_PMP = 56 + 5, //111 101
        D_ZYX_PPM = 56 + 6, //111 110
        D_ZYX_PPP = 56 + 7, //111 111
    };*/

    enum Orientation {
        _M_ = 0, //
        _P_ = 1, //
        _C_ = 2, //
    };

    enum Axes {
        _X_ = 0, //
        _Y_ = 1, //
        _Z_ = 2, //
    };

    enum Plane {
        _XY_ = 24,
        _YZ_ = 48,
        _ZX_ = 40,
    };


    typedef unsigned short Direction;

    inline Direction ToDirection(const Plane &p,
                                 const Orientation &o1,
                                 const Orientation &o2) {
        ASSERT(o1 != _C_);
        ASSERT(o2 != _C_);
        switch (p) {
            case _XY_:
                return p + o1 + (o2 << 1);
            case _YZ_:
                return p + (o1 << 1) + (o2 << 2);
            case _ZX_:
                return p + (o1 << 2) + o2;
        };
    }

    inline Direction ToDirection(const Orientation &x,
                                 const Orientation &y,
                                 const Orientation &z) {
        ASSERT(x != _C_);
        ASSERT(y != _C_);
        ASSERT(z != _C_);
        return 56 + x + (y << 1) + (z << 2);
    }

    inline unsigned short HI(const Direction &d) {
        return d >> 3;
    }

    inline unsigned short LO(const Direction &d) {
        return d & 7;
    }

    static const short COUNT_1[8] = {
            0, 1, 1, 2, 1, 2, 2, 3
    };

    inline bool IsFaceDirection(const Direction &d) {
        return COUNT_1[(d >> 3)] == 1;
    }

    inline bool IsFacePDirection(const Direction &d) {
        unsigned short hi = d >> 3;
        unsigned short low = d & 7;
        return (COUNT_1[hi] == 1) && (hi & low) != 0;
    }

    inline bool IsXYDirection(const Direction &d) {
        return (d >> 3) == 3;
    }

    inline bool IsYZDirection(const Direction &d) {
        return (d >> 3) == 6;
    }

    inline bool IsZXDirection(const Direction &d) {
        return (d >> 3) == 5;
    }

    inline bool IsPlaneDirection(const Direction &d) {
        return COUNT_1[(d >> 3)] == 2;
    }

    inline bool IsXYZDirection(const Direction &d) {
        return (d >> 3) == 7;
    }

    inline Direction FaceDirectionInOrder(const size_t i) {
        ASSERT(i < 6);
        static const Direction ARR_FD[] = {8, 9, 16, 18, 32, 36};
        return ARR_FD[i];
    }

    inline Direction XYDirectionInOrder(const size_t i) {
        ASSERT(i < 4);
        return 24 + i;
    }

    inline Direction YZDirectionInOrder(const size_t i) {
        ASSERT(i < 4);
        return 48 + i * 2;
    }

    inline Direction ZXDirectionInOrder(const size_t i) {
        ASSERT(i < 4);
        static const Direction ARR_ZXD[] = {
                40 + 0, //101 000
                40 + 1, //101 001
                40 + 4, //101 100
                40 + 5, //101 101
        };
        return ARR_ZXD[i];
    }

    inline Direction XYZDirectionInOrder(const size_t i) {
        ASSERT(i < 8);
        return 56 + i;
    }

    inline Direction DirectionInOrder(const size_t i){
        ASSERT(i < 26);
        if (i<6) return FaceDirectionInOrder(i);
        if (i<10) return XYDirectionInOrder(i-6);
        if (i<14) return YZDirectionInOrder(i-10);
        if (i<18) return ZXDirectionInOrder(i-14);
        return XYZDirectionInOrder(i-18);
    }


//default type
    typedef double CooValueType;
    typedef double ValueType;


}

#endif /* GRIDDEF_H_ */
