#ifndef CELL_H_
#define CELL_H_

#include "domain_define.hpp"

namespace carpio {

    template<typename VALUE, st DIM>
    class Cell_ {
    public:
        static const st Dim = DIM;
        static const st NumFaces = DIM + DIM;
        static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);

        typedef VALUE vt;

        typedef Cell_<VALUE, DIM> Self;

        typedef void (*pfunction)(Self *, utPointer);

    protected:
        vt _center[Dim];
        vt _hd[Dim];
    public:
        /*
         *  constructor
         */
        Cell_() {
            for (st i = 0; i < Dim; ++i) {
                _center[i] = 0.0;
                _hd[i] = 0.0;
            }
        }

        Cell_(const vt &x, const vt &dhx, //
             const vt &y = 0.0, const vt &dhy = 0.0, //
             const vt &z = 0.0, const vt &dhz = 0.0) {
            for (st i = 0; i < Dim; ++i) {
                if (i == 0) {
                    _center[i] = x;
                    ASSERT(dhx > 0.0);
                    _hd[i] = dhx;
                } else if (i == 1) {
                    _center[i] = y;
                    ASSERT(dhy > 0.0);
                    _hd[i] = dhy;
                } else if (i == 2) {
                    _center[i] = z;
                    ASSERT(dhz > 0.0);
                    _hd[i] = dhz;
                }
            }
        }

        /*
         *  get
         */
        inline vt get(const Orientation &ori, const Axes &axes) const {
            vt res = 0.0;
            if (axes >= Dim) {
                return 0.0;
            }
            switch (ori) {
                case _M_: {
                    res = _center[axes] - _hd[axes];
                    break;
                }
                case _C_: {
                    res = _center[axes];
                    break;
                }
                case _P_: {
                    res = _center[axes] + _hd[axes];
                    break;
                }
                default: {
                    break;
                }
            }
            return res;
        }

        inline vt get_d(const Axes &axes) const {
            return 2.0 * _hd[axes];
        }

        inline vt get_hd(const Axes &axes) const {
            return _hd[axes];
        }

        inline vt volume() const {
            vt res = 1.0;
            for (st i = 0; i < Dim; ++i) {
                res *= 2.0 * _hd[i];
            }
            return res;
        }

        bool is_in_on(const vt x,
                      const vt y = 0,
                      const vt z = 0) const {
            return (IsInRange(this->get(_M_, _X_), x, this->get(_P_, _X_), _cc_)
                    && ((Dim >= 2) ?
                        IsInRange(this->get(_M_, _Y_), y, this->get(_P_, _Y_), _cc_) : true)
                    && ((Dim == 3) ?
                        IsInRange(this->get(_M_, _Z_), z, this->get(_P_, _Z_), _cc_) : true)
            );
        }

        void show(pfunction fun = nullptr, utPointer utp = nullptr) const {
            fun(this, utp);
        }
    };


}

#endif /* CELL_H_ */
