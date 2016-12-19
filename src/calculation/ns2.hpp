#ifndef _NS_H_
#define _NS_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "../utility/Clock.h"
#include "../io/mmio.h"

#include <vector>
#include <memory>

//#define __Debug__

#define __x__  -0.4
#define __y__  0.45

namespace carpio {

template<typename COO_VALUE, typename VALUE, int DIM>
class NS_: public Equation_<COO_VALUE, VALUE, DIM> {
public:
	static const st Dim = DIM;
	static const st NumFaces = DIM + DIM;
	static const st NumVertexes = (DIM == 3) ? 8 : (DIM + DIM);
	static const st NumNeighbors = NumFaces;

	typedef Equation_<COO_VALUE, VALUE, DIM> Equation;
	typedef Equation_<COO_VALUE, VALUE, DIM>& ref_Equation;

	typedef Poisson_<COO_VALUE, VALUE, DIM> Self;
	typedef Poisson_<COO_VALUE, VALUE, DIM>& ref_Self;
	typedef Poisson_<COO_VALUE, VALUE, DIM>* pSelf;
	typedef const Poisson_<COO_VALUE, VALUE, DIM>* const_pSelf;

	typedef typename Equation::vt vt;
	typedef typename Equation::cvt cvt;

	typedef typename Equation::pDomain pDomain;
	typedef typename Equation::spEvent spEvent;

	typedef typename Equation::Grid Grid;
	typedef typename Equation::pGrid pGrid;
	typedef typename Equation::const_pGrid const_pGrid;

	typedef typename Equation::Node Node;
	typedef typename Equation::pNode pNode;
	typedef typename Equation::const_pNode const_pNode;

	typedef typename Equation::Exp Exp;
	typedef typename Equation::pExp pExp;
	typedef typename Equation::spExp spExp;
	typedef typename Equation::spFace spFace;

	typedef typename Equation::Face Face;
	typedef typename Equation::pFace pFace;
	typedef typename Equation::const_pFace const_pFace;

	typedef typename Equation::Flag Flag;

	typedef typename Equation::Mat Mat;
	typedef typename Equation::Arr Arr;

	typedef typename Equation::BoundaryCondition BoundaryCondition;
	typedef typename Equation::pBoundaryCondition pBoundaryCondition;
	typedef typename Equation::const_pBoundaryCondition const_pBoundaryCondition;

	typedef typename Equation::Function Function;
	typedef typename Equation::Fun_get_spExp Fun_get_spExp;

	typedef ArrayListT<spExp> Arr_spExp;
	typedef typename Equation::Map_FE Map_FE;
	typedef typename Equation::pMap_FE pMap_FE;
	typedef typename Equation::Iter_FE Iter_FE;
	typedef typename Equation::const_Iter_FE const_Iter_FE;
	typedef typename Equation::Pair_FE Pair_FE;
	typedef typename Equation::Ret_FE Ret_FE;

protected:

	struct FaceData {
		spExp exp_val;
		spExp exp_gra;
		vt veo;
		vt veot;
	};

	typedef std::map<spFace, FaceData, comp_spFace> Map_FD;
	typedef std::map<spFace, FaceData, comp_spFace>* pMap_FD;
	typedef typename std::map<spFace, FaceData, comp_spFace>::iterator Iter_FD;
	typedef typename std::map<spFace, FaceData, comp_spFace>::const_iterator const_Iter_FD;
	typedef std::pair<spFace, FaceData> Pair_FD;
	typedef std::pair<Iter_FE, bool> Ret_FD;

	vt _tol;
	int _max_iter;

	//matrx
	typedef ArrayListV<st> Arr_st;

	st _scheme;
	std::vector<Limiter> _v_limiter;

#ifdef __Debug__
	std::fstream debug_log;
#endif
public:
	/*
	 * constructor
	 *   2D
	 *
	 * pressure       0
	 * velocity  u    1
	 * velocity  v    2
	 * velocity  u*   3
	 * velocity  v*   4
	 * source    sx   5
	 * source    sy   6
	 * Density   rho  7
	 * viscocity mu   8
	 *
	 * utp
	 * utp map          0
	 */
	NS_(pDomain pf) :
			_pdomain(pf) {
		_push_limiter();
		this->_scheme = 0;
	}

	~NS_() {
	}

	int initial() {
		if (this->is_stand_alone()) {
			this->_construct_stand_alone();
			// set time step
			if (this->has_time_term()) {

				this->_pdomain->resize_data(this->max_idx(this->_vars_c) + ss,
						0, 0, this->max_idx(this->_vars_ut) + 1);
			} else {
				std::cerr << ">! no time term";
				SHOULD_NOT_REACH;
			}
			this->init_mu();
			this->init_rho();

		} else {
			SHOULD_NOT_REACH;
		}

		this->_build_FD_on_leaf();
		this->_build_NE_on_leaf(1);
		this->_build_fun_get();

		return -1;
	}

	void init_mu() {
		if (this->has("uniform_mu")) {
			return;
		}

		if (this->has_function("set_mu")) {
			Function f = this->_functions["set_mu"];
			this->set_val(f, this->_vars_c["mu"]);
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void init_rho() {
		if (this->has("uniform_rho")) {
			return;
		}

		if (this->has_function("set_rho")) {
			Function f = this->_functions["set_rho"];
			this->set_val(f, this->_vars_c["rho"]);
		} else {
			SHOULD_NOT_REACH;
		}
	}

	/*
	 * build face exp to utp
	 */
	int solve(vt tol = 1e-6, int max_iter = 1000, int info = 0) {
		std::list<vt> lr;
		int rcode = this->solve(tol, max_iter, lr, info);
		//std::cout<< "max iter "<<max_iter<< " tol "<< tol <<"\n";
		return rcode;
	}

	int solve(vt& tol, int& max_iter, std::list<vt>& lr, int info = 0) {
		Clock clock;
#ifdef __Debug__
		debug_log << "->solve at line " << __LINE__ << "\n";
#endif

#ifdef __Debug__
		debug_log << "  solve: 1 build u star on center\n";
#endif
		for (st i = 0; i < Dim; ++i) {
#ifdef __Debug__
			debug_log << "  solve: on axe" << ToString(ToAxes(i)) << "\n";
#endif
			_u_star_on_center(ToAxes(i), _dt);
		}
		clock.break_point("u* on center");
#ifdef __Debug__
		debug_log << "  solve: 2 build u star on boundary\n";
#endif
		// set u star bc
		for (st i = 0; i < Dim; ++i) {
#ifdef __Debug__
			debug_log << "  solve: on axe "<< ToString(ToAxes(i)) <<"\n";
#endif
			_u_star_on_ghost(ToAxes(i), _dt);
		}
		clock.break_point("u* on ghost");

		Mat mat;
		Arr b;
		_build_pressure_equation(mat, b);
		clock.break_point("build pressure equation");
		Arr x(b.size());
		// initial x
		_grid_to_arr(x, _p_idx);
		clock.break_point("copy array to x");
#ifdef __Debug__

//		b.show();
//		std::list<Gnuplot_actor> lga;
//		Gnuplot_actor ga;
//		GnuplotActor_MatrixSCR(ga, mat);
//		lga.push_back(ga);
//		Gnuplot gp;
//		gp.set_equal_ratio();
//		gp.plot(lga);

#endif
		//std::cout<<" diagnally dominant " << mat.is_diagonally_dominant() <<"\n";
		//mm_write_mtx_sparse("mat", mat);
		int sf = Dia_BiCGSTAB(mat, x, b, max_iter, tol, lr);

		clock.break_point("solve");
		//int sf = Jacobi(mat, x, b, max_iter, tol, lr);
		if (info == 1) {
			std::cout << "iter n " << lr.size() << " residual " << tol << "\n";
		}
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
		}
		//x.show();
		_arr_to_grid(x, _p_idx);
		clock.break_point("copy x to grid");
		// set Pressure star bc
		for (st i = 0; i < Dim; ++i) {
			this->_pdomain->set_val_ghost_by_bc(_p_idx);
		}
		clock.break_point("set pressure on ghost");
		_correct_veo_on_center();
		clock.break_point("correct veo");
		if (info == 1) {
			clock.show();
		}
#ifdef __Debug__
		debug_log << "<-solve " << "\n";
#endif
		return 1;
	}

	int advance(st max_step, vt tol = 1e-6, int max_iter = 1000, int info = 0) {
#ifdef __Debug__
		debug_log << "->advance at line " << __LINE__ << "\n";
		debug_log << "  tol " << tol << ", max_iter " << max_iter << ", info "
		<< info << "\n";
#endif
		_build_FE_on_leaf();
		tick_t tick_cpu = Clock::SystemTime();
		tick_t tick_wall = Clock::Tick();
		for (st i = 0; i < max_step; ++i) {
			if (info == 1) {
				fmt::print(
						"step: {:>8} dt: {:>8.5f} cpu: {:>8.5f} wall: {:>8.5f}\n",
						i, i * _dt,
						Clock::TimespanToSecondsD(tick_cpu,
								Clock::SystemTime()),
						Clock::TimespanToSecondsD(tick_wall, Clock::Tick()));
				tick_cpu = Clock::SystemTime();
				//tick_wall = Clock::Tick();
			}
#ifdef __Debug__
			debug_log << "  advance step = " << i << "\n";
#endif
			solve(tol, max_iter, info);
		}
		_delete_FE_on_leaf();
		return 1;
	}

	/*
	 * Boundary condition
	 */
	void set_boundary_condition(st si, st segi, st vi, pBoundaryCondition pbc) {
		_pdomain->_pbindex->insert(si, segi, vi, pbc);
	}
	void set_BC_velocity(st si, st segi, Axes a, pBoundaryCondition pbc) {
		this->set_boundary_condition(si, segi, _v_idx[a], pbc);
		// the boundary condition of velo star is the same as the velo
		//this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
	}
	void set_BC_velocity_t(st si, st segi, Axes a, pBoundaryCondition pbc) {
		this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
		// the boundary condition of velo star is the same as the velo
		//this->set_boundary_condition(si, segi, _vt_idx[a], pbc);
	}

	/*
	 * set
	 */
	void set_v_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_v_ghost(Axes a) {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_v_idx[a]);
	}
	void set_velo(Function pfun, Axes a) {
		this->set_v_grid(pfun, a);
		this->set_v_ghost(a);
	}
	void set_source_grid(Function pfun, Axes a) {
		ASSERT(a < Dim);
		_pdomain->set_val_grid(_v_idx[a], pfun);
	}

	void set_source_ghost(Function pfun, Axes a) {
		_pdomain->set_val_ghost(_v_idx[a], pfun);
	}
	void set_source(Function pfun, Axes a) {
		this->set_source_grid(pfun, a);
		this->set_source_ghost(pfun, a);
	}

	void set_rho_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_rho_idx);
	}
	void set_rho_grid(Function pfun) {
		_pdomain->set_val(_rho_idx, pfun);
	}
	void set_rho(Function pfun) {
		this->set_rho_grid(pfun);
		this->set_rho_ghost();
	}
	void set_mu_ghost() {
		// the boundary must be set
		_pdomain->set_val_ghost_by_bc(_mu_idx);
	}
	// can be used to set initial boundary condition
	void set_mu_grid(Function pfun) {
		_pdomain->set_val(_mu_idx, pfun);
	}
	void set_mu(Function pfun) {
		this->set_mu_grid(pfun);
		this->set_mu_ghost();
	}

protected:
	/**
	 *  @brief map fd on utp array
	 *
	 */
	void _new_map_fd(pNode pn, st idx) {
		utPointer& utp = pn->utp(idx);
		if (utp == nullptr) {
			utp = new Map_FE();
		} else {
			Map_FD& mfd = CAST_REF(pMap_FD, utp);
			mfd.clear();
		}
	}

	void _delete_map_fd(pNode pn, st idx) {
		utPointer& utp = pn->utp(idx);
		if (utp != nullptr) {
			pMap_FE mfd = CAST(pMap_FD, utp);
			mfd->clear();
			delete mfd;
			utp = nullptr;
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void _push_limiter() {
		Limiter fou = [](vt r) {return 0;};

		Limiter minmod = [](vt r) {
			return Max(0.0, Min(1.0, r));
		};

		Limiter superbee = [](vt r) {
			return Max(0.0, Min(2.0 * r, 1.0), Min(r, 2.0));
		};

		Limiter vanLeer = [](vt r) {
			return (r + Abs(r)) / (1 + r);
		};

		_v_limiter.push_back(fou);
		_v_limiter.push_back(minmod);
		_v_limiter.push_back(superbee);
		_v_limiter.push_back(vanLeer);
	}

	vt _CFL_number(pNode pn, vt dt) {
		vt veo[Dim];
		vt cfl[Dim];
		for (st i = 0; i < Dim; i++) {
			veo[i] = pn->cdva(_v_idx[i]);
			cfl[i] = veo[i] * dt / pn->d(ToAxes(i));
		}
		return Max(cfl, Dim);
	}

	pNode _find_C(pFace pface, Float veo_f) {
		pNode pori = pface->pori();
		pNode pnei = pface->pnei();
		Direction dir = pface->dir();
		if (IsFacePDirection(dir)) {
			if (veo_f > 0) {
				return pori;  //o
			} else {
				return pnei;
			}
		} else {
			if (veo_f > 0) {
				return pnei;
			} else {
				return pori;
			}
		}
		SHOULD_NOT_REACH;
		return nullptr;
	}

}
;

}

#endif
