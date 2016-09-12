#ifndef _POISSON_H_
#define _POISSON_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "equation.hpp"
#include "event.hpp"

namespace carpio {
//This file use to solve poisson equation
//
//   ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//

/*
 * the Poisson class
 */
template<typename COO_VALUE, typename VALUE, int DIM>
class Poisson_: public Equation_<COO_VALUE, VALUE, DIM> {
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

protected:

	vt _tol;
	st _max_iter;

public:
	/**
	 * @brief constructor
	 *
	 * @param  [pDomain] a point of domain you want to use.
	 *
	 * @param  [dt]      the difference of the time, default value is -1 which means
	 *                   the equation doesn't have time term.
	 * @param  [maxstep] the max time step. default value is 0.
	 *                   the equation doesn't have time term.
	 *
	 **/
	Poisson_(
			pDomain pd,      //
			vt dt = -1,      //
			st maxstep = 0,  //
			int alphai = -1, int betai = -1, int phii = -1, int fi = -1,
			int fei = -1) :
			Equation(pd, dt, maxstep) {
		if (alphai == -1) {
			//stand alone
			this->_construct_stand_alone();
			this->_pdomain->resize_data(this->max_idx(this->_vars_c) + 1, 0, 0,
					this->max_idx(this->_vars_ut) + 1);
		} else {
			this->_construct_depend(alphai, betai, phii, fi, fei);
		}

		this->set_solve_para();
	}

	/**
	 * @brief override function solve
	 *        only one step, no time term
	 */
	int solve() {
		this->_build_FE_on_leaf();


		this->_delete_FE_on_leaf();
		return 1;
	}



	void set_solve_para(vt tol = 1e-6, int max_iter = 1000) {
		this->_tol = tol;
		this->_max_iter = max_iter;
	}

	int run_one_step(st step) {

		return -1;
	}

	/**
	 *  @brief set time term
	 */
	void set_time_term(){
		spEvent pse(new EventFlag_<cvt,vt, Dim>("time", 1));
		this->_events["time"] = pse;
	}
	void unset_time_term(){
		this->_events.erase("time");
	}

protected:
	/*
	 * stand alone constructor
	 * constructor a Poisson class just used for solver only one poisson equation
	 */
	void _construct_stand_alone() {
		// variables on the center of node
		this->_vars_c["phi"] = 0;
		this->_vars_c["beta"] = 1;
		this->_vars_c["alpha"] = 2;
		this->_vars_c["f"] = 3;
		// untype variables on the node
		this->_vars_ut["FE"] = 0;
	}
	/*
	 * depend constructor
	 * constructor a Poisson class as a part of a group of equations
	 */
	void _construct_depend(st ai, st bi, st pi, st fi, st fei) {
		// variables on the center of node
		this->_vars_c["phi"] = pi;
		this->_vars_c["beta"] = bi;
		this->_vars_c["alpha"] = ai;
		this->_vars_c["f"] = fi;
		// untype variables on the node
		this->_vars_ut["FE"] = fei;
	}

	/*
	 * set
	 */
}
;
}

#endif
