#ifndef _POISSON_H_
#define _POISSON_H_

#include "calculation_define.hpp"
#include "expression.hpp"
#include "equation.hpp"
#include "event.hpp"
#include "../algebra/solver_matrix.hpp"

namespace carpio {
//This file use to solve poisson equation
//
//   ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
//
//   tau (d phi(x,y) / dt) + alpha(x,y)•phi(x,y) + ▽•(beta(x,y)  ▽phi(x,y)  ) = f(x,y)     2D --version
//   tau (d phi(x,y) / dt) + alpha(x,y)•phi(x,y) + ▽•(beta(x,y,z)▽phi(x,y,z)) = f(x,y,z)   3D --version
/*
 * the Poisson class
 */

//#define DEGUG


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

	typedef ArrayListT<spExp> Arr_spExp;
	typedef typename Equation::Map_FE Map_FE;
	typedef typename Equation::pMap_FE pMap_FE;
	typedef typename Equation::Iter_FE Iter_FE;
	typedef typename Equation::const_Iter_FE const_Iter_FE;
	typedef typename Equation::Pair_FE Pair_FE;
	typedef typename Equation::Ret_FE Ret_FE;

	typedef typename Equation::Mat Mat;
	typedef typename Equation::Arr Arr;

	typedef typename Equation::BoundaryCondition BoundaryCondition;
	typedef typename Equation::pBoundaryCondition pBoundaryCondition;
	typedef typename Equation::const_pBoundaryCondition const_pBoundaryCondition;

	typedef typename Equation::Function Function;
	typedef typename Equation::Fun_pNode Fun_pNode;
	typedef typename Equation::Fun_get_spExp Fun_get_spExp;

protected:

	vt _tol;
	int _max_iter;

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
			Equation(pd) {
		if (alphai == -1) {
			//stand alone
			//this->_construct_stand_alone();
			//this->_pdomain->resize_data(this->max_idx(this->_vars_c) + 1, 0, 0,
			//		this->max_idx(this->_vars_ut) + 1);
			this->set_uniform_beta();
		} else {
			this->_construct_depend(alphai, betai, phii, fi, fei);
		}
		this->set_solve_para();
	}

	Exp _face_exp(Face& f, const Exp& exp_val, const Exp& exp_gra) {
		// beta * ▽(phi(x,y))
		vt beta_f;
		if (!this->has("uniform_beta")) {
			beta_f = exp_val.substitute(this->_vars_c["beta"]);
		} else {
			// beta is uniform, no interpolation
			beta_f = this->_values["uniform_beta"];
			//std::cout << "beta " << beta_f << std::endl;
			//beta_f = f.pori()->cda(this->_vars_c["beta"]);
		}
		Exp res(exp_gra);
		int sign = IsFacePDirection(f.dir()) ? 1 : -1;
		res.times(sign * beta_f);
		return std::move(res);
	}

	Exp _f_term(pNode pn) {
		Exp exp;
		exp.insert(-(pn->cda(this->_vars_c["f"]) * pn->volume()), pn, 0);
		return std::move(exp);
	}

	Exp _alpha_term(pNode pn) {
		Exp exp;
		Vt a_f = pn->cda(this->_vars_c["alpha"]);
		exp.insert(a_f * pn->volume(), pn, 1);
		return std::move(exp);
	}

	void _node_exp(pNode pn, Exp& nexp) {
		nexp.clear();
		utPointer& utp = pn->utp(this->_vars_ut["FE"]);
		ASSERT(utp != nullptr);
		Map_FE& mfe = CAST_REF(pMap_FE, utp);
		for (Iter_FE iter = mfe.begin(); iter != mfe.end(); ++iter) {
			const spFace& spf = iter->first;
			Arr_spExp & arr_spexp = iter->second;
			Exp fexp = this->_face_exp(*spf, *(arr_spexp[0]), *(arr_spexp[1]));
			fexp.times(spf->area());
			nexp.plus(fexp);
		}
		if (this->has("alpha_term")) {
			nexp.plus(_alpha_term(pn));
		}
		nexp.plus(_f_term(pn));
		//
	}

	/**
	 * @brief get diffusion number
	 */
	void get_diffusion_number() const {
		ASSERT(this->has_time_term());

	}

	/**
	 * @brief override function solve
	 *        only one step, no time term
	 */
	int solve() {

		Mat mat;
		Arr b;

		Fun_get_spExp fun = [this](pNode pn) {
			spExp node_exp = this->_get_node_spexp(pn);
			this->_exp_substitute_ghost(pn, (*node_exp),this->_vars_c["phi"]);
			return node_exp;
		};

		this->_build_matrix(mat, b, fun);

		Arr x(b.size());  // initial x
		this->_grid_to_arr(x, this->_vars_c["phi"]);

		std::list<Float> lr;	//list residual
		//solver =======================
		int sf = Dia_BiCGSTAB(mat, x, b, this->_max_iter, this->_tol, lr);
		if (sf != 0) {
			std::cerr << " >! Poisson solve failed \n";
			//return -1;
		}
		//put the value back
		this->_arr_to_grid(x, this->_vars_c["phi"]);

		return 1;
	}

	void set_solve_para(vt tol = 1e-6, int max_iter = 1000) {
		this->_tol = tol;
		this->_max_iter = max_iter;
	}

	int initial() {
		if (this->is_stand_alone()) {
			this->_construct_stand_alone();
			// set time step
			if (this->has_time_term()) {
				this->_timestep->set_unknow_idx(this->_vars_c["phi"]);
				this->_timestep->set_idx(this->max_idx(this->_vars_c));
				//this->_timestep->set_dt(this->);
				// resize center data array for inner time
				st ss = this->_timestep->size_inner_step();
				this->_pdomain->resize_data(this->max_idx(this->_vars_c) + ss,
						0, 0, this->max_idx(this->_vars_ut) + 1);
			} else {
				this->_pdomain->resize_data(this->max_idx(this->_vars_c) + 1, 0,
						0, this->max_idx(this->_vars_ut) + 1);
			}
			this->init_beta();
			this->init_f();
			this->init_phi();
		} else {
			SHOULD_NOT_REACH;
		}

		this->_build_FE_on_leaf();
		this->_build_NE_on_leaf(1);
		this->_build_fun_get();

		std::function<void(pNode, Exp&)> f = std::bind(&Self::_node_exp, this,
				std::placeholders::_1, std::placeholders::_2);
		this->_build_node_exp(f, 0);  //exp on 0 of arr_spexp

		return -1;
	}

	int finalize() {
		this->_delete_FE_on_leaf();
		this->_delete_NE_on_leaf();
		return -1;
	}

	int run_one_step(st step) {
		std::cout << "poisson one step " << step << "\n";
		do {
			if (this->_timestep->do_solve()) {
				std::cout << "solve \n";
				// build a spexp on the node to solve
				Mat mat;
				Arr b;

				Fun_get_spExp fun =
						[this](pNode pn) {
							spExp node_exp(new Exp(*(this->_get_node_spexp(pn))));
							this->_exp_substitute_ghost(pn, (*node_exp),this->_vars_c["phi"]);

							spExp res = this->_timestep->new_exp(pn, node_exp);
							//this->_timestep.set_ex(pn, node_exp);
							return res;
						};

				this->_build_matrix(mat, b, fun);

				Arr x(b.size());  // initial x
				this->_grid_to_arr(x, this->_timestep->idx_old());

				std::list<Float> lr;	//list residual
				//solver =======================
				int sf = Dia_BiCGSTAB(mat, x, b, this->_max_iter, this->_tol,
						lr);
				if (sf != 0) {
					std::cerr << " >! solve failed \n";
					//return -1;
				}
				//put the value back
				this->_arr_to_grid(x, this->_timestep->idx_new());



			} else {
				std::cout << "solve explicit\n";

				Fun_pNode fun =
						[this](pNode pn) {

							spExp node_exp(new Exp(*(this->_get_node_spexp(pn))));
							this->_exp_substitute_ghost(pn, (*node_exp),this->_vars_c["phi"]);

							// coe = tau * volume / dt
							vt coe = this->_timestep->tau() * pn->volume() / this->_timestep->dt();
							st idx_n = this->_timestep->idx_new();
							st idx_o = this->_timestep->idx_old();

							vt v_o = pn->cd(idx_o);
							vt vex = node_exp->substitute(idx_o);

							vt vexp = v_o + node_exp->substitute(idx_o)/coe;
							// set value back to node
							pn->cd(idx_n) = vexp;
							if(pn->is_in_on(__X__,__Y__)) {
								std::cout<<"---------------------------\n";
								node_exp->show();
								std::cout<<"idx n   "<<idx_n<<std::endl;
								std::cout<<"idx o   "<<idx_o<<std::endl;

								std::cout<<"coe   = "<< coe<<std::endl;
								std::cout<<"vex   = "<< vex<<std::endl;
								std::cout<<"v o   = "<< v_o<<std::endl;
								std::cout<<"v n   = "<< pn->cd(idx_n)<<std::endl;
							}

						};

				pGrid pgrid = this->_pdomain->pgrid();

				pgrid->for_each_leaf(fun);

			}
			this->_timestep->inner_advance();  // advance inner step
		} while (!this->_timestep->is_end());
		// copy new to old
		Fun_pNode fun = [this](pNode pn) {
			vt idx_b = this->_timestep->idx_begin();
			vt idx_e = this->_timestep->idx_end();
			if(pn->is_in_on(__X__,__Y__)) {
				std::cout<<"idx b "<<idx_b<<" idx e "<<idx_e<<std::endl;
				std::cout<<"idx b "<<idx_b<<" idx e "<<idx_e<<std::endl;
				std::cout<<"val b "<< pn->cd(idx_b)<<std::endl;
			}
			pn->cd(idx_b) = pn->cd(idx_e);
			if(pn->is_in_on(__X__,__Y__)) {
				std::cout<<"val b "<< pn->cd(idx_b)<<std::endl;
			}
		};

		pGrid pgrid = this->_pdomain->pgrid();
		pgrid->for_each_leaf(fun);

		this->_timestep->next_step();

		return -1;
	}

	/**
	 *  @brief set alpha term
	 *         default --> no alpha term
	 *         default --> uniform alpha
	 */
	void set_alpha_term(Function pfun) {
		spEvent pse(new Flag("alpha_term", 1));
		this->_events["alpha_term"] = pse;
		this->set_val(pfun, this->_vars_c["alpha"]);
	}
	void unset_alpha_term() {
		this->_events.erase("alpha_term");
	}
	bool has_alpha_term() const {
		return this->has("alpha_term", 1);
	}
	/**
	 *  @brief set uniform beta
	 */
	void set_uniform_beta(vt beta = 1) {     //default
		spEvent pse(new Flag("uniform_beta", 1));
		this->_events["uniform_beta"] = pse;
		this->_values["uniform_beta"] = beta;
	}
	void unset_uniform_beta() {
		this->_events.erase("uniform_beta");
		this->_values.erase("uniform_beta");
	}

	/**
	 *  @brief set boundary condition phi
	 */
	void set_bc_phi(st si, st segi, pBoundaryCondition pbc) {
		this->set_boundary_condition(si, segi, this->_vars_c["phi"], pbc);
	}

	/**
	 * @brief set value on center of the node
	 */
	void set_beta(Function fun) {
		this->unset_uniform_beta();

		//
		this->_functions["set_beta"] = fun;

		//this->set_val(fun, this->_vars_c["beta"]);
	}

	void set_f(Function fun) {
		this->_functions["set_f"] = fun;
		//this->set_val(fun, this->_vars_c["f"]);
	}

	void set_phi(Function fun) {
		this->_functions["set_phi"] = fun;
	}

	void init_beta() {
		if (this->has("uniform_beta")) {
			return;
		}

		if (this->has_function("set_beta")) {
			Function f = this->_functions["set_beta"];
			this->set_val(f, this->_vars_c["beta"]);
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void init_f() {
		if (this->has_function("set_f")) {
			Function f = this->_functions["set_f"];
			this->set_val(f, this->_vars_c["f"]);
		} else {
			SHOULD_NOT_REACH;
		}
	}

	void init_phi() {
		if (this->has_function("set_phi")) {
			Function f = this->_functions["set_phi"];
			this->set_val(f, this->_vars_c["phi"]);
		} else {
			SHOULD_NOT_REACH;
		}
	}

	st phi_idx() const {
		return this->_vars_c.at("phi");
	}

protected:
	/*
	 * stand alone constructor
	 * constructor a Poisson class just used for solver only one poisson equation
	 */
	void _construct_stand_alone() {
		// variables on the center of node
		if (!this->has("uniform_beta")) {
			this->_vars_c["phi"] = 0;
			this->_vars_c["beta"] = 1;
			this->_vars_c["f"] = 2;
		} else {
			this->_vars_c["phi"] = 0;
			this->_vars_c["f"] = 1;
		}
		// untype variables on the node
		this->_vars_ut["FE"] = 0;
		this->_vars_ut["NE"] = 1;
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

	void _build_fun_get() {
		this->_get_node_spexp = [this](pNode pn) {
			ASSERT(pn!=nullptr);
			utPointer utp = pn->utp(this->_vars_ut["NE"]);
			Arr_spExp& arr_spexp = CAST_REF(Arr_spExp*, utp);
			return arr_spexp[0];
		};

	}

	/*
	 * set
	 */
}
;
}

#endif
