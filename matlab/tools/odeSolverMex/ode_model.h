/*
 *  ode_model.h
 *  odeSolve
 *
 *  Created by Vassilios Stathopoulos on 28/10/2011.
 *  Copyright 2011 Computing Science. All rights reserved.
 *
 */

#ifndef __ODE_MODEL_H__
#define __ODE_MODEL_H__

#ifdef __cplusplus  /* wrapper to enable C++ usage */
extern "C" {
#endif

#include <cvodes/cvodes.h>
#include <cvodes/cvodes_dense.h>
#include <nvector/nvector_serial.h>
#include <sundials/sundials_types.h>
#include <sundials/sundials_math.h>

#define	ODE_SOLVER_REL_ERR		RCONST(1.0e-6)
#define	ODE_SOLVER_ABS_ERR		RCONST(1.0e-7)
#define	ODE_SOLVER_MX_STEPS		50000
	
typedef int (*rhs_f)(realtype, N_Vector, N_Vector, void *);
typedef int (*jac_f)(int, realtype,
						N_Vector, N_Vector,
						DlsMat, void*,
						N_Vector, N_Vector, N_Vector);
	
typedef int (*rhs_sens)(int, realtype, N_Vector, N_Vector,
						int, N_Vector, N_Vector,
						void *, N_Vector, N_Vector);
	
typedef int (*func)(realtype, N_Vector, realtype *, void *);
typedef int (*func_sens)(realtype, N_Vector, N_Vector *, double *, void *);
	
typedef struct{
	int			N;				/* number of variables */
	int			P;				/* number of parameters */
	int			F;				/* number of functions */
	double*		v;				/* initial conditions */
	double*		p;				/* default parameters */
	rhs_f		vf_eval;		/* function pointer for ode RHS */
	jac_f		vf_jac;			/* function pointer for jacobian */
	rhs_sens	vf_sens;		/* function pointer for ode sensitivities */
	func		vf_func;		/* function pointer for functions RHS */
	func_sens	vf_func_sens;	/* function pointer for functions sensitivities */
	char**		v_names;		/* string array of variable names */
	char**		p_names;		/* string array of parameter names */
	char**		f_names;		/* string array of function names */
	void*		dylib;			/* pointer to dynamicaly linked library with the ode model */
} ode_model;
	

typedef struct{
	void*		cvode_mem;
	ode_model*	odeModel;
	N_Vector    y;
	N_Vector*	yS;
	double*		params;
	int			Psens;
}ode_solver;
		

/* Loads shared library with user defined functions and ode model data */
ode_model*	ode_model_loadFromFile(const char *filename);

/* inlined functions to allow access to the internals of the ode_model structure */
extern inline double ode_model_getN(const ode_model* model);

extern inline double ode_model_getP(const ode_model* model);

extern inline double ode_model_getF(const ode_model* model);

extern inline int ode_model_has_sens(const ode_model* model);

extern inline int ode_model_has_funcs(const ode_model* model);

extern inline int ode_model_has_funcs_sens(const ode_model* model);

extern inline const char** ode_model_get_var_names(const ode_model* model);

extern inline const char** ode_model_get_param_names(const ode_model* model);

extern inline const char** ode_model_get_func_names(const ode_model* model);

extern inline void ode_solver_disable_sens(ode_solver* solver);


/* returns initial conditions in y0 */
void	ode_model_get_initial_conditions(const ode_model* model, double* y0, int lenY0);
/* returns default parameter values in params */
void	ode_model_get_default_params(const ode_model* model, double* params, int lenParams);

	
/* Frees resources used by ode_model */
void		ode_model_free(ode_model* model);

/* Creates a new ode_solver for the given model */
ode_solver*	ode_solver_alloc(ode_model* model);
	
/* Initialises the ode solver. y0, yS0 and p can be NULL.
 */
void		ode_solver_init(ode_solver* solver, const double t0, double* y0, int lenY, double* p, int lenP );
void		ode_solver_init_sens(ode_solver* solver,  double* yS0, int lenP, int lenY);

/* Sets error tollerances, NOTE this function should be called only after an ode_solver_init or ode_solver_reinit.
*/
void		ode_solver_setErrTol(ode_solver* solver, const double rel_tol, double* abs_tol, const int abs_tol_len);
	
/* Re-initialises the ode solver. y0, yS0 and p can be NULL in which case the default values are used.
 */	
void		ode_solver_reinit(ode_solver* solver, const double t0,  double* y0, int lenY, const double* p, int lenP );
void		ode_solver_reinit_sens(ode_solver* solver, double* yS0, int lenP, int lenY);
	
/* Solves the ode system until time t and rerurns solution in y.
 */	
int		ode_solver_solve(ode_solver* solver, const double t, double* y, double* tout);

/* Solves the ode system until a steady state solution is found or t is reached.
*/
int		ode_solver_solveSteadyState(ode_solver* solver, const long int maxIter, double* y, double* tout);
    
/* Returns sensitivities for the current solution at t. This function must be called only after ode_solver_solve.
*/	
void		ode_solver_get_sens(ode_solver* solver, double t, double* yS);
/* Returns functions of the carrent solution.
*/	
void		ode_solver_get_func(ode_solver* solver, const double t, double* y, double* fy);
/* Returns sensitivities of functions of the carrent solution.
*/	
void		ode_solver_get_func_sens(ode_solver* solver, const double t, double* y, double* yS, double* fyS);

	
/* Prints to out solver statistics.
*/	
void		ode_solver_print_stats(const ode_solver* solver, FILE* outF);
	
/* Frees resources used by the ode_solver */
void		ode_solver_free(ode_solver* solver);

	
#ifdef __cplusplus
}
#endif

#endif