/*
 *  odeSolverMex.c
 *  sde_mcmc
 *
 *  Created by Vassilios Stathopoulos on 18/01/2013.
 *  Copyright 2013 Computing Science. All rights reserved.
 *
 *	[Y dY/dP] = odeSolverMex('model.so',[timePoints],[params],[y0],[dy0],[tol],[t0],[initSens])
 *				
 *	Y			- [N,tps] array of doubles. Solution of the ode.
 *	dY/dP		- [N*P,tps] array of doubles]. Sensitivities of Y w.r.t. parameters.
 *	'model.so'	- char array. Name of the ode model shared library.
 *	timePoints  - [tps,1] vector of doubles. Time points for obtaining solutions. If empty the system is solved until a steady state.
 *	params		- [P,1] vector of doubles. Parameters for the ode system. [Optional] default 
 *				   values are specified in the model shared library.
 *	y0			- [N,1] vector of doubles. Initial conditions. [Optional] default 
 *				   values are specified in the model shared library.
 *  dY0			- [N,P] vector of doubles. Initial conditions for sensitivities. 
 *				  [Optional] default is 0.
 *	tol			- [2 1] vector of doubles. Relative and absolute tolerance for the solver. 
 *				  [Optional] default values relTol = 1.e-5 absTol = 1.e-7.
 *	t0			-  scalar double. Integration starting time. [Optional] default value t0 = 0.
 *
 */

#include "mex.h"
#include "matrix.h"
#include "ode_model.h"

void printUsage(){
    mexPrintf("Usage: [Y Ys] = odeSolverMex(ode_system, [time_points], [params],[Y0],[dy0],[tolerances],[t0]). Arguments in [] are optional.\n");
    mexPrintf("Y			- [N,ntps] array of doubles. Solution of the ode.\n");
    mexPrintf("Ys           - [N*P,ntps] array of doubles]. Sensitivities of Y w.r.t. parameters.\n");
    mexPrintf("ode_system	- char array. Name of the ode model shared library.\n");
    mexPrintf("time_points	- [ntps,1] vector of doubles. Time points for obtaining solutions. [Optional] If empty the system is solved until a steady state.\n");
    mexPrintf("params       - [P,1] vector of doubles. Parameters for the ode system. [Optional] default values are specified in the model shared library.\n");
    mexPrintf("Y0           - [N,1] [N,1] vector of doubles. Initial conditions. [Optional] default values are specified in the model shared library.\n");
    mexPrintf("dY0          - [N,P] vector of doubles. Initial conditions for sensitivity equations. [Optional] default is 0.\n");
    mexPrintf("             - To estimate sensitivities w.r.t to initial conditions concatenate dY0 with an NxN identity matrix, e.g. [dy0 eye(N)].\n");
    mexPrintf("             - In that case Ys will be an [N*N*P, ntps] array. \n");
    mexPrintf("tolerances   - [2 1] vector of doubles. Relative and absolute tolerance for the solver. [Optional] default values relTol = 1e-5, absTol = 1e-7.\n");
    
}

void checkParams(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	if ( nrhs > 7 || nrhs < 1 ) {
        printUsage();
		mexErrMsgTxt("Incorect number of input arguments");
	}
	if (nlhs > 2) {
        printUsage();
		mexErrMsgTxt("Too many output arguments.");
	}
	if ( !mxIsChar(prhs[0]) ) {
        printUsage();
		mexErrMsgTxt("First argument must be a char array with the file name of the ode system shared library created by vfgen.");
	}
    
	if ( nrhs > 1 && !mxIsEmpty(prhs[1]) ) {
        int ntps = mxGetN(prhs[1]);
        if ( !mxIsDouble(prhs[1]) || ntps != 1 ){
            printUsage();
            mexErrMsgTxt("Second argument must be a vector of doubles with time points.");
        }
        double* tps = mxGetPr(prhs[1]);
        int i;
        for ( i = 0; i < ntps-1; i++ ) {
            if ( tps[i] >= tps[i+1] ){
                printUsage();
                mexErrMsgTxt("Second argument must be a strictly increasing vector of doubles with time points.");
            }
        }
    }
    
	if ( nrhs > 2 && !mxIsEmpty(prhs[2]) ){
		if ( !mxIsDouble(prhs[2]) || mxGetN(prhs[2]) != 1 ){
            printUsage();
			mexErrMsgTxt("Third argument must a vector of doubles with the parameters.");
		}
	}
	if ( nrhs > 3 && !mxIsEmpty(prhs[3]) ){
		if ( !mxIsDouble(prhs[3]) || mxGetN(prhs[3]) != 1 ){
            printUsage();
			mexErrMsgTxt("Fourth argument must a vector of doubles with the initial conditions.");
		}
	}
	if ( nrhs > 4 && !mxIsEmpty(prhs[4]) ){
		if ( !mxIsDouble(prhs[4]) ){
            printUsage();
			mexErrMsgTxt("Fifth argument must a matrix of doubles with the initial conditions for the sensitivities.");
		}
	}
	
	if ( nrhs > 5 && !mxIsEmpty(prhs[5]) ){
		if ( !mxIsDouble(prhs[5]) || mxGetN(prhs[5]) != 1 || mxGetM(prhs[5]) != 2 ){
            printUsage();
			mexErrMsgTxt("Sixth argument must a vector of 2 doubles with the relative and absolute tollerance.");
		}
	}
	if ( nrhs > 6 && !mxIsEmpty(prhs[6]) ) {
		if ( !mxIsDouble(prhs[6]) || mxGetN(prhs[6]) != 1 || mxGetM(prhs[6]) != 1 ){
			printUsage();
            mexErrMsgTxt("Seventh argument must a double scalar with the initial time of integration, i.e. t0.");
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]){
	
	checkParams(nlhs, plhs, nrhs, prhs);
	
	/* default values */
	double relTol = 1e-5;
	double absTol = 1e-7;
	double t0	  = 0;
	size_t ntps = 0;
	
	const mxArray* modelName_mx			= prhs[0];			/* char array */
    
    const mxArray* timePoints_mx = prhs[1];
    ntps = 1;
    if (nrhs > 1 && !mxIsEmpty(prhs[1])){
        /* time points */
        timePoints_mx		= prhs[1];                      /* [ntps 1] of doubles */
        ntps = mxGetM(timePoints_mx);
    }
    
	const mxArray* params_mx = NULL;
	if ( nrhs > 2 ){
		/* user parameters */
		params_mx						= prhs[2];			/* [P 1] of doubles */
	}
	const mxArray* initY0_mx = NULL;
	if ( nrhs > 3 ) {
		/* user initial conditions */
		initY0_mx						= prhs[3];			/* [N 1] of doubles */
	}
	const mxArray* initDY0_mx = NULL;
	if (nrhs > 4) {
		/* user initial conditions */
		initDY0_mx						= prhs[4];			/* [N P] of doubles */
	}
	if ( nrhs > 5 && !mxIsEmpty(prhs[5]) ) {
		/* user tolerances */
		const mxArray* tol_mx			= prhs[5];			/* [2 1] of doubles */
		double* tol = mxGetPr(tol_mx);
		relTol = tol[0];
		absTol = tol[1];
	}
	if ( nrhs > 6 && !mxIsEmpty(prhs[6]) ) {
		/* user t0 */
		const mxArray* t0_mx			= prhs[6];			/* [1 1] of doubles */
		double* t0_c = mxGetPr(t0_mx);
		t0 = t0_c[0];
	}
	
	/* load shared library model */
	const char* modelName = mxArrayToString(modelName_mx);
	ode_model* odeSystem = ode_model_loadFromFile(modelName);
	if(odeSystem == NULL){
		/* memory allocated with mxMalloc is cleaned by mexErrMsgTxt */
		mexErrMsgTxt("Could not load ode model shared library.");
	}
	mxFree(modelName);
	
	int N = ode_model_getN(odeSystem);
	int P = ode_model_getP(odeSystem);
	int Psens = P;
	
    
	double* theta = NULL;
	if ( nrhs > 2 && !mxIsEmpty(params_mx) ){
		if ( mxGetM(params_mx)!= P ){
			ode_model_free(odeSystem);
			mexErrMsgTxt("Parameters should be of equal size with the parametes in the ode system.");
		}
		theta = mxGetPr(params_mx);
	}else {
		theta = mxMalloc(sizeof(double)*P);
		ode_model_get_default_params(odeSystem, theta, P);
	}
	
	double* Y0 = NULL;
	if ( nrhs > 3 && !mxIsEmpty(initY0_mx) ){
		if ( mxGetM(initY0_mx)!= N ){
			ode_model_free(odeSystem);
			mexErrMsgTxt("Initial conditions should be of equal size with the ode system equations.");
		}
		Y0 = mxGetPr(initY0_mx);
	}else{
		Y0 = mxMalloc(sizeof(double)*N);
		ode_model_get_initial_conditions(odeSystem, Y0, N);
	}
	double* dY0 = NULL;
	if (nrhs > 4 && !mxIsEmpty(initDY0_mx)) {
		if ( mxGetM(initDY0_mx)!= N ){
			ode_model_free(odeSystem);
			mexErrMsgTxt("Sensitivity initial conditions should be of equal size with the ode system equations.");
		}
		dY0 = mxGetPr(initDY0_mx);
		Psens = mxGetN(initDY0_mx);
	}

	/* allocate memmory for outputs */
	mxArray* Y_mx = mxCreateDoubleMatrix(0, 0, mxREAL);
	mxArray* Ys_mx = NULL;
	mxSetM(Y_mx, N);
	mxSetN(Y_mx, ntps);
	mxSetData(Y_mx, mxMalloc(sizeof(double)*N*ntps));
	if ( ode_model_has_sens(odeSystem) &&  nlhs > 1 ){
		Ys_mx = mxCreateDoubleMatrix(0, 0, mxREAL);
		mxSetM(Ys_mx, N*Psens);
		mxSetN(Ys_mx, ntps);
		mxSetData(Ys_mx, mxMalloc(sizeof(double)*N*Psens*ntps));
	}
	
	/* TODO: error checking */
	ode_solver* solver = ode_solver_alloc(odeSystem);
	ode_solver_init(solver, t0, Y0, N, theta, P);
	ode_solver_setErrTol(solver, relTol, &absTol, 1);
	if ( ode_model_has_sens(odeSystem) && nlhs > 1 ) {
		if (dY0 != NULL)
			ode_solver_init_sens(solver, dY0, Psens, N);
		else
			ode_solver_init_sens(solver, NULL, 0, 0);
	}

    double* Y = mxGetPr(Y_mx);
    
    if (nrhs > 1 && !mxIsEmpty(timePoints_mx)){
        
        double* tps = mxGetPr(timePoints_mx);
        
        double t_new, tout;
        int i;
        for ( i = 0; i < ntps; i++ ) {
            t_new = tps[i];

            int CVerror =  ode_solver_solve(solver, t_new, &Y[N*i], &tout);
            if ( CVerror ) {
                mexPrintf("ODE solver failed.\n");
			
                /* clean up outputs */
                /*mxFree(Y_mx);
                 if ( ode_model_has_sens(odeSystem) &&  nlhs > 1 ){
                 mxFree(Ys_mx);
                 }*/
                break;
            }
		
            if ( ode_model_has_sens(odeSystem) && nlhs > 1 ){
                double* Ys = mxGetPr(Ys_mx);
                ode_solver_get_sens(solver, tout, &Ys[N*Psens*i]);
            }
        }
    }
    else{
        //default maxT = 500.0;
        long int maxIter = 10000;
        double tout = 0.0;
        int CVerror =  ode_solver_solveSteadyState(solver, maxIter, &Y[0], &tout);
        if ( CVerror ){
             mexPrintf("ODE solver failed to reach a steady state.\n");
        }
        else
        {
            if ( ode_model_has_sens(odeSystem) && nlhs > 1 ){
                double* Ys = mxGetPr(Ys_mx);
                ode_solver_get_sens(solver, tout, &Ys[0]);
            }
        }
    }
    
	/* Set outputs */
	if (Y_mx != NULL) {
		plhs[0] =  Y_mx;		/* [N ntps] of doubles */
	}
	if ( nlhs > 1 && Ys_mx != NULL){
		plhs[1] = Ys_mx;		/* [N*p ntps] of doubles */
	}
	
	/* clean up */
	if ( (nrhs < 3) || (mxIsEmpty(params_mx)) ){
		mxFree(theta);
	}
	if ( (nrhs < 4) || (mxIsEmpty(initY0_mx)) ){
		mxFree(Y0);
	}
	
	ode_solver_free(solver);
	ode_model_free(odeSystem);
}

