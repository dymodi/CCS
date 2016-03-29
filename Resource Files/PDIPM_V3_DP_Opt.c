 // *
 // * PDIPM_V3_DP_Opt.c
 // *
 // *  Created on: 2014-8-14
 // *      Author: Yi
 // *	The 4th version of QP solving routines.
 // *
 // *	Algorithm:		Predictor Corrector Primal Dual Interior-Point Method
 // *	Gobal Check:	Yes
 // *	Starting Point:	Determined in the outer MPC function
 // *	Step Length:	Accurate formulation (The subfunction alpha_decreas() )
 // *	LS Solving:		Cholesky factorization
 // *	Termination:	Himmelblau Termination Criteria
 // *  

#include <xdc/std.h>

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <xdc/runtime/Timestamp.h>
#include "mpc_lib_DP.h"

extern double Time_LineSolve1,Time_LineSolve2;
extern double Time_Section1,Time_Section2,Time_Section3,Time_Section4,Time_Section5,Time_Section6,Time_Section7,Time_Section8;
extern double QP_Iteration;

extern int MaxIteration;

extern double *x_old, *rd, *rp, *rp_y,*inv_y_lambda, *y_lambda_A,*equation_F, *equation_b;
extern double *delta_x, *delta_y, *delta_lambda,*QP_temp_cond,*x_try,*sig_mu, *center_part,*QP_temp1, *QP_temp2, *QP_temp3;
extern double *c_rev;

extern double *line_solve_L,*line_solve_p,*equation_F_scaling,*equation_b_scaling;

extern double *L,*B,*G,*GL,*GL_T,*F,*Phi,*OMEGA_L,*N_OMEGA_L,*N_OMEGA_L_T;
extern double *delta_x_aff,*delta_y_aff,*delta_lambda_aff,*sig_mu_lam_y,*aff_part;

extern double b_max,c_max,G_max,A_max,p_max;

//V4 modification:
//	1 Himmelblau Termination Criteria
//	2 Initial Value Reset
void priduip_v4_DP(double *G, double *GL,double *GL_T, double *c, double *A, double *A_t, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u)
{

	double Fx_old, Fx, cons_obey;
	double res = 0.1;
	double tau = 0.9;
	int Iter = 0;
	int k, k_in;
	double mu,mu_aff;
	double alpha_aff_pri = 0.999995;
	double alpha_aff_dual = 0.999995;
	double alpha_aff = 0.99995;
	double alpha_cond;	//Use to check whether the alpha is valid.
	double sigma;
	double alpha_tau_pri = 0.999995;
	double alpha_tau_dual = 0.999995;
	double alpha;
	double res_1,res_2;
	double tau_part;
	int flag = 0;
	double scale_F, scale_b, scale_linesolve;

	UInt32 start_time,end_time,dead_time,total_time;
	UInt32 start_time_ls1,end_time_ls1,total_time_ls1;
	UInt32 start_time_ls2,end_time_ls2,total_time_ls2;
	UInt32 start_time_section,end_time_section,total_time_section;

	start_time = Timestamp_get32();
	end_time = Timestamp_get32();
	dead_time = end_time - start_time;

	//Show the problem paremeters.
	//show_matrix(G,n,n);putchar('\n');show_matrix(c,1,n);putchar('\n');
	//show_matrix(A,n,m);putchar('\n');show_matrix(b,1,m);putchar('\n');
	//show_matrix(x,1,n);putchar('\n');show_matrix(y,1,m);putchar('\n');show_matrix(lambda,1,m);putchar('\n');
	
	//tau preparation
	//tau_part = 0.099999/MaxIteration;

	//Initial Fx_old value
	mat_vec(G,x,n,n,QP_temp1);
	Fx_old = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
	//printf("This initial F(x) value is: %f\n",Fx_old);

	//Initial x_old value
	#pragma MUST_ITERATE (3);
	for(k=0;k<n;k++)
		x_old[k] = x[k];

	//Check if the Uncontrained solution meet the constraints
	//vec_rev(c,n,c_rev);
	luEvaluate(GL,GL_T,c_rev,n,x_try);
	mat_vec(A,x_try,m,n,QP_temp1);
	vec_sub(QP_temp1,b,m,QP_temp2);
	if((vec_min(QP_temp2,m))>0)
	{
		memcpy(delta_u,x_try,n*sizeof(double));
		QP_Iteration = Iter;
	}
	else
	{
		// If not global optimal, do iterative process.
		for(k=1;k<=MaxIteration;k++)
		{
			start_time_section = Timestamp_get32();

			mat_vec(G,x,n,n,QP_temp1);		
			mat_vec(A_t,lambda,n,m,QP_temp2);	
			vec_sub(QP_temp1, QP_temp2, n,QP_temp3);
			vec_add(QP_temp3,c,n,rd);
			//printf("rd is:\n");show_matrix(rd,1,n);putchar('\n');

			mat_vec(A,x,m,n,QP_temp1);
			vec_sub(b,QP_temp1,m,rp_y);
			vec_sub(QP_temp1,y,m,QP_temp2);
			vec_sub(QP_temp2,b,m,rp);
			//printf("rp is:\n");show_matrix(rp,1,m);putchar('\n');

			vec_div(lambda,y,m,inv_y_lambda);
			bsxfun(A,inv_y_lambda,m,n,y_lambda_A);
			//printf("inv_y_lambda is :\n");show_matrix(inv_y_lambda,1,m);putchar('\n');
			//printf("y_lambda_A is :\n");show_matrix(y_lambda_A,n,m);putchar('\n');

			//end_time_section = Timestamp_get32();
			////total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section1 +=(double)total_time_section;


			//start_time_section = Timestamp_get32();
			mat_mul(A_t,n,m,y_lambda_A,n,QP_temp1);	//Time Consuming! Much better if custom!
			mat_add(G,QP_temp1,n,n,equation_F);
			//printf("equation_F is :\n");show_matrix(equation_F,n,n);putchar('\n');

			vec_mul(inv_y_lambda,rp_y,m,QP_temp1);
			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
			vec_sub(QP_temp2,rd,n,equation_b);
			//printf("equation_b is: \n");show_matrix(equation_b,1,n);

			////Line_solve scaling
			//scale_F = scalar_min(vec_min(equation_F,n*n),vec_max(equation_F,n*n));
			//scale_b = scalar_min(vec_min(equation_b,n),vec_max(equation_b,n));
			//scale_F = scalar_abs(scale_F);
			//scale_b = scalar_abs(scale_b);
			//scale_linesolve = scalar_min(scale_F,scale_b);
			//for(k_in=0;k_in<n*n;k_in++)
			//	equation_F_scaling[k_in] = equation_F[k_in]/scale_linesolve;
			//for(k_in=0;k_in<n;k_in++)
			//	equation_b_scaling[k_in] = equation_b[k_in]/scale_linesolve;
			//printf("scale_linesolve is:%f \n",scale_linesolve);
			//printf("equation_F_scaling is: \n");show_matrix(equation_F_scaling,n,n);putchar('\n');
			//printf("equation_b_scaling is: \n");show_matrix(equation_b_scaling,1,n);putchar('\n');

			
			//start_time_ls1 = Timestamp_get32();
			//printf("equation_F is: \n");show_matrix(equation_F,n,n);putchar('\n');
			//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');
			//line_solve(equation_F,equation_b,n,delta_x_aff);
			memcpy(line_solve_L,equation_F,n*n*sizeof(double));
			chol_NRC(line_solve_L,n,line_solve_p);
			luEvaluate_NRC(line_solve_L,equation_b,line_solve_p,n,delta_x_aff);
			
			//end_time_ls1 = Timestamp_get32();
			//total_time_ls1 = end_time_ls1 - start_time_ls1 - dead_time;
			//Time_LineSolve1+=(double)total_time_ls1;
			//printf("delta_x_aff is: \n");show_matrix(delta_x_aff,1,n);putchar('\n');

			mat_vec(A,delta_x_aff,m,n,QP_temp1);
			vec_add(QP_temp1,rp,m,delta_y_aff);
			//printf("delta_y_aff is: \n");show_matrix(delta_y_aff,1,m);putchar('\n');

			vec_mul(inv_y_lambda,delta_y_aff,m,QP_temp1);
			vec_add(QP_temp1,lambda,m,QP_temp2);
			vec_rev(QP_temp2,m,delta_lambda_aff);
			//printf("delta_lambda_aff is: \n");show_matrix(delta_lambda_aff,1,m);putchar('\n');

			//end_time_section = Timestamp_get32();
			//total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section2 +=(double)total_time_section;

			//start_time_section = Timestamp_get32();

			mu = dot_product(y,lambda,m)/m;
			//printf("mu is: %f \n",mu);

			alpha_aff_pri = alpha_decreas(y,delta_y_aff,m,1);
			alpha_aff_dual = alpha_decreas(lambda,delta_lambda_aff,m,1);
			//alpha_aff_pri = alpha_decreas_old(y,delta_y_aff,m,1);
			//alpha_aff_dual = alpha_decreas_old(lambda,delta_lambda_aff,m,1);
			alpha_aff = scalar_min(alpha_aff_pri,alpha_aff_dual);
			//printf("alpha_aff is : %f\n",alpha_aff);

			//end_time_section = Timestamp_get32();
			//total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section3 +=(double)total_time_section;

			//start_time_section = Timestamp_get32();

			sca_vec_mut(alpha_aff,delta_y_aff,m,QP_temp3);
			vec_add(y,QP_temp3,m,QP_temp1);
			sca_vec_mut(alpha_aff,delta_lambda_aff,m,QP_temp3);
			vec_add(lambda,QP_temp3,m,QP_temp2);
			mu_aff = dot_product(QP_temp1,QP_temp2,m)/m;

			sigma = mu_aff/mu;
			sigma = sigma*sigma*sigma;

			//start_time = clock();

			vec_mul(delta_lambda_aff,delta_y_aff,m,QP_temp1);
			sca_vec_sub(sigma*mu,QP_temp1,m,sig_mu_lam_y);
			vec_div(sig_mu_lam_y,lambda,m,QP_temp2);
			vec_add(rp_y,QP_temp2,m,aff_part);
			vec_mul(inv_y_lambda,aff_part,m,QP_temp1);
			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
			vec_sub(QP_temp2,rd,n,equation_b);
			//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');

			////Line_solve scaling
			//scale_b = scalar_min(vec_min(equation_b,n),vec_max(equation_b,n));
			//scale_b = scalar_abs(scale_b);
			//scale_linesolve = scalar_min(scale_F,scale_b);
			//for(k_in=0;k_in<n*n;k_in++)
			//	equation_F_scaling[k_in] = equation_F[k_in]/scale_linesolve;
			//for(k_in=0;k_in<n;k_in++)
			//	equation_b_scaling[k_in] = equation_b[k_in]/scale_linesolve;
			//printf("scale_linesolve is:%f \n",scale_linesolve);
			//printf("equation_F_scaling is: \n");show_matrix(equation_F_scaling,n,n);putchar('\n');
			//printf("equation_b_scaling is: \n");show_matrix(equation_b_scaling,1,n);putchar('\n');
			
			//start_time_ls2 = clock();
			//line_solve(equation_F,equation_b,n,delta_x);
			memcpy(line_solve_L,equation_F,n*n*sizeof(double));
			chol_NRC(line_solve_L,n,line_solve_p);
			luEvaluate_NRC(line_solve_L,equation_b,line_solve_p,n,delta_x);
			
			//end_time_ls2 = clock();
			//total_time_ls2 = end_time_ls2 - start_time_ls2 - dead_time;
			//Time_LineSolve2+=(double)total_time_ls2;

			//printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');

			mat_vec(A,delta_x,m,n,QP_temp1);
			vec_add(QP_temp1,rp,m,delta_y);
			//printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');

			vec_mul(lambda,delta_y,m,QP_temp1);
			vec_mul(delta_lambda_aff,delta_y_aff,m,QP_temp2);
			vec_add(QP_temp1,QP_temp2,m,QP_temp3);
			sca_vec_sub(sigma*mu,QP_temp3,m,QP_temp1);
			vec_div(QP_temp1,y,m,QP_temp2);
			vec_sub(QP_temp2,lambda,m,delta_lambda);
			//printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');

			//end_time_section = Timestamp_get32();
			//total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section4 +=(double)total_time_section;

			//start_time_section = Timestamp_get32();

			if(res<0.1)
				tau = 1-res;
			//tau = tau + tau_part * Iter;
			//printf("tau is %f \n\n",tau);
			
			//end_time = Timestamp_get32();
			//total_time = end_time - start_time - dead_time;
			//Time_QP_step3[0]+=(double)total_time;

			//start_time = clock();

			alpha_tau_pri = alpha_decreas(y,delta_y,m,tau);
			alpha_tau_dual = alpha_decreas(lambda,delta_lambda,m,tau);
			//alpha_tau_pri = alpha_decreas_old(y,delta_y,m,tau);
			//alpha_tau_dual = alpha_decreas_old(lambda,delta_lambda,m,tau);
			alpha = scalar_min(alpha_tau_pri,alpha_tau_dual);
			//printf("alpha is : %f\n",alpha);

			//end_time_section = Timestamp_get32();
			//total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section5 +=(double)total_time_section;

			//start_time_section = Timestamp_get32();

			//Obtain new iteration value.
			sca_vec_mut(alpha,delta_x,n,QP_temp1);
			vec_add(x,QP_temp1,n,x);
			sca_vec_mut(alpha,delta_y,m,QP_temp1);
			vec_add(y,QP_temp1,m,y);
			sca_vec_mut(alpha,delta_lambda,m,QP_temp1);
			vec_add(lambda,QP_temp1,m,lambda);
			//printf("x is :\n");show_matrix(x,1,n);putchar('\n');
			//printf("y is :\n");show_matrix(y,1,m);putchar('\n');
			//printf("lambda is :\n");show_matrix(lambda,1,m);putchar('\n');

			// Iteration number + 1
			Iter = Iter + 1;
			//printf("Iteration: %d\n",Iter);

			//end_time_section = Timestamp_get32();
			//total_time_section = end_time_section - start_time_section - dead_time;
			//Time_Section6 +=(double)total_time_section;

			//start_time_section = Timestamp_get32();

			//// KKT residual (Wright) Termination Criteria
			//flag = TC_KKT_Wright(rd, rp, b_max,c_max, mu, 0.0001,m, n);
			//if(flag==1)
			//	break;

			//// Himmelblau Termination Criteria
			mat_vec(G,x,n,n,QP_temp1);
			Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
			res_1 = scalar_abs(Fx_old-Fx);
			vec_sub(x_old,x,n,QP_temp1);
			vec_abs(QP_temp1,n,QP_temp2);
			res_2 = vec_max(QP_temp2,n);
			mat_vec(A,x,m,n,QP_temp1);
			vec_sub(QP_temp1,b,m,QP_temp2);
			cons_obey = vec_min(QP_temp2,m);
			res = scalar_max(res_1,res_2);
			//printf("The Fx_old is: %f \n",Fx_old);
			//printf("The Fx is: %f \n",Fx);
			// Update
			Fx_old = Fx;
			for(k_in=0;k_in<n;k_in++)
				x_old[k_in] = x[k_in];
			//printf("The res_1 is: %f \n",res_1);
			//printf("The res_2 is: %f \n",res_2);
			//printf("cons_obey is: %f \n",cons_obey);
			//Terminal condition(Himmelblau)
			if(res<0.0001&&cons_obey>=-0.0001)
				break;


			//end_time = clock();
			//total_time = end_time - start_time - dead_time;
			//Time_QP_step4[0]+=(double)total_time;

			//// Show new Fx value
			//mat_mul(x,1,n,G,n,QP_temp1);
			//Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
			//printf("This iteration F(x) value is: %f\n\n",Fx);

			////Terminal condition(KKT Residuel)
			//if(res<0.00001)
			//	break;


		}
		memcpy(delta_u,x,n*sizeof(double));
		QP_Iteration = Iter;

		//end_time_section = Timestamp_get32();
		//total_time_section = end_time_section - start_time_section - dead_time;
		//Time_Section7 +=(double)total_time_section;
	}
}

//================================priduip end======================================================

