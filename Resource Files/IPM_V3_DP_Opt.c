 // *
 // * IPM_V3_DP.c
 // *
 // *  Created on: 2014-7-20
 // *      Author: Yi
 // *
 // *	ALL memory processing put in the main function.
 // *	ALL variables are double precision.
 // *
 // *	Algorithm:		Primal Dual Interior-Point Method
 // *	Center path:	Yes
 // *	Pre.-Cor.:		None
 // *	Gobal Check:	Yes
 // *	Starting Point:	Determined in the outer MPC function
 // *	Step Length:	AUT method (can choose While loop form)
 // *	LS Solving:		Cholesky factorization
 // *	Termination:	Himmelblau Termination Criteria
 // *
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
extern int QP_Iteration;
extern int LS_solve_times;
//extern int alpha_compute_times;
extern double Time_QP_non_Iter,Time_QP_Iter;

extern int MaxIteration;

extern double *x_old, *rd, *rp, *rp_y,*inv_y_lambda, *y_lambda_A,*equation_F, *equation_b;
extern double *delta_x, *delta_y, *delta_lambda,*temp_cond,*x_try,*sig_mu, *center_part,*QP_temp1, *QP_temp2, *QP_temp3;
extern double *c_rev;

extern double *line_solve_L,*line_solve_p,*line_solve_U,*line_solve_D,*line_solve_LD;

extern double *L,*B,*G,*GL,*GL_T,*F,*Phi,*OMEGA_L,*N_OMEGA_L,*N_OMEGA_L_T;

extern double b_max,c_max,G_max,A_max,p_max;

extern FILE *fp_Fx,*fp_cons_obey;


void IPM_V3_DP(double *G, double *GL,double *GL_T, double *c, double *A, double *A_t, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u)
{

	double Fx_old, Fx, cons_obey, cons_obey_old;
	int Iter = 0;
	double mu,mu_old;
	//double alpha_aff_pri = 0.999995;
	//double alpha_aff_dual = 0.999995;
	//double alpha_aff = 0.99995;
	double sigma = 0.1;
	double res = 0.3;
	//double tau = 0.7;
	double alpha_tau_pri = 0.999995;
	double alpha_tau_dual = 0.999995;
	double alpha;
	double res_1,res_2;

	int flag,CDC_flag;
	int k, k_in;
	int i,j;

	double epsilon = 0.00000001;
	//double epsilon = 0.1;

	UInt32 start_time,end_time,dead_time,total_time;
	UInt32 start_time_section,end_time_section,total_time_section;
	UInt32 QP_non_Iter_start_time,QP_non_Iter_end_time;
	UInt32 QP_Iter_start_time,QP_Iter_end_time;

	start_time = Timestamp_get32();
	end_time = Timestamp_get32();
	dead_time = end_time - start_time;

	//Show the problem paremeters.
	//printf("QP informations.\n");
	//show_matrix(G,n,n);putchar('\n');show_matrix(c,1,n);putchar('\n');
	//show_matrix(A,n,m);putchar('\n');show_matrix(b,1,m);putchar('\n');
	//show_matrix(x,1,n);putchar('\n');show_matrix(y,1,m);putchar('\n');show_matrix(lambda,1,m);putchar('\n');

	QP_non_Iter_start_time = Timestamp_get32();

	//Initial Fx_old value
	mat_vec(G,x,n,n,QP_temp1);
	Fx_old = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
	//printf("This initial F(x) value is: %f\n",Fx_old);

	//Initial cons_obey_old  value
	mat_vec(A,x,m,n,QP_temp1);
	vec_sub(b,QP_temp1,m,QP_temp2);
	cons_obey = vec_max(QP_temp2,m);

	mu = dot_product(y,lambda,m)/m;

	//Initial x_old value
	//#pragma MUST_ITERATE (3);
	for(k=0;k<n;k++)
		x_old[k] = x[k];

	//Check if the Uncontrained solution meet the constraints

	//printf("GL is:\n");show_matrix(GL,n,n);putchar('\n');
	//printf("GL_T is:\n");show_matrix(GL_T,n,n);putchar('\n');
	//printf("c_rev is:\n");show_matrix(c_rev,1,n);putchar('\n');
	luEvaluate(GL,GL_T,c_rev,n,x_try);
	//printf("x_try is:\n");show_matrix(x_try,1,n);putchar('\n');
	mat_vec(A,x_try,m,n,QP_temp1);
	//printf("QP_temp1 is:\n");show_matrix(QP_temp1,1,m);putchar('\n');
	vec_sub(QP_temp1,b,m,QP_temp2);
	//printf("QP_temp2 is:\n");show_matrix(QP_temp2,1,m);putchar('\n');

	QP_non_Iter_end_time = Timestamp_get32();
	Time_QP_non_Iter = (double)QP_non_Iter_end_time-QP_non_Iter_start_time-dead_time;

	QP_Iter_start_time = Timestamp_get32();

	//printf("(vec_min(QP_temp2,m):%f\n",vec_min(QP_temp2,m));
	if((vec_min(QP_temp2,m))>0)
	{
		memcpy(delta_u,x_try,n*sizeof(double));
		QP_Iteration = 0;
	}
	else
	{
		// If not global optimal, do iterative process.
		//printf("MaxIteration is : %d \n",MaxIteration);
		for(k=1;k<=MaxIteration;k++)
		{
			//start_time = Timestamp_get32();
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

			mat_mul(A_t,n,m,y_lambda_A,n,QP_temp1);	//Time Consuming! Much better if custom!
			vec_add(G,QP_temp1,n*n,equation_F);
			//printf("equation_F is: \n");show_matrix(equation_F,n,n);putchar('\n');

			sca2vec(sigma*mu,m,sig_mu);

			vec_div(sig_mu,lambda,m,QP_temp2);
			vec_add(rp_y,QP_temp2,m,center_part);
			vec_mul(inv_y_lambda,center_part,m,QP_temp1);
			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
			vec_sub(QP_temp2,rd,n,equation_b);
			//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');

			end_time_section = Timestamp_get32();
			total_time_section = end_time_section - start_time_section - dead_time;
			Time_Section1 += (double)total_time_section;

			start_time_section = Timestamp_get32();

			//////line_solve(equation_F,equation_b,n,delta_x);
			memcpy(line_solve_L,equation_F,n*n*sizeof(double));
			flag = chol_NRC(line_solve_L,n,line_solve_p);
			if(flag==0)
			{
				//Factorization Failed solution 1
				//Below is Modified Cholesky Factorization
				mcf(equation_F,n,line_solve_L,line_solve_D);
				memcpy(line_solve_U,line_solve_L,n*n*sizeof(double));
				transpose(line_solve_U, n, n);
				Lower_diag_vec(line_solve_L,line_solve_D,n,line_solve_LD);
				luEvaluate(line_solve_LD,line_solve_U,b,n,delta_x);

				//Factorization Failed solution 2
				//break;
			}
			else
				luEvaluate_NRC(line_solve_L,equation_b,line_solve_p,n,delta_x);

			//printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');

			mat_vec(A,delta_x,m,n,QP_temp1);
			vec_add(QP_temp1,rp,m,delta_y);
			//printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');

			vec_mul(lambda,delta_y,m,QP_temp2);
			vec_sub(QP_temp2,sig_mu,m,QP_temp1);
			vec_div(QP_temp1,y,m,QP_temp2);
			vec_rev(QP_temp2,m,QP_temp3);
			vec_sub(QP_temp3,lambda,m,delta_lambda);
			//printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');

			end_time_section = Timestamp_get32();
			total_time_section = end_time_section - start_time_section - dead_time;
			Time_Section2 +=(double)total_time_section;


			start_time_section = Timestamp_get32();

			////alpha_compute_times++;

			alpha_tau_pri = alpha_decreas(y,delta_y,m,1);
			alpha_tau_dual = alpha_decreas(lambda,delta_lambda,m,1);
			alpha = scalar_min(alpha_tau_pri,alpha_tau_dual);
			alpha = scalar_min(0.999999,alpha);
			//printf("alpha is : %f\n",alpha);

			//////alpha = alpha_decreas_while(y, lambda, delta_y, delta_lambda, m);

			end_time_section = Timestamp_get32();
			total_time_section = end_time_section - start_time_section - dead_time;
			Time_Section3 +=(double)total_time_section;

			start_time_section = Timestamp_get32();

			mu_old = mu;

			////Obtain new iteration value.
			sca_vec_mut(alpha,delta_x,n,QP_temp1);
			vec_add(x,QP_temp1,n,x);
			sca_vec_mut(alpha,delta_y,m,QP_temp1);
			vec_add(y,QP_temp1,m,y);
			sca_vec_mut(alpha,delta_lambda,m,QP_temp1);
			vec_add(lambda,QP_temp1,m,lambda);

			//printf("x is :\n");show_matrix(x,1,n);putchar('\n');
			//printf("y is :\n");show_matrix(y,1,m);putchar('\n');
			//printf("lambda is :\n");show_matrix(lambda,1,m);putchar('\n');

			//// Iteration number + 1
			Iter = Iter + 1;
			//printf("Iteration: %d\n",Iter);

			mu = dot_product(y,lambda,m)/m;

			sigma = mu/mu_old;
			sigma = sigma*sigma*sigma;
			sigma = scalar_min(sigma,0.99999);

			end_time_section = Timestamp_get32();
			total_time_section = end_time_section - start_time_section - dead_time;
			Time_Section4 +=(double)total_time_section;

			start_time_section = Timestamp_get32();

			//// Calculate KKT residual.
			//mat_vec(G,x,n,n,QP_temp1);
			//mat_vec(A_t,lambda,n,m,QP_temp3);
			//vec_sub(QP_temp3,c,n,QP_temp2);
			//vec_sub(QP_temp1,QP_temp2,n,QP_temp3);
			//res_1 = scalar_max(scalar_abs(vec_min(QP_temp3,n)),scalar_abs(vec_max(QP_temp3,n)));
			//mat_vec(A,x,m,n,QP_temp1);
			//vec_add(y,b,m,QP_temp2);
			//vec_sub(QP_temp1,QP_temp2,n,QP_temp3);
			//res_2 = scalar_max(scalar_abs(vec_min(QP_temp3,n)),scalar_abs(vec_max(QP_temp3,n)));
			//res = scalar_max(scalar_max(res_1,res_2),scalar_abs(mu));
			////printf("Res_1: %f\n",res_1);
			////printf("Res_2: %f\n",res_2);
			////printf("mu : %f\n",abs(mu));
			////printf("Residual: %f\n",res);

			////Terminal condition(KKT Residuel)
			//if(res<0.00001)
			//	break;

			////KKT residual (Wright) Termination Criteria
			//flag = TC_KKT_Wright(rd, rp, b_max,c_max, mu, 0.0001,m, n);
			//if(flag==1)
			//	break;


			 cons_obey_old = cons_obey;
			//// Himmelblau Termination Criteria
			mat_vec(G,x,n,n,QP_temp1);
			Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
			res_1 = scalar_abs(Fx_old-Fx);
			vec_sub(x_old,x,n,QP_temp1);
			vec_abs(QP_temp1,n,QP_temp2);
			res_2 = vec_max(QP_temp2,n);
			mat_vec(A,x,m,n,QP_temp1);
			vec_sub(b,QP_temp1,m,QP_temp2);
			cons_obey = vec_max(QP_temp2,m);
			res = scalar_max(res_1,res_2);
			res = scalar_max(res,cons_obey);
			//printf("The Fx_old is: %f \n",Fx_old);
			//printf("The Fx is: %f \n",Fx);
			// Update
			Fx_old = Fx;
			for(k_in=0;k_in<n;k_in++)
				x_old[k_in] = x[k_in];
			//printf("The res_1 is: %f \n",res_1);
			//printf("The res_2 is: %f \n",res_2);
			//printf("cons_obey is: %f \n",cons_obey);
			//fprintf(fp_Fx,"%f\t",Fx);
			//fprintf(fp_cons_obey,"%f\n",cons_obey);

			//Terminal condition(Himmelblau)
			//if(res<epsilon)
			//	break;

			////CDC
			CDC_flag = CDC(cons_obey, cons_obey_old, Fx, Fx_old, x, delta_x, G, c, mu, 0.5, 0.9);
			if(CDC_flag==1)
				{
				//	printf("CDC converged result\n");
					break;
				}
			else if(CDC_flag==2)
				{
					printf("Underconverged result\n");
					break;
				}

			//CDC_flag = CDC(cons_obey, cons_obey_old, Fx, Fx_old, x, delta_x, G, c, mu, 0.9, 0.9);
			//fprintf(fp_Fx,"%d\n",CDC_flag);


			//// Show new Fx value
			//mat_mul(x,1,n,G,n,QP_temp1);
			//Fx = 0.5*dot_product(temp1,x,n)+dot_product(c,x,n);
			//printf("This iteration F(x) value is: %f\n\n",Fx);


			end_time_section = Timestamp_get32();
			total_time_section = end_time_section - start_time_section - dead_time;
			Time_Section5 +=(double)total_time_section;

		}
		memcpy(delta_u,x,n*sizeof(double));
		QP_Iteration = Iter;
		//printf("QP_Iter:%d\n",QP_Iteration);

	}

	//printf("Total: Iteration: %d\n",Iter);


	QP_Iter_end_time = Timestamp_get32();
	Time_QP_Iter = (double)QP_Iter_end_time-QP_Iter_start_time-dead_time;

}
