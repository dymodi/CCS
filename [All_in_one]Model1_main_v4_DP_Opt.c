///*
// * main.c
// *
// */
//
//
///*
// *  ======== main.c ========
// */
//
//#include <stdio.h>
//#include <stdlib.h>
//#include <math.h>
//#include <limits.h>
////#include <time.h>
//
//#include <xdc/std.h>
//
//#include <xdc/runtime/Error.h>
//#include <xdc/runtime/System.h>
//#include <xdc/runtime/Types.h>
//#include <xdc/runtime/Memory.h>
//#include <xdc/runtime/Timestamp.h>
//
//#include <ti/sysbios/BIOS.h>
//
//#include <ti/sysbios/knl/Task.h>
//
////#include "DSPF_sp_mat_mul.h"
//
////#include "mpc_lib_DP.h"
////#include "PDIPM_V4.h"
////#include "IPM_V3_DP.h"
//
/////* Defines */
////#if defined(__TI_EABI__)
////#define kernel_size _kernel_size
////#endif
//
////extern char kernel_size;
//
////Here are some global variables.
//int P = 20;		//Prediction horizon
//int M = 3;		//Control horizon
//double Q = 1;	//
//double R = 0.1;
//
//int nu = 1;		//The number of inputs
//int ny = 1;		//The number of outputs
//int nx = 2;		//The number of original states.(Augment model has one more state)
//
//int ndec;		//Number of decision variables in QP solving.
//int mc;			//?Number of inequality constraints in QP solving.
//
//// Model parameters
//// double *A_e, *B_e, *C_e;
//double A_e[] = {1,0.1,0,0,0.9,0,1,0.1,1};
//double B_e[] = {0,0.0787,0};
//double C_e[] = {0,0,1};
//
////Constraint initialization
//double delta_U_p = 1;
//double delta_U_n = -1;
//double U_p = 2;
//double U_n = -2;
//double Y_p = 2;
//double Y_n = -2;
//
////VARIABLES DECLARARION
////Variable declaration for main
////Constant in online computing.
//double *L,*B,*G,*GL,*GL_T,*F,*Phi,*OMEGA_L,*N_OMEGA_L,*N_OMEGA_L_T;
////Constant in online computing.
//double *delta_u_ini,*y_ini,*lambda_ini;
//
////Variable declaration for mpc_DP
//double *aug_u_k_1;
//double *omega_r,*n_omega_r;
//double *c,*c_rev;
//double *h_x;
//double *delta_u_M; //delta_u in the M control horizon
//double *MPC_temp1,*MPC_temp2,*MPC_temp3,*MPC_temp4;
//double b_max,c_max,G_max,A_max,p_max;
//
////Variable declaration for IPM_v3_DP
//int MaxIteration = 20;
//double *x_old;
//double *rd, *rp, *rp_y;
//double *inv_y_lambda, *y_lambda_A;
//double *equation_F, *equation_b;
//double *delta_x, *delta_y, *delta_lambda;
//double *temp_cond;	//To check condition in the choice of alpha_tau_pri & alpha_tau_dual.
//double *x_try;
//double *sig_mu, *center_part;
//double *QP_temp1, *QP_temp2, *QP_temp3;
//
////Additional variable delcartion for priduip_v4
//double *delta_x_aff, *delta_y_aff, *delta_lambda_aff;
//double *sig_mu_lam_y, *aff_part;
//
////Variable declaration for feedback_v3
//double *fb_temp1, *fb_temp2;
//
////Variable declaration for mcf
//double *mcf_C;
//double *mcf_temp_vec;
//
////Variable declaration for alpha_decreas
//double *alpha_decreas_neg_delta_value, *alpha_decreas_neg_vec_value;
//double *alpha_decreas_temp1,*alpha_decreas_temp2;
//
////Variable declaration for line_solve
//double *line_solve_L,*line_solve_p;
//
////Time consuming record
//double Time_MPC = 0,Time_QP = 0,QP_Iteration = 0,Time_QP_non_Iter = 0,Time_QP_Iter = 0;
//double Time_Section1 = 0,Time_Section2 = 0,Time_Section3 = 0,Time_Section4 = 0,Time_Section5 = 0;
////Time consuming for subfunctions
//int mat_mul_times = 0;		double Time_Mat_mul = 0;
//int Mat_vec_times = 0;
//int vec_add_times = 0;		double Time_vec_add = 0;
//int vec_sub_times = 0;		double Time_vec_sub = 0;
//int vec_div_times = 0;		double Time_vec_div = 0;
//int vec_mul_times = 0;		double Time_vec_mul = 0;
//int bsxfun_times = 0; 		double Time_bsxfun = 0;
//int sca2vec_times = 0;		double Time_sca2vec = 0;
//int vec_rev_times = 0;		double Time_vec_rev = 0;
//int sca_vec_mut_times = 0;	double Time_sca_vec_mut = 0;
//int dot_product_times = 0;	double Time_dot_product = 0;
//int alpha_compute_times = 0;
//
//
////======================mpc_lib_DP.c==START======================================
//// Print the matrix to the consolo. Matrix stored in a array. w column, h rows.
//// Rosetta.
//void show_matrix(double *m, int w, int h)
//{
//	int i, j;
//	for (i = 0; i < h; i++) {
//		for (j = 0; j < w; j++)
//			printf("%2g ", m[i * w + j]);
//		putchar('\n');
//	}
//}
//
//
////Find the absloute value of a double parameter.
////2014.4.10 Yi Ding.
//double scalar_abs(double x)
//{
//	if(x<0)
//		return -x;
//	return x;
//}
//
//
////Find the minimun one of two scalars.
////2014.4.9 Yi Ding.
//double scalar_min(double x,double y)
//{
//    double z;
//    z=(x<y)?x:y;
//    return z;
//}
//
////Find the maximun one of two scalars.
////2014.4.11 Yi Ding.
//double scalar_max(double x,double y)
//{
//    double z;
//    z=(x>y)?x:y;
//    return z;
//}
//
////Filling a vector with a specific scalar.
////double *sca2vec(double sca, int n)
//void sca2vec(double sca, int n,double *vec)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	for(i=0;i<n;i++)
//		vec[i] = sca;
//
//	//end_time = Timestamp_get32();
//	//sca2vec_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_sca2vec+=(double)total_time;
//
//}
//
////Find the minimum one of a vector.
//////2014.4.9 Yi Ding.
//double vec_min(double *x, int n)
//{
//	double min;
//	int i;
//	min = x[0];
//	for(i=1;i<n;i++)
//	{
//		if(x[i]<min)
//			min = x[i];
//	}
//	return min;
//}
//
////Find the maximum one of a vector.
//////2014.4.10 Yi Ding.
//double vec_max(double *x, int n)
//{
//	double max;
//	int i;
//	max = x[0];
//	for(i=1;i<n;i++)
//	{
//		if(x[i]>max)
//			max = x[i];
//	}
//	return max;
//}
//
////Vector augment!
//void vec_aug(double *vec,int m, int n,double *aug)
//{
//	int i,j;
//	//m: vector element?
//	//n: how many vector?
//
//	for(i=0;i<n;i++)
//	{
//		for(j=0;j<m;j++)
//		{
//			aug[i*m+j] = vec[j];
//		}
//	}
//}
//
//
////Vector addation
////2014.4.9 Yi Ding.
////void vec_add(double *a, double *b, int n,double *c)
//void vec_add(double * restrict a, double * restrict b, int n,double *c)
//{
//	int i;
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	#pragma MUST_ITERATE (3);
//	for(i=0;i<n;i++)
//	{
//		c[i]=a[i]+b[i];
//	}
//
//	//end_time = Timestamp_get32();
//	//vec_add_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_vec_add+=(double)total_time;
//}
//
////Vector subtraction
////2014.4.9 Yi Ding.
////void vec_sub(double *a, double *b, int n,double *c)
//void vec_sub(double * restrict a, double * restrict b, int n,double * restrict c)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	#pragma MUST_ITERATE (3);
//	for(i=0;i<n;i++)
//	{
//		c[i]=a[i]-b[i];
//	}
//
//	//end_time = Timestamp_get32();
//	//vec_sub_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_vec_sub+=(double)total_time;
//
//	//return c;
//}
//
//// Return the product of a scalar and vector;
////void sca_vec_mut(double sca, double *vec, double n, double *pro)
//void sca_vec_mut(double sca, double * restrict vec, double n, double * restrict pro)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	//#pragma MUST_ITERATE (3);
//	for(i=0;i<n;i++)
//	{
//		pro[i] = sca*vec[i];
//	}
//
//	//end_time = Timestamp_get32();
//	//sca_vec_mut_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_sca_vec_mut+=(double)total_time;
//
//}
//
////This is a function to make a scalar[a] into a vecotr and substract another vector vec_b to get vec_x
////Yi Ding. 2014.5.21
//void sca_vec_sub(double a, double *vec_b, int length, double *vec_x)
//{
//	int i;
//
//	for(i=0;i<length;i++)
//		vec_x[i] = a - vec_b[i];
//
//}
//
////This is a function to do .* operation to 2 vectors.
////Not dot product, but the corresponding element to do multiplication and make a new vector.
////Yi Ding. 2014.5.21
////void vec_mul(double *a, double *b, int m, double *x)
//void vec_mul(double * restrict a, double * restrict b, int m, double * restrict x)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	#pragma MUST_ITERATE (3);
//	for(i=0;i<m;i++)
//		x[i] = a[i]*b[i];
//
//	//end_time = Timestamp_get32();
//	//vec_mul_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_vec_mul+=(double)total_time;
//
//}
//
////This is a function to do ./ operation to 2 vectors.
////Let the corresponding element to do division and make a new vector.
////Yi Ding. 2014.5.21
////void vec_div(double *a, double *b, int m, double *x)
//void vec_div(double * restrict a, double * restrict b, int m, double * restrict x)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	#pragma MUST_ITERATE (1);
//	for(i=0;i<m;i++)
//		x[i] = a[i]/b[i];
//
//	//end_time = Timestamp_get32();
//	//vec_div_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_vec_div+=(double)total_time;
//
//}
//
////This function implement the m-fun bsxfun() is Matlab, but only implement the mat-vec times function.
////Not the same usage with that in Matlab, be careful to use.
////Yi Ding. 2014.5.21
////void bsxfun(double *A, double *x, int row, int col, double *B)
//void bsxfun(double * restrict A, double * restrict x, int row, int col, double * restrict B)
//{
//	int i,j;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	#pragma MUST_ITERATE (12);
//	for(i=0;i<row;i++)
//		#pragma MUST_ITERATE (3);
//		for(j=0;j<col;j++)
//			B[i*col+j] = A[i*col+j]*x[i];
//
//	//end_time = Timestamp_get32();
//	//bsxfun_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_bsxfun += (double)total_time;
//}
//
////Calculate -a of a
////2014.4.9 Yi Ding.
////void vec_rev(double *a, int n, double *b)
//void vec_rev(double * restrict a, int n, double * restrict b)
//{
//	int i;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//start_time = Timestamp_get32();
//
//	#pragma MUST_ITERATE (1);
//	for(i=0;i<n;i++)
//	{
//		b[i] = -a[i];
//	}
//
//	//end_time = Timestamp_get32();
//	//vec_rev_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_vec_rev += (double)total_time;
//}
//
////Calculate |a| of a
////2014.5.14 Yi Ding.
////void vec_abs(double * a, int n, double *b)
//void vec_abs(double * restrict a, int n, double *restrict b)
//{
//	int i;
//	#pragma MUST_ITERATE (1);
//	for(i=0;i<n;i++)
//	{
//		if(a[i]<0)
//			b[i] = -a[i];
//		else
//			b[i] = a[i];
//	}
//}
//
////Make the vector into a diagonal form.
////2014.4.9 Yi Ding.
////Modified @2014.4.16
//void vec_diag(double *a, int n, double *diag)
//{
//	int i;
//	int j=0;
//
//	//Here we should clear the diag first! @2014.4.16
//	for(i=0;i<n*n;i++)
//		diag[i]=0;
//
//	//Insert the diagonal value.
//	for(i=0;i<n*n;i+=n)
//	{
//		diag[i+j]=a[j];
//		j++;
//	}
//}
//
//// Extract the diagonal elements of a diagonal matrix and form a vector.
//void diag2vec(double *diag, int n , double *vec)
//{
//	int i;
//	int j=0;
//	for(i=0;i<n*n;i+=n)
//	{
//		vec[j] = diag[i+j];
//		j++;
//	}
//}
//
////Calculate the inverse of a diagonal matrix.
////2014.4.9 Yi Ding.
//void inv_diag(double *diag, int n ,double *inv_diag)
//{
//	int i;
//	int j=0;
//	//Here's a trick to access diagnose element using a single cycle but not a double cycle.
//	for(i=0;i<n*n;i+=n)
//	{
//		inv_diag[i+j]=1/diag[i+j];
//		j++;
//	}
//}
//
//
////Add two matrices
////2014.4.9 Yi Ding.
//void mat_add(double *a, double *b, int row, int column, double *sum)
//{
//	int i;
//	for(i=0;i<row*column;i++)
//	{
//		sum[i] = a[i]+b[i];
//	}
//}
//
////Do the multiplication A*x = b;
////Yi Ding. 2014.4.25
////void mat_vec(double * A, double * x, int row, int col, double *b)
//void mat_vec(double * restrict A, double * restrict x, int row, int col, double *restrict b)
//{
//	int i,j;
//
//	#pragma MUST_ITERATE (3);
//	for(i=0;i<row;i++)
//	{
//		b[i]=0;
//		#pragma MUST_ITERATE (3);
//		for(j=0;j<col;j++)
//		{
//				b[i]+=A[i*col+j]*x[j];				//Attention!2014.5.5
//		}
//	}
//
//}
//
////Do the multiplication A*x = b;
////Yi Ding. 2014.4.25
////void mat_vec(double * A, double * x, int row, int col, double *b)
//void mat_vec_offline(double * restrict A, double * restrict x, int row, int col, double *restrict b)
//{
//	int i,j;
//
//	for(i=0;i<row;i++)
//	{
//		b[i]=0;
//		for(j=0;j<col;j++)
//		{
//				b[i]+=A[i*col+j]*x[j];				//Attention!2014.5.5
//		}
//	}
//
//}
//
//
////Do the multiplication x'*A = b;
////Yi Ding. 2014.7.23
//void vec_mat(double *A, double *x, int row, int col, double *b)
//{
//	int i,j;
//	#pragma MUST_ITERATE (3);
//	for(j=0;j<col;j++)
//	{
//		b[j]=0;
//		#pragma MUST_ITERATE (3);
//		for(i=0;i<row;i++)
//		{
//				b[j]+=A[i*col+j]*x[i];				//Attention!2014.5.5
//		}
//	}
//}
//
////Do the multiplication A*diag(x) = b;
////Yi Ding. 2014
//void mat_diag(double *A, double *diag, int row, int col, double *b)
//{
//	int i,j;
//
//	for(i=0;i<row;i++)
//	{
//		for(j=0;j<col;j++)
//			b[i*col+j] = A[i*col+j]*diag[j*col+j];
//	}
//}
//
////Used for LU Evaluate!
////2014.5.15
//void Lower_diag_vec(double *A, double *diag_vec, int m, double *b)
//{
//	int i,j;
//	double temp;
//	for(j=0;j<m;j++)
//	{
//		temp = diag_vec[j];
//		for(i=j;i<m;i++)
//		{
//			b[i*m+j] = A[i*m+j]*temp;
//			//show_matrix(b,m,m);
//		}
//	}
//}
//
////Do the A*X = B, here A, X, B are all diagnal matrix
////Yi Ding 2014.4.26
//void diag_diag(double *A, double *X, int n, double *B)
//{
//	int i;
//
//	for(i=0;i<n;i++)
//		B[i*n+i] = A[i*n+i]*X[i*n+i];
//}
//
////Do the A*x = b, here A is diag matrix, x is a vector
////Yi Ding 2014.4.26
//void diag_vec(double *A, double *x, int n, double *b)
//{
//	int i;
//	for(i=0;i<n;i++)
//		b[i] = A[i*n+i]*x[i];
//}
//
////Matrix multiplication. Just like DSP_sp_mat_mul.c in ti DSPlib.
////2014.4.8 Yi Ding.
//void mat_mul(double * restrict x1, int r1, int c1, double * restrict x2, int c2, double * restrict y)
//{
//	int i,j,k;
//	int length;
//
//	//Uint32 start_time,end_time,total_time,dead_time;
//	//start_time = Timestamp_get32();
//	//end_time = Timestamp_get32();
//	//dead_time = end_time - start_time;
//
//	//length = r1*c2;
//
//	//for(i=0;i<length;i++)
//	//{
//	//	y[i]=0;
//	//}
//
//	#pragma MUST_ITERATE (3);
//	for(i=0;i<r1;i++)
//	{
//	#pragma MUST_ITERATE (3);
//		for(j=0;j<c2;j++)
//		{
//			y[(i)*c2+(j)]=0;
//			#pragma MUST_ITERATE (12);
//			for(k=0;k<c1;k++)
//			{
//				y[(i)*c2+(j)]+=x1[(i)*c1+(k)]*x2[(k)*c2+(j)];
//			}
//		}
//	}
//	//end_time = Timestamp_get32();
//	//mat_mul_times++;
//	//total_time = end_time - start_time - dead_time;
//	//Time_Mat_mul+=(double)total_time;
//}
//
////Matrix multiplication. For off line use. For the accurate running time of on line mat_mul.
////2014.7.18 Yi Ding.
//void mat_mul_offline(double *x1, int r1, int c1, double *x2, int c2, double *y)
//{
//	int i,j,k;
//
//	for(i=0;i<r1;i++)
//	{
//		for(j=0;j<c2;j++)
//		{
//			y[(i)*c2+(j)]=0;
//			for(k=0;k<c1;k++)
//			{
//				y[(i)*c2+(j)]+=x1[(i)*c1+(k)]*x2[(k)*c2+(j)];
//			}
//		}
//	}
//}
//
//
////Compute the power of a matrix
////2014.4.16
//void mat_pow(double *mat, int row, int pow, double *MatPow)
//{
//	int i;
//	double *temp1;
//
//	temp1 = (double*)calloc(row*row,sizeof(double));
//	//temp2 = (double*)calloc(row*row,sizeof(double));
//
//	if(pow==0)
//	{
//		sca2vec(1,row,temp1);
//		vec_diag(temp1,row,MatPow);
//	}
//	else
//	{
//		memcpy(temp1,mat,row*row*sizeof(double));
//		memcpy(MatPow,mat,row*row*sizeof(double));
//		for(i=1;i<pow;i++)
//		{
//			mat_mul_offline(temp1,row,row,mat,row,MatPow);
//			memcpy(temp1,MatPow,row*row*sizeof(double));
//		}
//	}
//	free(temp1);
//}
//
//// Cholesky decomposition.
//// 2014.3.29 Yi Ding. Rosetta.
//void cholesky(double *A, int n, double *L)
//{
//	int i,j,k;
//	double s;
//	//double *L = (double*)calloc(n * n, sizeof(double));
// //   if (L == NULL)
// //       exit(EXIT_FAILURE);
//
//    for (i = 0; i < n; i++)
//        for (j = 0; j < (i+1); j++) {
//        	s = 0;
//            for (k = 0; k < j; k++)
//                s += L[i * n + k] * L[j * n + k];
//				L[i * n + j] = (i == j) ?
//            			sqrt(A[i * n + i] - s) :
//						//sqrtsp(A[i * n + i] - s) :
//                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
//        }
//}
//
////A Cholesky Factorization routine from the book Numerical Recipes in C
//// 2014.7.15 Yi Ding
////void chol_NRC(double *A, int n,double *p)
//void chol_NRC(double * restrict A, int n,double * restrict p)
//{
//	int i,j,k;
//	double sum;
//	#pragma MUST_ITERATE (3);
//	for (i=1;i<=n;i++)
//	{
//		#pragma MUST_ITERATE (3);
//		for (j=i;j<=n;j++)
//		{
//			for (sum=A[(i-1)*n+(j-1)],k=i-1;k>=1;k--)
//			{
//				sum -= A[(i-1)*n+(k-1)]*A[(j-1)*n+(k-1)];
//				//printf("sum: %f\n",sum);
//				//show_matrix(A,n,n);putchar('\n');
//				//printf("A[%d][%d]: %f\n",i,k,A[(i-1)*n+(k-1)]);
//				//printf("A[%d][%d]: %f\n",j,k,A[(j-1)*n+(k-1)]);
//			}
//			//show_matrix(A,n,n);putchar('\n');
//			//printf("A[%d][%d]: %f\n",i,j,A[(i-1)*n+(j-1)]);
//			if (i == j)
//			{
//				//printf("sum: %f\n",sum);
//				if (sum <= 0.0)
//					printf("Factorization Failed.\n");
//				p[i-1]=sqrt(sum);
//				//printf("p[%d]: %f\n",i,p[i-1]);
//			}
//			else A[(j-1)*n+(i-1)]=sum/p[i-1];
//			//show_matrix(A,n,n);putchar('\n');
//		}
//	}
//}
//
////Here, the diagnal D is stored in a vector.
////Modified Cholesky Factorization 2014.5.15
//void mcf(double *A, int m, double *L, double *D)
//{
//	int i,j,k;
//
//	double beta = 100;
//	double sigma = 0.1;
//	double temp;
//
//
//	for(i=0;i<m;i++)
//	{
//		L[i] = 0;
//		D[i] = 0;
//	}
//	for(i=0;i<m;i++)
//		L[i*m+i] = 1;
//
//	mcf_C[0] = A[0];
//	D[0] = scalar_max(scalar_max(scalar_abs(mcf_C[0]),(mcf_C[0]/beta)*(mcf_C[0]/beta)),sigma);
//	//printf("D[0] is: %f\n",D[0]);
//
//	for(k=1;k<m;k++)
//	{
//		L[k*m] = A[k*m]/D[0];
//	}
//	//show_matrix(L,m,m);putchar('\n');
//
//	for(j=1;j<m;j++)
//	{
//		temp = 0;
//		for(k=0;k<i;k++)
//			temp += D[k]*L[j*m+k]*L[j*m+k];
//		mcf_C[j*m+j] = A[j*m+j]-temp;
//		//printf("C is:\n");
//		//show_matrix(C,m,m);putchar('\n');
//
//		for(k=j-1;k<m;k++)
//			mcf_temp_vec[k-j+1] = mcf_C[k*m+j];
//		vec_abs(mcf_temp_vec,m-j+1,mcf_temp_vec);
//		//show_matrix(temp_vec,1,m-j+1);putchar('\n');
//		temp = vec_max(mcf_temp_vec,m-j+1);
//		temp = temp*temp;
//		D[j] = scalar_max(scalar_max(scalar_abs(mcf_C[j*m+j]),temp),sigma);
//		//printf("D is:\n");
//		//show_matrix(D,1,m);putchar('\n');
//
//		for(i=j+1;i<m;i++)
//		{
//			temp = 0;
//			for(k=0;k<j;k++)
//			{
//				temp += D[k]*L[i*m+k]*L[j*m+k];
//				//printf("L is:\n");
//				//show_matrix(L,m,m);
//				//printf("i is:%d\n",i);
//				//printf("j is:%d\n",j);
//				//printf("k is:%d\n",k);
//				//printf("D[k] is:%f\n",D[k]);
//				//printf("L[i*m+k] is:%f\n",L[i*m+k]);
//				//printf("L[j*m+k] is:%f\n",L[j*m+k]);
//				//printf("temp is:%f\n",temp);
//				//putchar('\n');putchar('\n');
//			}
//
//			mcf_C[i*m+j] = A[i*m+j]-temp;
//			//printf("C is:\n");
//			//show_matrix(C,m,m);putchar('\n');
//
//			L[i*m+j] = mcf_C[i*m+j]/D[j];
//			//printf("L is:\n");
//			//show_matrix(L,m,m);putchar('\n');
//		}
//	}
//}
//
//
//// Do the dot product for two array. n is the number of elements
//// 2014.3.29 Yi Ding. Rosetta.
////double dot_product(double * a, double * b, int n)
//double dot_product(double * restrict a, double * restrict b, int n)
//{
//        double sum = 0;
//        int i;
//
//    	//Uint32 start_time,end_time,total_time,dead_time;
//        //start_time = Timestamp_get32();
//        //end_time = Timestamp_get32();
//        //dead_time = end_time - start_time;
//
//    	//start_time = Timestamp_get32();
//
//		#pragma MUST_ITERATE (3);
//        for (i = 0; i < n; i++) {
//                sum += a[i] * b[i];
//        }
//
//        //end_time = Timestamp_get32();
//        //dot_product_times++;
//        //total_time = end_time - start_time - dead_time;
//        //Time_dot_product += (double)total_time;
//
//        return sum;
//
//}
//
////This function implement a special searching method for alpha.
////The idea comes from the quad_wright() and symbol meanings see algorithms.
//// Yi Ding. 2014.5.21
//double alpha_decreas(double *vec, double *delta_vec,int length, double tau)
//{
//	int i;
//	int neg_num = 0;
//	//double *neg_index;
//	double alpha = 0.9999995;
//
//	//neg_index = (double*)calloc(length,sizeof(double));
//
//	#pragma MUST_ITERATE (12);
//	for(i=0;i<length;i++)
//	{
//		if(delta_vec[i]<0)
//		{
//			//neg_index[neg_num] = i;
//			alpha_decreas_neg_delta_value[neg_num] = delta_vec[i];
//			alpha_decreas_neg_vec_value[neg_num] = vec[i];
//			neg_num = neg_num + 1;
//		}
//	}
//
//	if(neg_num>0)
//	{
//		if(tau==1)
//			vec_rev(alpha_decreas_neg_vec_value,neg_num,alpha_decreas_temp1);
//		else
//			sca_vec_mut(-tau,alpha_decreas_neg_vec_value,neg_num,alpha_decreas_temp1);
//		vec_div(alpha_decreas_temp1,alpha_decreas_neg_delta_value,neg_num,alpha_decreas_temp2);
//		alpha = 0.9995*vec_min(alpha_decreas_temp2,neg_num);
//	}
//
//	return alpha;
//
//}
//
//// Alpha calculation method (while loop)
//double alpha_decreas_while(double *y, double *lambda, double *delta_y, double *delta_lambda, int m)
//{
//
//	double alpha_cond;	//Use to check whether the alpha is valid.
//	double alpha = 1;
//	double *temp1,*temp2,*temp3;
//
//	temp1 = (double*)calloc(m,sizeof(double));
//	temp2 = (double*)calloc(m,sizeof(double));
//	temp3 = (double*)calloc(m,sizeof(double));
//
//	sca_vec_mut(alpha,delta_y,m,temp1);
//	vec_add(y,temp1,m,temp2);
//	sca_vec_mut(alpha,delta_lambda,m,temp1);
//	vec_add(lambda,temp1,m,temp3);
//	alpha_cond = scalar_min(vec_min(temp2,m),vec_min(temp3,m));
//	while(alpha_cond<0)
//	{
//		alpha = alpha - 0.01;		//Attention!
//		//alpha_aff = alpha_aff/2;
//		if(alpha<=0)
//			break;
//		sca_vec_mut(alpha,delta_y,m,temp1);
//		vec_add(y,temp1,m,temp2);
//		sca_vec_mut(alpha,delta_lambda,m,temp1);
//		vec_add(lambda,temp1,m,temp3);
//		alpha_cond = scalar_min(vec_min(temp2,m),vec_min(temp3,m));
//	}
//
//	free(temp1);free(temp2);free(temp3);
//
//	return alpha;
//}
//// Test accomplished! 2014.4.9
//
//
//// Do back substatution to solve linear equation using LU.
//// 2013.3.30 Yi Ding.
////Modified @2014.4.18
//void luEvaluate(double *L,double *U, double*b,int n,double *x)
//{
//	double *y = (double*)calloc(n,sizeof(double));
//	int i,j;
//	double temp = 0;
//	if(x == NULL || y == NULL)
//		exit(0);
//
//	//Foward solve Ly = b;
//	y[0] = b[0]/L[0];
//
//	#pragma MUST_ITERATE (3);
//	for(i=1;i<n;i++)
//	{
//		for(j=0;j<i;j++)
//		{
//			temp += L[i*n+j]*y[j];
//		}
//		y[i] = b[i] - temp;
//		y[i] = y[i]/L[i*n+i];
//		temp = 0;
//	}
//	//show_matrix(y,1,n);
//
//	//Backward solve Ux = y
//	x[n-1] = y[n-1]/U[n*n-1];
//	temp = 0;
//	#pragma MUST_ITERATE (1);
//	for(i=n-2;i>=0;i--)
//	{
//		for(j=i+1;j<n;j++)
//		{
//			temp += U[i*n+j]*x[j];
//		}
//		x[i] = y[i] - temp;
//		x[i] = x[i]/U[i*n+i];
//		temp = 0;
//	}
//	free(y);
//}
//
//
//// Do back substatution to solve linear equation using LU.
//// A routine from the book Numerical Recipes in C
//// Yi Ding 2014.7.15
////void luEvaluate_NRC(double *L, double*b,double *p,int n,double *x)
//void luEvaluate_NRC(double * restrict L, double* restrict b,double * restrict p,int n,double * restrict x)
//{
//int i,k;
//double sum;
//
//	#pragma MUST_ITERATE (3);
//	for (i=1;i<=n;i++)
//	{
//		//#pragma MUST_ITERATE (3);
//		for (sum=b[i-1],k=i-1;k>=1;k--)
//			sum -= L[(i-1)*n+(k-1)]*x[k-1];
//		x[i-1]=sum/p[i-1];
//	}
//
//	#pragma MUST_ITERATE (3);
//	for (i=n;i>=1;i--)
//	{
//		//#pragma MUST_ITERATE (3);
//		for (sum=x[i-1],k=i+1;k<=n;k++)
//			sum -= L[(k-1)*n+(i-1)]*x[k-1];
//		x[i-1]=sum/p[i-1];
//	}
//}
//
////Transpose the matrix m and store the matrix in w.
////2014.4.8 Yi Ding.
//void transpose(double *m, int w, int h)
//{
//	int start, next, i;
//	double tmp;
//
//	for (start = 0; start <= w * h - 1; start++) {
//		next = start;
//		i = 0;
//		do {	i++;
//			next = (next % h) * w + next / h;
//		} while (next > start);
//		if (next < start || i == 1) continue;
//
//		tmp = m[next = start];
//		do {
//			i = (next % h) * w + next / h;
//			m[next] = (i == start) ? tmp : m[i];
//			next = i;
//		} while (next > start);
//	}
//}
//
////Solve linear equation using cholesky decomposition.
////Simply calling two other functions: cf() and luEvaluate()
////2014.3.30 Yi Ding.
////double *line_solve(double *A, double *b, int n)
//void line_solve(double *A, double *b, int n, double *x)
//{
//	//double *U = (double*)calloc(n * n, sizeof(double));
//	//double *D = (double*)calloc(n, sizeof(double));
//	//double *LD = (double*)calloc(n * n, sizeof(double));
//
//
//	////Below is Standard Cholesky Factorization
//	//cholesky(A,n,L);
//	//memcpy(U,L,n*n*sizeof(double));
//	//transpose(U, n, n);
//	//luEvaluate(L,U,b,n,x);
//
//	//Below is the Linear solve process from the book Numerical Recipes in C
//	memcpy(line_solve_L,A,n*n*sizeof(double));
//	chol_NRC(line_solve_L,n,line_solve_p);
//	luEvaluate_NRC(line_solve_L,b,line_solve_p,n,x);
//
//	////Below is Modified Cholesky Factorization
//	//mcf(A,n,L,D);
//	//memcpy(U,L,n*n*sizeof(double));
//	//transpose(U, n, n);
//	//Lower_diag_vec(L,D,n,LD);
//	//luEvaluate(LD,U,b,n,x);
//
//
//	//free(U);
//	//free(D);
//	//free(LD);
//
//}
//
//
////The feedback step in MPC, a transformation of funciton fankui.m
//void feedback_v3(double *x_k_1, double *y_k, double *u_k_1, double *u_k_2, double *x_k,double *L)
//{
//	int i;
//
//	//vec_sub(u_k_1,u_k_2,nu,fb_temp1);
//	for(i=0;i<nu;i++)
//	{
//		fb_temp1[i] = u_k_1[i] - u_k_2[i];
//	}
//	mat_vec_offline(B_e,fb_temp1,nx+ny,nu,fb_temp2);
//	mat_vec_offline(A_e,x_k_1,nx+ny,nx+ny,fb_temp1);
//	vec_add(fb_temp1,fb_temp2,nx+ny,x_k);
//
//	mat_vec_offline(C_e,x_k,ny,nx+ny,fb_temp1);
//	//vec_sub(y_k,fb_temp1,ny,fb_temp2);
//	for(i=0;i<ny;i++)
//	{
//		fb_temp2[i] = y_k[i] - fb_temp1[i];
//	}
//	mat_vec_offline(L,fb_temp2,nx+ny,ny,fb_temp1);
//	vec_add(x_k,fb_temp1,nx+ny,fb_temp2);
//
//	memcpy(x_k,fb_temp2,(nx+ny)*sizeof(double));
//
//}
////Feedbcak test accomplished @2014.4.14
//
////The formation of omega_r.
////This version has y contraints.
//void omega_r_form(double *aug_u_k_1, double *x_k, double *omega_r, double *F)
//{
//	int i;
//	double *temp1;
//	temp1 = (double*)calloc(ny*P,sizeof(double));
//	mat_vec_offline(F,x_k,ny*P,nx+ny,temp1);
//
//	for(i=0;i<nu*M;i++)
//		omega_r[i] = delta_U_p;
//	for(i=nu*M;i<2*nu*M;i++)
//		omega_r[i] = -delta_U_n;
//	for(i=2*nu*M;i<3*nu*M;i++)
//		omega_r[i] = U_p - aug_u_k_1[i-2*nu*M];
//	for(i=3*nu*M;i<4*nu*M;i++)
//		omega_r[i] = -U_n + aug_u_k_1[i-3*nu*M];
//	for(i=4*nu*M;i<4*nu*M+ny*P;i++)
//		omega_r[i] = Y_p - temp1[i-4*nu*M];
//	for(i=4*nu*M+ny*P;i<4*nu*M+2*ny*P;i++)
//		omega_r[i] = -Y_n + temp1[i-4*nu*M-ny*P];
//	free(temp1);
//}
//
////The formation of omega_r.
////This version does not have y contraints.
//void omega_r_form2(double *aug_u_k_1, double *x_k, double *omega_r, double *F)
//{
//	int i;
//	//double *temp1;
//	//temp1 = (double*)calloc(ny*P,sizeof(double));
//	////mat_mul(F,ny*P,nx+ny,x_k,1,temp1);
//	//mat_vec(F,x_k,ny*P,nx+ny,temp1);
//
//	for(i=0;i<nu*M;i++)
//		omega_r[i] = delta_U_p;
//	for(i=nu*M;i<2*nu*M;i++)
//		omega_r[i] = -delta_U_n;
//	for(i=2*nu*M;i<3*nu*M;i++)
//		omega_r[i] = U_p - aug_u_k_1[i-2*nu*M];
//	for(i=3*nu*M;i<4*nu*M;i++)
//		omega_r[i] = -U_n + aug_u_k_1[i-3*nu*M];
//	//for(i=4*nu*M;i<4*nu*M+ny*P;i++)
//	//	omega_r[i] = Y_p - temp1[i-4*nu*M];
//	//for(i=4*nu*M+ny*P;i<4*nu*M+2*ny*P;i++)
//	//	omega_r[i] = -Y_n + temp1[i-4*nu*M-ny*P];
//	//free(temp1);
//}
//
//
////Find the parameter F and Phi, not in the MPC loop. Relaxation on performance.
////Yi Ding. 2014.4.16
//void fphi(double *F, double *Phi)
//{
//	int i,j,k,kk;
//	double *temp1, *temp2;
//	temp1 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));
//	temp2 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));
//
//	mat_mul_offline(C_e,ny,nx+ny,A_e,nx+ny,temp2);
//	for(i=0;i<P;i++)
//	{
//		for(j=0;j<ny*(ny+nx);j++)
//		{
//			F[i*(ny*(ny+nx))+j] = temp2[j];
//		}
//		memcpy(temp1,temp2,ny*(nx+ny)*sizeof(double));
//		mat_mul_offline(temp1,ny,nx+ny,A_e,nx+ny,temp2);
//	}
//
//	for(i=0;i<M;i++)
//	{
//		for(j=0;j<P;j++)
//		{
//			if(j<i)
//			{
//				sca2vec(0,nu*ny,temp1);
//			}
//			else
//			{
//				mat_pow(A_e,nx+ny,j-i,temp1);
//				mat_mul_offline(C_e,ny,nx+ny,temp1,nx+ny,temp2);
//				mat_mul_offline(temp2,ny,nx+ny,B_e,nu,temp1);
//			}
//			for(k=0;k<nu;k++)
//			{
//				for(kk=0;kk<ny;kk++)
//				{
//					Phi[j*ny*nu*M+kk*nu*M+i*nu+k] = temp1[kk*nu+k];
//				}
//			}
//		}
//	}
//	free(temp1);free(temp2);
//}
////fphi() Test accomplished @2014.4.16
////Attention! Find a bug in vec_diag(), modified!
//
////Form the B used in the OMEGA_L
//void form_B(double *B)
//{
//	int i,j,k;
//	for(i=0;i<nu*M;i++)
//	{
//		for(j=0;j<nu*M;j++)
//		{
//			for(k=0;k<M;k++)
//				if(i==(j+k*nu))
//					B[i*nu*M+j]=1;
//		}
//	}
//}
//
////Form the OMEGA_L for the constrains in the optimization problems.
////This version has y constraints.
//void OMEGA_L_form(double *B, double *Phi, double *OMEGA_L)
//{
//	int i,j;
//
//	for(i=0;i<nu*M;i++)
//	{
//		for(j=0;j<nu*M;j++)
//		{
//			if(i==j)
//				OMEGA_L[i*nu*M+j] = 1;
//		}
//	}
//	for(i=nu*M;i<2*nu*M;i++)
//	{
//		for(j=0;j<nu*M;j++)
//		{
//			if((i-nu*M)==j)
//				OMEGA_L[i*nu*M+j] = -1;
//		}
//	}
//	for(i=(2*nu*M)*nu*M;i<(3*nu*M)*nu*M;i++)
//		OMEGA_L[i] = B[i-(2*nu*M)*nu*M];
//	for(i=(3*nu*M)*nu*M;i<(4*nu*M)*nu*M;i++)
//		OMEGA_L[i] = -B[i-(3*nu*M)*nu*M];
//	for(i=(4*nu*M)*nu*M;i<(4*nu*M+ny*P)*nu*M;i++)
//		OMEGA_L[i] = Phi[i-(4*nu*M)*nu*M];
//	for(i=(4*nu*M+ny*P)*nu*M;i<(4*nu*M+2*ny*P)*nu*M;i++)
//		OMEGA_L[i] = -Phi[i-(4*nu*M+ny*P)*nu*M];
//}
//
////Form the OMEGA_L for the constrains in the optimization problems.
////This version does not have y constraints.
//void OMEGA_L_form2(double *B, double *Phi, double *OMEGA_L)
//{
//	int i,j;
//
//	for(i=0;i<nu*M;i++)
//	{
//		for(j=0;j<nu*M;j++)
//		{
//			if(i==j)
//				OMEGA_L[i*nu*M+j] = 1;
//		}
//	}
//	for(i=nu*M;i<2*nu*M;i++)
//	{
//		for(j=0;j<nu*M;j++)
//		{
//			if((i-nu*M)==j)
//				OMEGA_L[i*nu*M+j] = -1;
//		}
//	}
//	for(i=(2*nu*M)*nu*M;i<(3*nu*M)*nu*M;i++)
//		OMEGA_L[i] = B[i-(2*nu*M)*nu*M];
//	for(i=(3*nu*M)*nu*M;i<(4*nu*M)*nu*M;i++)
//		OMEGA_L[i] = -B[i-(3*nu*M)*nu*M];
//	//for(i=(4*nu*M)*nu*M;i<(4*nu*M+ny*P)*nu*M;i++)
//	//	OMEGA_L[i] = Phi[i-(4*nu*M)*nu*M];
//	//for(i=(4*nu*M+ny*P)*nu*M;i<(4*nu*M+2*ny*P)*nu*M;i++)
//	//	OMEGA_L[i] = -Phi[i-(4*nu*M+ny*P)*nu*M];
//}
//
//
//
////Form the G for the opimization problem.
//void form_G(double *Phi,double *G)
//{
//	double *temp1,*temp2,*temp3;
//	int max_mem;
//
//	max_mem = scalar_max(P*ny*P*ny,M*nu*P*ny);
//	temp1 = (double*)calloc(max_mem,sizeof(double));
//	temp2 = (double*)calloc(max_mem,sizeof(double));
//	temp3 = (double*)calloc(max_mem,sizeof(double));
//	sca2vec(Q,ny*P,temp1);
//
//	vec_diag(temp1,ny*P,temp2);
//	memcpy(temp1,Phi,M*nu*P*ny*sizeof(double));
//	transpose(temp1, M*nu, P*ny);
//
//	mat_mul_offline(temp1,M*nu,P*ny,temp2,P*ny,temp3);
//	mat_mul_offline(temp3,M*nu,P*ny,Phi,M*nu,temp1);
//	sca2vec(R,nu*M,temp2);
//	vec_diag(temp2,nu*M,temp3);
//	mat_add(temp1,temp3,nu*M,nu*M,G);
//
//	free(temp1);free(temp2);free(temp3);
//}
//
////This is a routine to mimic the model using the disturbed state-space model.
//void process_model(double* xm, double *u_k, double *y_k)
//{
//	double Ac[] = {1,0.1,0,0.9};
//	double Bc[] = {0,0.0787};
//	double Cc[] = {1,0};
//
//	int i;
//
//	double *temp1,*temp2;
//	temp1 = (double*)calloc(nx,sizeof(double));
//	temp2 = (double*)calloc(nx,sizeof(double));
//
//	mat_vec_offline(Ac,xm,nx,nx,temp1);
//	mat_vec_offline(Bc,u_k,nx,nu,temp2);
//	//vec_add(temp1,temp2,nx,xm);
//	for(i=0;i<nx;i++)
//	{
//		xm[i]=temp1[i]+temp2[i];
//	}
//
//	mat_vec_offline(Cc,xm,ny,nx,y_k);
//
//	free(temp1);free(temp2);
//}
//
//void SP_wright(double *delta_u_ini, double *y_bar, double *lambda_bar, double *G, double *c, double *A,double *A_t, double *b,double *y_ini, double *lambda_ini)
//{
//
//	int i;
//
//	double *delta_x_aff = (double*)calloc(nu*M,sizeof(double));
//	double *delta_y_aff = (double*)calloc(4*nu*M,sizeof(double));
//	double *delta_lambda_aff = (double*)calloc(4*nu*M,sizeof(double));
//
//	mat_vec(G,delta_u_ini,nu*M,nu*M,MPC_temp1);
//	mat_vec(A_t,lambda_bar,nu*M,4*nu*M,MPC_temp2);
//	vec_sub(MPC_temp1, MPC_temp2,nu*M,MPC_temp3);
//	vec_add(MPC_temp3,c,nu*M,rd);
//	//printf("rd is:\n");show_matrix(rd,1,n);putchar('\n');
//
//	mat_vec(A,delta_u_ini,4*nu*M,nu*M,MPC_temp1);
//	vec_sub(b,MPC_temp1,4*nu*M,rp_y);
//	vec_sub(MPC_temp1,y_bar,4*nu*M,MPC_temp2);
//	vec_sub(MPC_temp2,b,4*nu*M,rp);
//	//printf("rp is:\n");show_matrix(rp,1,m);putchar('\n');
//
//	vec_div(lambda_bar,y_bar,4*nu*M,inv_y_lambda);
//	bsxfun(A,inv_y_lambda,4*nu*M,nu*M,y_lambda_A);
//	//printf("inv_y_lambda is :\n");show_matrix(inv_y_lambda,1,m);putchar('\n');
//	//printf("y_lambda_A is :\n");show_matrix(y_lambda_A,n,m);putchar('\n');
//
//	mat_mul(A_t,nu*M,4*nu*M,y_lambda_A,nu*M,MPC_temp1);
//	vec_add(G,MPC_temp1,nu*M*nu*M,equation_F);
//	//printf("equation_F is: \n");show_matrix(equation_F,n,n);putchar('\n');
//
//	vec_mul(inv_y_lambda,rp_y,4*nu*M,MPC_temp1);
//	mat_vec(A_t,MPC_temp1,nu*M,4*nu*M,MPC_temp2);
//	vec_sub(MPC_temp2,rd,nu*M,equation_b);
//	//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');
//
//	line_solve(equation_F,equation_b,nu*M,delta_x_aff);
//
//	//printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');
//
//	mat_vec(A,delta_x_aff,4*nu*M,nu*M,MPC_temp1);
//	vec_add(MPC_temp1,rp,4*nu*M,delta_y_aff);
//	//printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');
//
//	vec_mul(inv_y_lambda,delta_y_aff,4*nu*M,MPC_temp2);
//	vec_rev(MPC_temp2,4*nu*M,MPC_temp3);
//	vec_sub(MPC_temp3,lambda_bar,4*nu*M,delta_lambda_aff);
//	//printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');
//
//	vec_add(y_bar,delta_y_aff,4*nu*M,MPC_temp1);
//	vec_abs(MPC_temp1,4*nu*M,y_ini);
//	for(i=0;i<4*nu*M;i++)
//	{
//		if(y_ini[i]<1)
//			y_ini[i]=1;
//	}
//
//	vec_add(lambda_bar,delta_lambda_aff,4*nu*M,MPC_temp1);
//	vec_abs(MPC_temp1,4*nu*M,lambda_ini);
//	for(i=0;i<4*nu*M;i++)
//	{
//		if(lambda_ini[i]<1)
//			lambda_ini[i]=1;
//	}
//
//	//printf("lambda_ini is: \n");show_matrix(lambda_ini,1,4*nu*M);putchar('\n');
//	//printf("y_ini is: \n");show_matrix(y_ini,1,4*nu*M);putchar('\n');
//	//printf("delta_u_bar is: \n");show_matrix(delta_u_ini,1,nu*M);putchar('\n');
//
//
//	free(delta_x_aff);free(delta_y_aff);free(delta_lambda_aff);
//
//}
//
//// KKT residual (Wright) Termination Criteria
//double TC_KKT_Wright(double *rd, double *rp, double b_max,double c_max, double mu, double epsilon,int mc, int ndec)
//{
//	int flag = 0;
//	double res_1, res_2, res;
//	//printf("rp:\n");show_matrix(rp,1,mc);putchar('\n');
//	res_1 = scalar_max(scalar_abs(vec_min(rp,mc)),scalar_abs(vec_max(rp,mc)));
//	//printf("res_1: %f\n",res_1);
//	res_1 = res_1/b_max;
//	//printf("res_1: %f\n\n",res_1);
//	//printf("rd:\n");show_matrix(rd,1,ndec);
//	res_2 = scalar_max(scalar_abs(vec_min(rd,ndec)),scalar_abs(vec_max(rd,ndec)));
//	//printf("res_2: %f\n",res_2);
//	res_2 = res_2/c_max;
//	//printf("res_2: %f\n\n",res_2);
//	res = scalar_max(scalar_max(res_1,res_2),scalar_abs(mu));
//	//printf("res: %f\n\n",res);
//	if(res<epsilon)
//		flag = 1;
//	return flag;
//	//printf("Res_1: %f\n",res_1);
//	//printf("Res_2: %f\n",res_2);
//	//printf("mu : %f\n",abs(mu));
//	//printf("Residual: %f\n",res);
//}
//
////============================mpc_lib_DP.c==END=============================================
//
//
//
///*
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
// */
//void IPM_V3_SP(double *G, double *GL,double *GL_T, double *c, double *A, double *A_t, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u)
//{
//
//	double Fx_old, Fx, cons_obey;
//	int Iter = 0;
//	double mu,mu_old;
//	//double alpha_aff_pri = 0.999995;
//	//double alpha_aff_dual = 0.999995;
//	//double alpha_aff = 0.99995;
//	double sigma = 0.1;
//	double res = 0.3;
//	//double tau = 0.7;
//	double alpha_tau_pri = 0.999995;
//	double alpha_tau_dual = 0.999995;
//	double alpha;
//	double res_1,res_2;
//
//	int flag;
//	int k, k_in;
//	int i,j;
//
//	UInt32 start_time,end_time,dead_time,total_time;
//	UInt32 start_time_section,end_time_section,total_time_section;
//	UInt32 QP_non_Iter_start_time,QP_non_Iter_end_time;
//	UInt32 QP_Iter_start_time,QP_Iter_end_time;
//
//	start_time = Timestamp_get32();
//	end_time = Timestamp_get32();
//	dead_time = end_time - start_time;
//
//	//Make A transpose as a variable
//
//	////Show the problem paremeters.
//	//show_matrix(G,n,n);putchar('\n');
//	//show_matrix(c,1,n);putchar('\n');
//	//show_matrix(A,n,m);putchar('\n');
//	//show_matrix(b,1,m);putchar('\n');
//	//show_matrix(x,1,n);putchar('\n');
//	//show_matrix(y,1,m);putchar('\n');
//	//show_matrix(lambda,1,m);putchar('\n');
//
//	QP_non_Iter_start_time = Timestamp_get32();
//
//
//	//Initial Fx_old value
//	//mat_vec(G,x,n,n,QP_temp1);
//	//#pragma MUST_ITERATE (3);
//	for(i=0;i<n;i++)
//	{
//		QP_temp1[i]=0;
//		//#pragma MUST_ITERATE (3);
//		for(j=0;j<n;j++)
//		{
//			QP_temp1[i]+=G[i*n+j]*x[j];				//Attention!2014.5.5
//		}
//	}
//
//
//	Fx_old = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
//	//printf("This initial F(x) value is: %f\n",Fx_old);
//
//	mu = dot_product(y,lambda,m)/m;
//
//	// Scale prepararion
//	//vec_abs(b,m,QP_temp1);
//	//b_max = vec_max(QP_temp1,m);
//	//vec_abs(c,n,QP_temp1);
//	//c_max = vec_max(QP_temp1,n);
//
//	//Initial x_old value
//	//#pragma MUST_ITERATE (3);
//	for(k=0;k<n;k++)
//		x_old[k] = x[k];
//
//	//Check if the Uncontrained solution meet the constraints
//
//	//printf("GL is:\n");show_matrix(GL,n,n);putchar('\n');
//	//printf("GL_T is:\n");show_matrix(GL_T,n,n);putchar('\n');
//	//printf("c_rev is:\n");show_matrix(c_rev,1,n);putchar('\n');
//	luEvaluate(GL,GL_T,c_rev,n,x_try);
//	//printf("x_try is:\n");show_matrix(x_try,1,n);putchar('\n');
//	mat_vec(A,x_try,m,n,QP_temp1);
//	//printf("QP_temp1 is:\n");show_matrix(QP_temp1,1,m);putchar('\n');
//	vec_sub(QP_temp1,b,m,QP_temp2);
//	//printf("QP_temp2 is:\n");show_matrix(QP_temp2,1,m);putchar('\n');
//
//	QP_non_Iter_end_time = Timestamp_get32();
//	Time_QP_non_Iter = (double)QP_non_Iter_end_time-QP_non_Iter_start_time-dead_time;
//
//	QP_Iter_start_time = Timestamp_get32();
//
//	//printf("(vec_min(QP_temp2,m):%f\n",vec_min(QP_temp2,m));
//	if((vec_min(QP_temp2,m))>0)
//	{
//		memcpy(delta_u,x_try,n*sizeof(double));
//		QP_Iteration = Iter;
//	}
//	else
//	{
//		// If not global optimal, do iterative process.
//		//printf("MaxIteration is : %d \n",MaxIteration);
//		for(k=1;k<=MaxIteration;k++)
//		{
//			//start_time = Timestamp_get32();
//			start_time_section = Timestamp_get32();
//
//			mat_vec(G,x,n,n,QP_temp1);
//			mat_vec(A_t,lambda,n,m,QP_temp2);
//			vec_sub(QP_temp1, QP_temp2, n,QP_temp3);
//			vec_add(QP_temp3,c,n,rd);
//			//printf("rd is:\n");show_matrix(rd,1,n);putchar('\n');
//
//			mat_vec(A,x,m,n,QP_temp1);
//			vec_sub(b,QP_temp1,m,rp_y);
//			vec_sub(QP_temp1,y,m,QP_temp2);
//			vec_sub(QP_temp2,b,m,rp);
//			//printf("rp is:\n");show_matrix(rp,1,m);putchar('\n');
//
//			vec_div(lambda,y,m,inv_y_lambda);
//			bsxfun(A,inv_y_lambda,m,n,y_lambda_A);
//			//printf("inv_y_lambda is :\n");show_matrix(inv_y_lambda,1,m);putchar('\n');
//			//printf("y_lambda_A is :\n");show_matrix(y_lambda_A,n,m);putchar('\n');
//
//			//printf("G is: \n");show_matrix(G,n,n);putchar('\n');
//			//printf("A_t is: \n");show_matrix(A_t,m,n);putchar('\n');
//			mat_mul(A_t,n,m,y_lambda_A,n,QP_temp1);	//Time Consuming! Much better if custom!
//			vec_add(G,QP_temp1,n*n,equation_F);
//			//printf("equation_F is: \n");show_matrix(equation_F,n,n);putchar('\n');
//
//			sca2vec(sigma*mu,m,sig_mu);
//
//			vec_div(sig_mu,lambda,m,QP_temp2);
//			vec_add(rp_y,QP_temp2,m,center_part);
//			vec_mul(inv_y_lambda,center_part,m,QP_temp1);
//			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
//			vec_sub(QP_temp2,rd,n,equation_b);
//			//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');
//
//			end_time_section = Timestamp_get32();
//			total_time_section = end_time_section - start_time_section - dead_time;
//			Time_Section1 += (double)total_time_section;
//
//			start_time_section = Timestamp_get32();
//
//			//line_solve(equation_F,equation_b,n,delta_x);
//			memcpy(line_solve_L,equation_F,n*n*sizeof(double));
//			chol_NRC(line_solve_L,n,line_solve_p);
//			luEvaluate_NRC(line_solve_L,equation_b,line_solve_p,n,delta_x);
//
//			//printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');
//
//			mat_vec(A,delta_x,m,n,QP_temp1);
//			vec_add(QP_temp1,rp,m,delta_y);
//			//printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');
//
//			vec_mul(lambda,delta_y,m,QP_temp2);
//			vec_sub(QP_temp2,sig_mu,m,QP_temp1);
//			vec_div(QP_temp1,y,m,QP_temp2);
//			vec_rev(QP_temp2,m,QP_temp3);
//			vec_sub(QP_temp3,lambda,m,delta_lambda);
//			//printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');
//
//			end_time_section = Timestamp_get32();
//			total_time_section = end_time_section - start_time_section - dead_time;
//			Time_Section2 +=(double)total_time_section;
//
//
//			start_time_section = Timestamp_get32();
//
//			//alpha_compute_times++;
//
//			alpha_tau_pri = alpha_decreas(y,delta_y,m,1);
//			alpha_tau_dual = alpha_decreas(lambda,delta_lambda,m,1);
//			alpha = scalar_min(alpha_tau_pri,alpha_tau_dual);
//			alpha = scalar_min(0.999999,alpha);
//			//printf("alpha is : %f\n",alpha);
//
//			//alpha = alpha_decreas_while(y, lambda, delta_y, delta_lambda, m);
//
//			end_time_section = Timestamp_get32();
//			total_time_section = end_time_section - start_time_section - dead_time;
//			Time_Section3 +=(double)total_time_section;
//
//			start_time_section = Timestamp_get32();
//
//			mu_old = mu;
//
//			//Obtain new iteration value.
//			sca_vec_mut(alpha,delta_x,n,QP_temp1);
//			vec_add(x,QP_temp1,n,x);
//			sca_vec_mut(alpha,delta_y,m,QP_temp1);
//			vec_add(y,QP_temp1,m,y);
//			sca_vec_mut(alpha,delta_lambda,m,QP_temp1);
//			vec_add(lambda,QP_temp1,m,lambda);
//
//			//printf("x is :\n");show_matrix(x,1,n);putchar('\n');
//			//printf("y is :\n");show_matrix(y,1,m);putchar('\n');
//			//printf("lambda is :\n");show_matrix(lambda,1,m);putchar('\n');
//
//			// Iteration number + 1
//			Iter = Iter + 1;
//			//printf("Iteration: %d\n",Iter);
//
//			mu = dot_product(y,lambda,m)/m;
//
//			sigma = mu/mu_old;
//			sigma = sigma*sigma*sigma;
//			sigma = scalar_min(sigma,0.99999);
//
//			end_time_section = Timestamp_get32();
//			total_time_section = end_time_section - start_time_section - dead_time;
//			Time_Section4 +=(double)total_time_section;
//
//			start_time_section = Timestamp_get32();
//
//			//if(res<0.3)
//			//	tau = 1-res;
//
//			//// Calculate KKT residual.
//			//mat_vec(G,x,n,n,QP_temp1);
//			//mat_vec(A_t,lambda,n,m,QP_temp3);
//			//vec_sub(QP_temp3,c,n,QP_temp2);
//			//vec_sub(QP_temp1,QP_temp2,n,QP_temp3);
//			//res_1 = scalar_max(scalar_abs(vec_min(QP_temp3,n)),scalar_abs(vec_max(QP_temp3,n)));
//			//mat_vec(A,x,m,n,QP_temp1);
//			//vec_add(y,b,m,QP_temp2);
//			//vec_sub(QP_temp1,QP_temp2,n,QP_temp3);
//			//res_2 = scalar_max(scalar_abs(vec_min(QP_temp3,n)),scalar_abs(vec_max(QP_temp3,n)));
//			//res = scalar_max(scalar_max(res_1,res_2),scalar_abs(mu));
//			////printf("Res_1: %f\n",res_1);
//			////printf("Res_2: %f\n",res_2);
//			////printf("mu : %f\n",abs(mu));
//			////printf("Residual: %f\n",res);
//
//
//			// KKT residual (Wright) Termination Criteria
//			//flag = TC_KKT_Wright(rd, rp, b_max,c_max, mu, 0.00001,m, n);
//			//if(flag==1)
//			//	break;
//
//			// Himmelblau Termination Criteria
//			mat_vec(G,x,n,n,QP_temp1);
//			Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
//			res_1 = scalar_abs(Fx_old-Fx);
//			vec_sub(x_old,x,n,QP_temp1);
//			vec_abs(QP_temp1,n,QP_temp2);
//			res_2 = vec_max(QP_temp2,n);
//			mat_vec(A,x,m,n,QP_temp1);
//			vec_sub(QP_temp1,b,m,QP_temp2);
//			cons_obey = vec_min(QP_temp2,m);
//			res = scalar_max(res_1,res_2);
//			//printf("The Fx_old is: %f \n",Fx_old);
//			//printf("The Fx is: %f \n",Fx);
//			// Update
//			Fx_old = Fx;
//			for(k_in=0;k_in<n;k_in++)
//				x_old[k_in] = x[k_in];
//			//printf("The res_1 is: %f \n",res_1);
//			//printf("The res_2 is: %f \n",res_2);
//			//printf("cons_obey is: %f \n",cons_obey);
//
//			//Terminal condition(Himmelblau)
//			if(res<0.00001&&cons_obey>=-0.00001)
//				break;
//
//
//			//// Show new Fx value
//			//mat_mul(x,1,n,G,n,QP_temp1);
//			//Fx = 0.5*dot_product(temp1,x,n)+dot_product(c,x,n);
//			//printf("This iteration F(x) value is: %f\n\n",Fx);
//
//			////Terminal condition(KKT Residuel)
//			//if(res<0.00001)
//			//	break;
//
//
//			end_time_section = Timestamp_get32();
//			total_time_section = end_time_section - start_time_section - dead_time;
//			Time_Section5 +=(double)total_time_section;
//
//		}
//		memcpy(delta_u,x,n*sizeof(double));
//		QP_Iteration = Iter;
//
//	}
//
//	//printf("Total: Iteration: %d\n",Iter);
//
//
//	QP_Iter_end_time = Timestamp_get32();
//	Time_QP_Iter = (double)QP_Iter_end_time-QP_Iter_start_time-dead_time;
//
//}
//
//
//
////===============================priduip begin===============================================
//
////V4 modification:
////	1 Himmelblau Termination Criteria
////	2 Initial Value Reset
//void priduip_v4(double *G, double *GL,double *GL_T, double *c, double *A, double *A_t, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u)
//{
//
//	double Fx_old, Fx, cons_obey;
//	double res = 0.3;
//	double tau = 0.7;
//	int Iter = 0;
//	int k, k_in;
//	double mu,mu_aff;
//	double alpha_aff_pri = 0.999995;
//	double alpha_aff_dual = 0.999995;
//	double alpha_aff = 0.99995;
//	double alpha_cond;	//Use to check whether the alpha is valid.
//	double sigma;
//	double alpha_tau_pri = 0.999995;
//	double alpha_tau_dual = 0.999995;
//	double alpha;
//	double res_1,res_2;
//	int flag = 0;
//
//	//double *diag;
//	//double *A_t;
//
//	UInt32 start_time,end_time,dead_time,total_time;
//	UInt32 start_time_ls1,end_time_ls1,total_time_ls1;
//	UInt32 start_time_ls2,end_time_ls2,total_time_ls2;
//	UInt32 start_time_section,end_time_section,total_time_section;
//
//	start_time = Timestamp_get32();
//	end_time = Timestamp_get32();
//	dead_time = end_time - start_time;
//
//	//Make A transpose as a variable
//
//	////Show the problem paremeters.
//	show_matrix(G,n,n);putchar('\n');
//	show_matrix(c,1,n);putchar('\n');
//	show_matrix(A,n,m);putchar('\n');
//	show_matrix(b,1,m);putchar('\n');
//	show_matrix(x,1,n);putchar('\n');
//	show_matrix(y,1,m);putchar('\n');
//	show_matrix(lambda,1,m);putchar('\n');
//
//	//Initial Fx_old value
//	mat_vec(G,x,n,n,QP_temp1);
//	Fx_old = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
//	//printf("This initial F(x) value is: %f\n",Fx_old);
//
//	//Initial x_old value
//	#pragma MUST_ITERATE (3);
//	for(k=0;k<n;k++)
//		x_old[k] = x[k];
//
//	//Check if the Uncontrained solution meet the constraints
//	//vec_rev(c,n,c_rev);
//	luEvaluate(GL,GL_T,c_rev,n,x_try);
//	mat_vec(A,x_try,m,n,QP_temp1);
//	vec_sub(QP_temp1,b,m,QP_temp2);
//	if((vec_min(QP_temp2,m))>0)
//	{
//		memcpy(delta_u,x_try,n*sizeof(double));
//		QP_Iteration = Iter;
//	}
//	else
//	{
//		// If not global optimal, do iterative process.
//		for(k=1;k<=MaxIteration;k++)
//		{
//			//start_time = clock();
//			start_time_section = Timestamp_get32();
//
//			mat_vec(G,x,n,n,QP_temp1);
//			mat_vec(A_t,lambda,n,m,QP_temp2);
//			vec_sub(QP_temp1, QP_temp2, n,QP_temp3);
//			vec_add(QP_temp3,c,n,rd);
//			printf("rd is:\n");show_matrix(rd,1,n);putchar('\n');
//
//			mat_vec(A,x,m,n,QP_temp1);
//			vec_sub(b,QP_temp1,m,rp_y);
//			vec_sub(QP_temp1,y,m,QP_temp2);
//			vec_sub(QP_temp2,b,m,rp);
//			printf("rp is:\n");show_matrix(rp,1,m);putchar('\n');
//
//			vec_div(lambda,y,m,inv_y_lambda);
//			bsxfun(A,inv_y_lambda,m,n,y_lambda_A);
//			printf("inv_y_lambda is :\n");show_matrix(inv_y_lambda,1,m);putchar('\n');
//			printf("y_lambda_A is :\n");show_matrix(y_lambda_A,n,m);putchar('\n');
//
//			//end_time_section = Timestamp_get32();
//			////total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section1 +=(double)total_time_section;
//
//
//			//start_time_section = Timestamp_get32();
//			mat_mul(A_t,n,m,y_lambda_A,n,QP_temp1);	//Time Consuming! Much better if custom!
//			mat_add(G,QP_temp1,n,n,equation_F);
//			printf("equation_F is :\n");show_matrix(equation_F,n,n);putchar('\n');
//
//			vec_mul(inv_y_lambda,rp_y,m,QP_temp1);
//			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
//			vec_sub(QP_temp2,rd,n,equation_b);
//			printf("equation_b is: \n");show_matrix(equation_b,1,n);
//
//			//start_time_ls1 = Timestamp_get32();
//			line_solve(equation_F,equation_b,n,delta_x_aff);
//
//			//end_time_ls1 = Timestamp_get32();
//			//total_time_ls1 = end_time_ls1 - start_time_ls1 - dead_time;
//			//Time_LineSolve1+=(double)total_time_ls1;
//			printf("delta_x_aff is: \n");show_matrix(delta_x_aff,1,n);putchar('\n');
//
//			mat_vec(A,delta_x_aff,m,n,QP_temp1);
//			vec_add(QP_temp1,rp,m,delta_y_aff);
//			printf("delta_y_aff is: \n");show_matrix(delta_y_aff,1,m);putchar('\n');
//
//			vec_mul(inv_y_lambda,delta_y_aff,m,QP_temp1);
//			vec_add(QP_temp1,lambda,m,QP_temp2);
//			vec_rev(QP_temp2,m,delta_lambda_aff);
//			printf("delta_lambda_aff is: \n");show_matrix(delta_lambda_aff,1,m);putchar('\n');
//
//			//end_time_section = Timestamp_get32();
//			//total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section2 +=(double)total_time_section;
//
//
//			//start_time_section = Timestamp_get32();
//
//			mu = dot_product(y,lambda,m)/m;
//
//			//end_time = clock();
//			//total_time = end_time - start_time - dead_time;
//			//Time_QP_step1[0]+=(double)total_time;
//
//			alpha_aff_pri = alpha_decreas(y,delta_y_aff,m,1);
//			alpha_aff_dual = alpha_decreas(lambda,delta_lambda_aff,m,1);
//			//alpha_aff_pri = alpha_decreas_old(y,delta_y_aff,m,1);
//			//alpha_aff_dual = alpha_decreas_old(lambda,delta_lambda_aff,m,1);
//			alpha_aff = scalar_min(alpha_aff_pri,alpha_aff_dual);
//			printf("alpha_aff is : %f\n",alpha_aff);
//
//			//end_time_section = Timestamp_get32();
//			//total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section3 +=(double)total_time_section;
//
//			//start_time_section = Timestamp_get32();
//
//			sca_vec_mut(alpha_aff,delta_y_aff,m,QP_temp3);
//			vec_add(y,QP_temp3,m,QP_temp1);
//			sca_vec_mut(alpha_aff,delta_lambda_aff,m,QP_temp3);
//			vec_add(lambda,QP_temp3,m,QP_temp2);
//			mu_aff = dot_product(QP_temp1,QP_temp2,m)/m;
//
//			sigma = mu_aff/mu;
//			sigma = sigma*sigma*sigma;
//
//			//end_time = clock();
//			//total_time = end_time - start_time - dead_time;
//			//Time_QP_step2[0]+=(double)total_time;
//
//			//start_time = clock();
//
//			vec_mul(delta_lambda_aff,delta_y_aff,m,QP_temp1);
//			sca_vec_sub(sigma*mu,QP_temp1,m,sig_mu_lam_y);
//			vec_div(sig_mu_lam_y,lambda,m,QP_temp2);
//			vec_add(rp_y,QP_temp2,m,aff_part);
//			vec_mul(inv_y_lambda,aff_part,m,QP_temp1);
//			mat_vec(A_t,QP_temp1,n,m,QP_temp2);
//			vec_sub(QP_temp2,rd,n,equation_b);
//			printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');
//
//			//start_time_ls2 = clock();
//			line_solve(equation_F,equation_b,n,delta_x);
//
//			//end_time_ls2 = clock();
//			//total_time_ls2 = end_time_ls2 - start_time_ls2 - dead_time;
//			//Time_LineSolve2+=(double)total_time_ls2;
//
//			printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');
//
//			mat_vec(A,delta_x,m,n,QP_temp1);
//			vec_add(QP_temp1,rp,m,delta_y);
//			printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');
//
//			vec_mul(lambda,delta_y,m,QP_temp1);
//			vec_mul(delta_lambda_aff,delta_y_aff,m,QP_temp2);
//			vec_add(QP_temp1,QP_temp2,m,QP_temp3);
//			sca_vec_sub(sigma*mu,QP_temp3,m,QP_temp1);
//			vec_div(QP_temp1,y,m,QP_temp2);
//			vec_sub(QP_temp2,lambda,m,delta_lambda);
//			printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');
//
//			//end_time_section = Timestamp_get32();
//			//total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section4 +=(double)total_time_section;
//
//			//start_time_section = Timestamp_get32();
//
//			if(res<0.3)
//				tau = 1-res;
//
//			//end_time = Timestamp_get32();
//			//total_time = end_time - start_time - dead_time;
//			//Time_QP_step3[0]+=(double)total_time;
//
//			//start_time = clock();
//
//			alpha_tau_pri = alpha_decreas(y,delta_y,m,tau);
//			alpha_tau_dual = alpha_decreas(lambda,delta_lambda,m,tau);
//			//alpha_tau_pri = alpha_decreas_old(y,delta_y,m,tau);
//			//alpha_tau_dual = alpha_decreas_old(lambda,delta_lambda,m,tau);
//			alpha = scalar_min(alpha_tau_pri,alpha_tau_dual);
//			printf("alpha is : %f\n",alpha);
//
//			//end_time_section = Timestamp_get32();
//			//total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section5 +=(double)total_time_section;
//
//
//			//start_time_section = Timestamp_get32();
//
//			//Obtain new iteration value.
//			sca_vec_mut(alpha,delta_x,n,QP_temp1);
//			vec_add(x,QP_temp1,n,x);
//			sca_vec_mut(alpha,delta_y,m,QP_temp1);
//			vec_add(y,QP_temp1,m,y);
//			sca_vec_mut(alpha,delta_lambda,m,QP_temp1);
//			vec_add(lambda,QP_temp1,m,lambda);
//
//			printf("x is :\n");show_matrix(x,1,n);putchar('\n');
//			printf("y is :\n");show_matrix(y,1,m);putchar('\n');
//			printf("lambda is :\n");show_matrix(lambda,1,m);putchar('\n');
//
//			// Iteration number + 1
//			Iter = Iter + 1;
//			//printf("Iteration: %d\n",Iter);
//
//			//end_time_section = Timestamp_get32();
//			//total_time_section = end_time_section - start_time_section - dead_time;
//			//Time_Section6 +=(double)total_time_section;
//
//			//start_time_section = Timestamp_get32();
//
//			// KKT residual (Wright) Termination Criteria
//			//printf("b_max,c_max is:%f,%f\n\n",b_max,c_max);
//			//flag = TC_KKT_Wright(rd, rp, b_max,c_max, mu, 0.0001,m, n);
//			//if(flag==1)
//			//	break;
//
//
//			//// Himmelblau Termination Criteria
//			//mat_vec(G,x,n,n,QP_temp1);
//			//Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
//			//res_1 = scalar_abs(Fx_old-Fx);
//			//vec_sub(x_old,x,n,QP_temp1);
//			//vec_abs(QP_temp1,n,QP_temp2);
//			//res_2 = vec_max(QP_temp2,n);
//			//mat_vec(A,x,m,n,QP_temp1);
//			//vec_sub(QP_temp1,b,m,QP_temp2);
//			//cons_obey = vec_min(QP_temp2,m);
//			//res = scalar_max(res_1,res_2);
//			printf("The Fx_old is: %f \n",Fx_old);
//			//printf("The Fx is: %f \n",Fx);
//			// Update
//			//Fx_old = Fx;
//			//for(k_in=0;k_in<n;k_in++)
//			//	x_old[k_in] = x[k_in];
//			//printf("The res_1 is: %f \n",res_1);
//			//printf("The res_2 is: %f \n",res_2);
//			//printf("cons_obey is: %f \n",cons_obey);
//			//Terminal condition(Himmelblau)
//			//if(res<0.00001&&cons_obey>=-0.00001)
//			//	break;
//
//
//			//end_time = clock();
//			//total_time = end_time - start_time - dead_time;
//			//Time_QP_step4[0]+=(double)total_time;
//
//			//// Show new Fx value
//			//mat_mul(x,1,n,G,n,QP_temp1);
//			//Fx = 0.5*dot_product(QP_temp1,x,n)+dot_product(c,x,n);
//			//printf("This iteration F(x) value is: %f\n\n",Fx);
//
//			////Terminal condition(KKT Residuel)
//			//if(res<0.00001)
//			//	break;
//
//
//		}
//		memcpy(delta_u,x,n*sizeof(double));
//		QP_Iteration = Iter;
//
//		//end_time_section = Timestamp_get32();
//		//total_time_section = end_time_section - start_time_section - dead_time;
//		//Time_Section7 +=(double)total_time_section;
//	}
//}
//
////================================priduip end======================================================
//
//
//
////MPC online computing!
////Yi Ding @2014.4.18
////For priduip_v4 @2014.7.9
//void mpc_DP(double *delta_u, double *delta_u_ini,double *y_ini,double *lambda_ini, double *u_k, double *y_k,double *x_k,double *r_k, double *u_k_1, double *u_k_2)
//{
//
//	int i;
//	double WarmVal = 0.5;
//
//	UInt32 start_time,end_time,dead_time;
//	UInt32 QP_start_time,QP_end_time,QP_total_time;
//
//    start_time = Timestamp_get32();
//    end_time = Timestamp_get32();
//    dead_time = end_time - start_time;
//
//	feedback_v3(x_k,y_k,u_k_1,u_k_2,x_k,L);
//	//show_matrix(x_k,1,nx+ny);
//
//	vec_aug(u_k_1,nu,M,aug_u_k_1);
//
//	////Below form the omega_r with y constraints.
//	//omega_r = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//	//n_omega_r = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//	//omega_r_form(aug_u_k_1, x_k, omega_r,F);
//	//vec_rev(omega_r,4*nu*M+2*ny*P,n_omega_r);
//	////show_matrix(omega_r,1,4*nu*M+2*ny*P);
//
//	//Below form the omega_r without y constraints.
//	omega_r_form2(aug_u_k_1, x_k, omega_r,F);
//	vec_rev(omega_r,4*nu*M,n_omega_r);
//	//show_matrix(omega_r,1,4*nu*M+2*ny*P);
//
//	//Below form the c in QP.
//	mat_vec(F,x_k,ny*P,nx+ny,MPC_temp1);
//	vec_sub(MPC_temp1,r_k,ny*P,MPC_temp2);
//	sca_vec_mut(Q,MPC_temp2,ny*P,MPC_temp3);
//	vec_mat(Phi,MPC_temp3,ny*P,nu*M,c);
//	vec_rev(c,nu*M,c_rev);
//	//printf("c is:\n");show_matrix(c,1,nu*M);putchar('\n');
//	//printf("c_rev is:\n");show_matrix(c_rev,1,nu*M);putchar('\n');
//
//	//QP Starting Point initializing (delta_u)
//	for(i=0;i<nu*M;i++)
//		delta_u_ini[i] = delta_u_ini[i+1];
//	delta_u_ini[nu*M-1] = 0;
//
//	// Below is a parameter used in an initializing method from LOQO
//	//mat_vec(OMEGA_L,delta_u_ini,4*nu*M+2*ny*P,nu*M,MPC_temp4);
//	//vec_sub(omega_r,temp4,4*nu*M+2*ny*P,h_x);
//
//	//QP Starting Point initializing (A heuristic method for y and lambda)
//	for(i=0;i<4*nu*M+2*ny*P;i++)
//	{
//		y_ini[i] = scalar_max(y_ini[i],0.5);
//		lambda_ini[i]  = scalar_max(lambda_ini[i],0.5);
//	}
//
//	//QP Starting Point initializing (A simple method for y and lambda)
//	//for(i=0;i<4*nu*M+2*ny*P;i++)
//	//{
//	//	y_ini[i] = 0.5;
//	//	lambda_ini[i]  = 0.5;
//	//}
//
//
//	////QP Starting Point initializing (AUT's method)
//	// Scale preparation
//	vec_abs(n_omega_r,4*nu*M,MPC_temp1);
//	b_max = vec_max(MPC_temp1,4*nu*M);
//	vec_abs(c,nu*M,MPC_temp1);
//	c_max = vec_max(MPC_temp1,nu*M);
//	vec_abs(G,nu*M*nu*M,MPC_temp1);
//	G_max = vec_max(MPC_temp1,nu*M*nu*M);
//	vec_abs(N_OMEGA_L,4*nu*M*nu*M,MPC_temp1);
//	A_max = vec_max(MPC_temp1,4*nu*M*nu*M);
//	p_max = scalar_max(scalar_max(b_max,c_max),scalar_max(G_max,A_max));
//	//printf("b_max, c_max is :%f,%f\n",b_max,c_max);
//	//printf("p_max: %f\n",p_max);
//	//if(p_max>1)
//	//	WarmVal = sqrt(p_max);
//	//for(i=0;i<4*nu*M+2*ny*P;i++)
//	//{
//	//	y_ini[i] = WarmVal;
//	//	lambda_ini[i]  = WarmVal;
//	//}
//
//	////QP Starting Point initializing (Wright's method)
//	//SP_wright(delta_u_ini, y_ini, lambda_ini, G, c, N_OMEGA_L,N_OMEGA_L_T, n_omega_r,y_ini, lambda_ini);
//	//printf("delta_u_ini is: \n");show_matrix(delta_u_ini,1,nu*M);putchar('\n');
//	//printf("y_ini is: \n");show_matrix(y_ini,1,4*nu*M);putchar('\n');
//	//printf("lambda_ini is: \n");show_matrix(lambda_ini,1,4*nu*M);putchar('\n');
//
//	////The below is a way to initialize lambda with a random value between 0 and 1.
//	////In the process of debug with Matlab, the lambda_ini is define just like which in Matlab
//	//lambda_ini = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//	//for(i=0;i<4*nu*M+2*ny*P;i++)
//	//	lambda_ini[i] = ((double)rand()/RAND_MAX);
//	//Show initial values.
//	//putchar('\n');show_matrix(delta_u_ini,1,nu*M);
//	//putchar('\n');show_matrix(y_ini,1,4*nu*M);
//	//putchar('\n');show_matrix(lambda_ini,1,4*nu*M);
//	// This random method can be used to compare with other method.
//
//
//	////QP solve
//	QP_start_time = Timestamp_get32();
//	//IPM_V3_SP(G,GL,GL_T,c, N_OMEGA_L,N_OMEGA_L_T, n_omega_r, delta_u_ini, y_ini, lambda_ini, 4*nu*M, nu*M, delta_u_M);
//	priduip_v4(G, GL,GL_T, c, N_OMEGA_L, N_OMEGA_L_T, n_omega_r, delta_u_ini, y_ini, lambda_ini, 4*nu*M, nu*M, delta_u_M);
//	//printf("delta_u is: \n");show_matrix(delta_u_M,1,nu*M);
//	memcpy(delta_u,delta_u_M,nu*sizeof(double));
//	QP_end_time = Timestamp_get32();
//	QP_total_time = QP_end_time - QP_start_time - dead_time;
//	Time_QP = (double)(QP_total_time);
//
//	//vec_add(u_k,delta_u,nu,MPC_temp1);
//	for(i=0;i<nu;i++)
//	{
//		MPC_temp1[i]=u_k[i]+delta_u[i];
//	}
//	memcpy(u_k,MPC_temp1,nu*sizeof(double));
//}
//
///*
// *  ======== taskFxn ========
// */
//Void taskFxn(UArg a0, UArg a1)
//{
//    printf("enter taskFxn()\n");
//    int i,j;
//	//double *u_k,*delta_u,*x_k,*y_k,*r_k;
//	//double *u_k,*u_k_1,*u_k_2,*delta_u,*x_k,*y_k,*xm;
//
//	//The variables below are used to test time consuming.
//	double u_k[] = {0};//double u_k[] = {1};
//	double u_k_1[] = {0};//double u_k_1[] = {1};
//	double u_k_2[] = {0};
//	double delta_u[] = {0.1};//double delta_u[] = {1};
//	double x_k[] = {0,0,0};
//	double y_k[] = {0};
//	double xm[] = {0,0};//double xm[] = {0,0.0787};
//
//	//Reference curve
//	double r_k1[] = {3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3,3};
//	double r_k2[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};
//
//	UInt32 start_time,end_time,dead_time;
//	UInt32 MPC_start_time,MPC_end_time,MPC_total_time;
//	Types_FreqHz freq;
//
//	FILE *fp_u, *fp_y, *fp_t, *fp_i,*fp_mpct;
//	FILE *fp_qpt_nonIter,*fp_qpt_Iter;
//	FILE *fp_s1t, *fp_s2t, *fp_s3t, *fp_s4t, *fp_s5t;
//
//	//Determine the number of inequality constraints and number of decision values
//	ndec = nu*M;
//	mc = 4*nu*M;	//Outputs constraints not included.
//	//mc = 4*nu*M+2*ny*P;	//Outputs constraints included.
//
//	// Memory allocation for main function
//	L = (double*)calloc((nx+ny)*ny,sizeof(double));
//	F = (double*)calloc(P*ny*(nx+ny),sizeof(double));
//	Phi = (double*)calloc(M*nu*P*ny,sizeof(double));
//	B = (double*)calloc(P*nu*M*nu*M,sizeof(double));
//	OMEGA_L = (double*)calloc((4*nu*M)*nu*M,sizeof(double));
//	N_OMEGA_L = (double*)calloc((4*nu*M)*nu*M,sizeof(double));
//	N_OMEGA_L_T = (double*)calloc((4*nu*M)*nu*M,sizeof(double));
//	G = (double*)calloc(nu*M*nu*M,sizeof(double));
//	GL = (double*)calloc(nu*M*nu*M,sizeof(double));
//	GL_T = (double*)calloc(nu*M*nu*M,sizeof(double));
//	delta_u_ini = (double*)calloc(nu*M,sizeof(double));
//	y_ini = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//	lambda_ini = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//
//	// Memory allocation for mpc_DP
//	aug_u_k_1 = (double*)calloc(nu*M,sizeof(double));
//	omega_r = (double*)calloc(4*nu*M,sizeof(double));
//	n_omega_r = (double*)calloc(4*nu*M,sizeof(double));
//	c = (double*)calloc(nu*M,sizeof(double));
//	c_rev = (double*)calloc(nu*M,sizeof(double));
//	//h_x = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//	delta_u_M = (double*)calloc(nu*M,sizeof(double));
//	MPC_temp1 = (double*)calloc(ny*P*ny*P,sizeof(double));
//	MPC_temp2 = (double*)calloc(ny*P*ny*P,sizeof(double));
//	MPC_temp3 = (double*)calloc(ny*P*ny*P,sizeof(double));
//	MPC_temp4 = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
//
//	//Memory allocation for IPM_V3_DP
//	//Error may emerge if the temp variable exceed the allocated area!
//	x_old = (double*)calloc(nu*M,sizeof(double));
//	rd = (double*)calloc(nu*M,sizeof(double));
//	rp = (double*)calloc(4*nu*M,sizeof(double));
//	rp_y = (double*)calloc(4*nu*M,sizeof(double));
//	inv_y_lambda = (double*)calloc(4*nu*M,sizeof(double));
//	y_lambda_A = (double*)calloc(4*nu*M*nu*M,sizeof(double));
//	equation_F = (double*)calloc(nu*M*nu*M,sizeof(double));
//	equation_b = (double*)calloc(nu*M,sizeof(double));
//	delta_x = (double*)calloc(nu*M,sizeof(double));
//	delta_y = (double*)calloc(4*nu*M,sizeof(double));
//	delta_lambda = (double*)calloc(4*nu*M,sizeof(double));
//	temp_cond = (double*)calloc(4*nu*M,sizeof(double));
//	x_try = (double*)calloc(nu*M,sizeof(double));
//	sig_mu = (double*)calloc(4*nu*M,sizeof(double));
//	center_part = (double*)calloc(4*nu*M,sizeof(double));
//	QP_temp1 = (double*)calloc(4*nu*M*4*nu*M,sizeof(double));
//	QP_temp2 = (double*)calloc(4*nu*M*4*nu*M,sizeof(double));
//	QP_temp3 = (double*)calloc(4*nu*M*4*nu*M,sizeof(double));
//
//	//Additional Memory allocation for priduip_v4
//	//Error may emerge if the temp variable exceed the allocated area!
//	delta_x_aff = (double*)calloc(nu*M,sizeof(double));
//	delta_y_aff = (double*)calloc(4*nu*M,sizeof(double));
//	delta_lambda_aff = (double*)calloc(4*nu*M,sizeof(double));
//	//diag = (double*)calloc(m*m,sizeof(double));
//	sig_mu_lam_y = (double*)calloc(4*nu*M,sizeof(double));
//	aff_part = (double*)calloc(4*nu*M,sizeof(double));
//
//	//Memory allocation for feedback_v3
//	fb_temp1 = (double*)calloc(nx+ny,sizeof(double));
//	fb_temp2 = (double*)calloc(nx+ny,sizeof(double));
//
//	//Memory allocation for mcf
//	mcf_C = (double*)calloc(4*nu*M*4*nu*M+1,sizeof(double));
//	mcf_temp_vec = (double*)calloc(4*nu*M+1,sizeof(double));
//
//	//Memory allocation for alpha_decreas
//	alpha_decreas_neg_delta_value = (double*)calloc(4*nu*M,sizeof(double));
//	alpha_decreas_neg_vec_value = (double*)calloc(4*nu*M,sizeof(double));
//	alpha_decreas_temp1 = (double*)calloc(4*nu*M,sizeof(double));
//	alpha_decreas_temp2 = (double*)calloc(4*nu*M,sizeof(double));
//
//	//Memory allocation for line_solve
//	line_solve_L = (double*)calloc(nu*M * nu*M, sizeof(double));
//	line_solve_p = (double*)calloc(nu*M, sizeof(double));
//
//	//printf("ndec is: %d\n",ndec);
//
//	Timestamp_getFreq(&freq);
//	//printf("freq = {0x%x : 0x%x}\n", freq.hi, freq.lo);
//	printf("Frequency is: %u\n",freq.lo);
//
//	start_time = Timestamp_get32();
//	end_time = Timestamp_get32();
//	dead_time = end_time - start_time;
//
//	//Initialization of Feedback gain coefficient L.
//	j = ny+1;
//	for(i=nx*ny;i<(nx+ny)*ny;i+=j){
//		L[i] = 1;
//		j++;
//	}
//	//show_matrix(L,1,(nx+ny)*ny);
//
//	//Below form the MPC parameters F, Phi and B.
//	fphi(F,Phi);
//	form_B(B);
//	//show_matrix(B,nu*M,nu*M);
//
//	////Below form the OMEGA_L with y constraints
//	//OMEGA_L = (double*)calloc((4*nu*M+2*ny*P)*nu*M,sizeof(double));
//	//N_OMEGA_L = (double*)calloc((4*nu*M+2*ny*P)*nu*M,sizeof(double));
//	//OMEGA_L_form(B,Phi,OMEGA_L);
//	//vec_rev(OMEGA_L,(4*nu*M+2*ny*P)*nu*M,N_OMEGA_L);
//	////show_matrix(OMEGA_L,nu*M,(4*nu*M+2*ny*P));
//
//	//Below form the OMEGA_L without y constraints
//	OMEGA_L_form2(B,Phi,OMEGA_L);
//	vec_rev(OMEGA_L,(4*nu*M)*nu*M,N_OMEGA_L);
//	memcpy(N_OMEGA_L_T,N_OMEGA_L,(4*nu*M)*nu*M*sizeof(double));
//	transpose(N_OMEGA_L_T, nu*M, (4*nu*M));
//	//show_matrix(OMEGA_L,nu*M,(4*nu*M+2*ny*P));
//
//	//Below form the QP parameters G, GL and GL_T.
//	form_G(Phi,G);
//	cholesky(G,nu*M,GL);
//	memcpy(GL_T,GL,nu*M*nu*M*sizeof(double));
//	transpose(GL_T, nu*M, nu*M);
//	//show_matrix(G,nu*M,nu*M);putchar('\n');
//	//show_matrix(GL,nu*M,nu*M);putchar('\n');
//	//show_matrix(GL_T,nu*M,nu*M);putchar('\n');
//
//	//Initialization of QP parameters in main function.
//	for(i=0;i<nu*M;i++)
//		delta_u_ini[i] = 0.5;
//	for(i=0;i<4*nu*M+2*ny*P;i++)
//	{
//		y_ini[i] = 0.5;
//		lambda_ini[i] = 0.5;
//	}
//
//	// Runtime data record preparation.
//	fp_u=fopen("DSP_u_k.txt","w");fp_y=fopen("DSP_y_k.txt","w");
//	fp_t=fopen("DSP_QPtime.txt","w");fp_i=fopen("DSP_QPIter.txt","w");
//	fp_mpct=fopen("DSP_MPCtime.txt","w");
//	fp_qpt_nonIter = fopen("Time_QP_non_Iter.txt","w");
//	fp_qpt_Iter = fopen("Time_QP_Iter.txt","w");
//	fp_s1t = fopen("Time_QP_Section1.txt","w");
//	fp_s2t = fopen("Time_QP_Section2.txt","w");
//	fp_s3t = fopen("Time_QP_Section3.txt","w");
//	fp_s4t = fopen("Time_QP_Section4.txt","w");
//	fp_s5t = fopen("Time_QP_Section5.txt","w");
//
//	for(i=0;i<200;i++)
//	{
//		MPC_start_time = Timestamp_get32();
//
//		if(i<100)
//			mpc_SP(delta_u,delta_u_ini,y_ini,lambda_ini,u_k,y_k,x_k,r_k1,u_k_1,u_k_2);
//		else
//			mpc_SP(delta_u,delta_u_ini,y_ini,lambda_ini,u_k,y_k,x_k,r_k2,u_k_1,u_k_2);
//		//Here the DSP communicate with Matlab to transmite u(k) and get y(k)
//
//		MPC_end_time = Timestamp_get32();
//		MPC_total_time = MPC_end_time - MPC_start_time - dead_time;
//		Time_MPC = (double)(MPC_total_time);
//		fprintf(fp_mpct,"%f\n",Time_MPC/freq.lo*1000);
//
//
//		//Here we use a simple model instead.
//		process_model(xm,u_k,y_k);
//
//		memcpy(u_k_2,u_k_1,nu*sizeof(double));
//		memcpy(u_k_1,u_k,nu*sizeof(double));
//
//		//printf("Manipulated Variables	(u(k)):\n");
//		//show_matrix(u_k,1,nu);putchar('\n');
//		fprintf(fp_u,"%f\n",u_k[0]);
//		//printf("Manipulated Variables	(u(k-1)):\n");
//		//show_matrix(u_k_1,1,nu);putchar('\n');
//		//printf("Manipulated Variables	(u(k-2)):\n");
//		//show_matrix(u_k_2,1,nu);putchar('\n');
//		//printf("Aug State Variables	(x(k)):\n");
//		//show_matrix(x_k,1,nx+ny);putchar('\n');
//		//printf("State Variables		(xm(k)):\n");
//		//show_matrix(xm,1,nx);putchar('\n');
//		//printf("Controlled Variables	(y(k)):\n");
//		//show_matrix(y_k,1,ny);putchar('\n');putchar('\n');
//		fprintf(fp_y,"%f\n",y_k[0]);
//		fprintf(fp_t,"%f\n",Time_QP/freq.lo*1000);
//		fprintf(fp_i,"%f\n",QP_Iteration);
//		fprintf(fp_qpt_nonIter,"%f\n",Time_QP_non_Iter/freq.lo*1000);
//		fprintf(fp_qpt_Iter,"%f\n",Time_QP_Iter/freq.lo*1000);
//		fprintf(fp_s1t,"%f\n",Time_Section1/freq.lo*1000);
//		fprintf(fp_s2t,"%f\n",Time_Section2/freq.lo*1000);
//		fprintf(fp_s3t,"%f\n",Time_Section3/freq.lo*1000);
//		fprintf(fp_s4t,"%f\n",Time_Section4/freq.lo*1000);
//		fprintf(fp_s5t,"%f\n",Time_Section5/freq.lo*1000);
//
//		Time_QP = 0;QP_Iteration = 0;Time_QP_non_Iter = 0;Time_QP_Iter = 0;
//		Time_Section1 = 0;Time_Section2 = 0;Time_Section3 = 0;Time_Section4 = 0;Time_Section5 = 0;
//
//		//
//	}
//
//	fclose(fp_u);fclose(fp_y);
//
//	//end_time = clock();
//	//total_time = end_time - start_time - dead_time;
//	//Time_MPC = (double)total_time;
//	//printf("%d times duration: %f\n",i,(double)total_time);
//
//	//printf("MPC time consuming: %f\n",Time_MPC);
//	//printf("QP solving time consuming(Clocks): %f\n",Time_QP/freq.lo);
//	//printf("QP solving time consuming(ms Per Single time): %f\n",((Time_QP/freq.lo)/i)*1000);
//	//printf("QP solving step 1 time consuming: %f\n",Time_QP_step1);
//	//printf("QP solving step 2 time consuming: %f\n",Time_QP_step2);
//	//printf("QP solving step 3 time consuming: %f\n",Time_QP_step3);
//	//printf("QP solving step 4 time consuming: %f\n",Time_QP_step4);
//	//printf("Linear system1 solving time consuming(ms): %f\n",(Time_LineSolve1/freq.lo)*1000);
//	//printf("Linear system1 solving times: %d\n",LS_solve_times);
//	//printf("Linear system1 average solving time(ms): %f\n",(Time_LineSolve1/freq.lo)*1000/LS_solve_times);
//	//printf("Linear system2 solving time cnsuming: %f\n\n",Time_LineSolve2);
//
//	//printf("QP Section 1 time consuming: %f\n",Time_Section1);
//	//printf("QP Section 2 time consuming: %f\n",Time_Section2);
//	//printf("QP Section 3 time consuming(ms): %f\n",(Time_Section3/freq.lo)*1000);
//	//printf("alpha_compute_times:%d\n",alpha_compute_times);
//	//printf("QP Section 4 time consuming: %f\n",Time_Section4);
//	//printf("QP Section 5 time consuming: %f\n",Time_Section5);
//	//printf("QP Section 6 time consuming: %f\n",Time_Section6);
//	//printf("QP Section 7 time consuming: %f\n",Time_Section7);
//	//printf("QP Section 8 time consuming: %f\n\n",(Time_Section8/freq.lo)*1000);
//	//printf("Mat_vec solving times: %d\n",Mat_vec_times);
//	//printf("Mat_vec_times average solving time(ms): %f\n",(Time_Section8/freq.lo)*1000/Mat_vec_times);
//
//	//printf("dot_product time consuming (ms): %f\n",(Time_bsxfun/freq.lo)*1000);
//	//printf("dot_product solving times: %d\n",dot_product_times);
//	//printf("dot_product average solving time(ms): %f\n",(Time_dot_product/freq.lo)*1000/dot_product_times);
//
//
//	//printf("vec_add time consuming (ms): %f\n",(Time_vec_add/freq.lo)*1000);
//	//printf("vec_add solving times: %d\n",vec_add_times);
//	//printf("vec_add average solving time(ms): %f\n\n",(Time_vec_add/freq.lo)*1000/vec_add_times);
//	//printf("vec_sub time consuming (ms): %f\n",(Time_vec_sub/freq.lo)*1000);
//	//printf("vec_sub solving times: %d\n",vec_sub_times);
//	//printf("vec_sub average solving time(ms): %f\n\n",(Time_vec_sub/freq.lo)*1000/vec_sub_times);
//
//	//printf("vec_div time consuming (ms): %f\n",(Time_vec_div/freq.lo)*1000);
//	//printf("vec_div solving times: %d\n",vec_div_times);
//	//printf("vec_div average solving time(ms): %f\n\n",(Time_vec_div/freq.lo)*1000/vec_div_times);
//	//printf("vec_mul time consuming (ms): %f\n",(Time_vec_mul/freq.lo)*1000);
//	//printf("vec_mul solving times: %d\n",vec_mul_times);
//	//printf("vec_mul average solving time(ms): %f\n\n",(Time_vec_mul/freq.lo)*1000/vec_mul_times);
//
//	//printf("bsxfun time consuming (ms): %f\n",(Time_bsxfun/freq.lo)*1000);
//	//printf("vec_div solving times: %d\n",bsxfun_times);
//	//printf("vec_div average solving time(ms): %f\n\n",(Time_bsxfun/freq.lo)*1000/bsxfun_times);
//
//	//printf("Mat_mul time consuming (ms): %f\n",(Time_Mat_mul/freq.lo)*1000);
//	//printf("Mat_mul solving times: %d\n",mat_mul_times);
//	//printf("Mat_mul average solving time(ms): %f\n\n",(Time_Mat_mul/freq.lo)*1000/mat_mul_times);
//
//	//printf("sca2vec time consuming (ms): %f\n",(Time_sca2vec/freq.lo)*1000);
//	//printf("sca2vec solving times: %d\n",sca2vec_times);
//	//printf("sca2vec average solving time(ms): %f\n\n",(Time_sca2vec/freq.lo)*1000/sca2vec_times);
//
//	//printf("vec_rev time consuming (ms): %f\n",(Time_vec_rev/freq.lo)*1000);
//	//printf("vec_rev solving times: %d\n",vec_rev_times);
//	//printf("vec_rev average solving time(ms): %f\n\n",(Time_vec_rev/freq.lo)*1000/vec_rev_times);
//
//	//printf("sca_vec_mut time consuming (ms): %f\n",(Time_sca_vec_mut/freq.lo)*1000);
//	//printf("sca_vec_mut solving times: %d\n",sca_vec_mut_times);
//	//printf("sca_vec_mut average solving time(ms): %f\n\n",(Time_sca_vec_mut/freq.lo)*1000/sca_vec_mut_times);
//
//
//	// Memory free
//	// main function memory free
//	free(L);free(F);free(Phi);free(B);free(OMEGA_L);free(N_OMEGA_L);free(N_OMEGA_L_T);
//	free(G);free(GL);free(GL_T);free(delta_u_ini);free(y_ini);free(lambda_ini);
//
//	// MPC function memory free
//	free(aug_u_k_1);free(omega_r);free(n_omega_r);free(c);free(c_rev);//free(h_x);
//	free(MPC_temp1);free(MPC_temp2);free(MPC_temp3);free(MPC_temp4);
//
//	// QP function memory free
//	free(x_old);free(rd);free(rp);free(rp_y);
//	//free(y_diag);free(inv_y_diag);free(lambda_diag);free(inv_lambda_diag);
//	free(inv_y_lambda);free(y_lambda_A);
//	free(equation_F);free(equation_b);//free(delta_x_aff);free(delta_y_aff);free(delta_lambda_aff);free(diag);
//	free(delta_x);free(delta_y);free(delta_lambda);free(temp_cond);
//	free(x_try);free(sig_mu);free(center_part);
//	free(QP_temp1);free(QP_temp2);free(QP_temp3);
//
//	// feedback function memory free
//	free(fb_temp1);free(fb_temp2);
//
//	// mcf function memory free
//	free(mcf_C);free(mcf_temp_vec);
//
//	// alpha_decreas function memory free
//	free(alpha_decreas_neg_delta_value);free(alpha_decreas_neg_vec_value);
//	free(alpha_decreas_temp1);free(alpha_decreas_temp2);
//
//	// line_solve function memory free
//	free(line_solve_L);
//	free(line_solve_p);
//
//    printf("exit taskFxn()\n");
//}
//
///*
// *  ======== main ========
// */
//Void main()
//{
////    Task_Handle task;
//    Error_Block eb;
//
//    System_printf("enter main()\n");
//
//    Error_init(&eb);
////    task = Task_create(taskFxn, NULL, &eb);
////    if (task == NULL) {
////        System_printf("Task_create() failed!\n");
////        BIOS_exit(0);
////    }
//
//    BIOS_start();     /* enable interrupts and start SYS/BIOS */
//}
