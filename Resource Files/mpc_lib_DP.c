/*
 * mpc_lib.c
 *
 *  Created on: 2014-3-29
 *      Author: Yi
 *  The c file for functions used on mpc project.
 *	2014.4.10 Test successful!
 *	2014.4.12 Modify the pointer relevent function.
 *  2014.4.12 All functions need for online MPC computing are programmed and tested. First stage accomplished.
 *  Version 2.0 Many function modified!
 *  
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <xdc/runtime/Timestamp.h>

//#include "Algorithm_test\PDIPM_V1.h"
//#include "ti\mathlib\src\sqrtsp\sqrtsp.h"

extern int nx,nu,ny,nyc,P,M,mc,ndec;
extern double R;
extern double Q[];
extern double delta_U_p,delta_U_n,U_p,U_n,Y_p,Y_n;
extern double A_e[],B_e[],C_e[],C_e1[];

extern double xm[],xm_old[];

extern double *mcf_C,*mcf_temp_vec,*alpha_decreas_neg_delta_value, *alpha_decreas_neg_vec_value,*alpha_decreas_temp1,*alpha_decreas_temp2;
extern double *line_solve_L,*line_solve_p;

extern double Time_MPC,Time_QP,Time_QP_step1,Time_QP_step2,Time_QP_step3,Time_QP_step4;
extern double Time_LineSolve1,Time_LineSolve2;
extern double Time_Section1,Time_Section2,Time_Section3,Time_Section4,Time_Section;

extern double QP_Iteration;

extern int Mat_vec_times,vec_add_times,vec_sub_times,vec_div_times,vec_mul_times,bsxfun_times,mat_mul_times,sca2vec_times,vec_rev_times,sca_vec_mut_times,dot_product_times;

extern double Time_Mat_mul,Time_vec_add,Time_vec_sub,Time_vec_div,Time_vec_mul,Time_bsxfun,Time_sca2vec,Time_vec_rev,Time_sca_vec_mut,Time_dot_product;

extern double  *MPC_temp1, *MPC_temp2, *MPC_temp3;
extern double *rd, *rp, *rp_y, *inv_y_lambda,*y_lambda_A, *equation_F, *equation_b;

extern double *fb_temp1,*fb_temp2;

extern double *Q_temp;

extern double *CDC_temp1, *CDC_temp2;

// Print the matrix to the consolo. Matrix stored in a array. w column, h rows.
// Rosetta.
void show_matrix(double *m, int w, int h)
{
	int i, j;
	for (i = 0; i < h; i++) {
		for (j = 0; j < w; j++)
			printf("%2g ", m[i * w + j]);
		putchar('\n');
	}
}


//Find the absloute value of a double parameter.
//2014.4.10 Yi Ding.
double scalar_abs(double x)
{
	if(x<0)
		return -x;
	return x;
}


//Find the minimun one of two scalars.
//2014.4.9 Yi Ding.
double scalar_min(double x,double y)
{
    double z;
    z=(x<y)?x:y;
    return z;
}

//Find the maximun one of two scalars.
//2014.4.11 Yi Ding.
double scalar_max(double x,double y)
{
    double z;
    z=(x>y)?x:y;
    return z;
}

//Filling a vector with a specific scalar.
//double *sca2vec(double sca, int n)
void sca2vec(double sca, int n,double *vec)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	for(i=0;i<n;i++)
		vec[i] = sca;

	//end_time = Timestamp_get32();
	//sca2vec_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_sca2vec+=(double)total_time;

}

//Find the minimum one of a vector.
////2014.4.9 Yi Ding.
double vec_min(double *x, int n)
{
	double min;
	int i;
	min = x[0];
	for(i=1;i<n;i++)
	{
		if(x[i]<min)
			min = x[i];
	}
	return min;
}

//Find the maximum one of a vector.
////2014.4.10 Yi Ding.
double vec_max(double *x, int n)
{
	double max;
	int i;
	max = x[0];
	for(i=1;i<n;i++)
	{
		if(x[i]>max)
			max = x[i];
	}
	return max;
}

//Vector augment!
void vec_aug(double *vec,int m, int n,double *aug)
{
	int i,j;
	//m: vector element?
	//n: how many vector?

	for(i=0;i<n;i++)
	{
		for(j=0;j<m;j++)
		{
			//printf("Q is:\n");show_matrix(vec,1,ny);
			//printf("i is: %d\n j is:%d\n",i,j);
			//printf("Q[] is:%f\n",vec[j]);
			//printf("temp[] is:%f\n",aug[i*m+j]);
			aug[i*m+j] = vec[j];
		}
	}
}


//Vector addation
//2014.4.9 Yi Ding.
//void vec_add(double *a, double *b, int n,double *c)
void vec_add(double * restrict a, double * restrict b, int n,double *c)
{
	int i;
	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	#pragma MUST_ITERATE (3);
	for(i=0;i<n;i++)
	{
		c[i]=a[i]+b[i];
	}

	//end_time = Timestamp_get32();
	//vec_add_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_vec_add+=(double)total_time;
}

//Vector subtraction
//2014.4.9 Yi Ding.
//void vec_sub(double *a, double *b, int n,double *c)
void vec_sub(double * restrict a, double * restrict b, int n,double * restrict c)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	#pragma MUST_ITERATE (3);
	for(i=0;i<n;i++)
	{
		c[i]=a[i]-b[i];
	}

	//end_time = Timestamp_get32();
	//vec_sub_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_vec_sub+=(double)total_time;

	//return c;
}

// Return the product of a scalar and vector;
//void sca_vec_mut(double sca, double *vec, double n, double *pro)
void sca_vec_mut(double sca, double * restrict vec, double n, double * restrict pro)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	//#pragma MUST_ITERATE (3);
	for(i=0;i<n;i++)
	{
		pro[i] = sca*vec[i];
	}

	//end_time = Timestamp_get32();
	//sca_vec_mut_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_sca_vec_mut+=(double)total_time;

}

//This is a function to make a scalar[a] into a vecotr and substract another vector vec_b to get vec_x
//Yi Ding. 2014.5.21
void sca_vec_sub(double a, double *vec_b, int length, double *vec_x)
{
	int i;

	for(i=0;i<length;i++)
		vec_x[i] = a - vec_b[i];

}

//This is a function to do .* operation to 2 vectors.
//Not dot product, but the corresponding element to do multiplication and make a new vector.
//Yi Ding. 2014.5.21
//void vec_mul(double *a, double *b, int m, double *x)
void vec_mul(double * restrict a, double * restrict b, int m, double * restrict x)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	#pragma MUST_ITERATE (3);
	for(i=0;i<m;i++)
		x[i] = a[i]*b[i];

	//end_time = Timestamp_get32();
	//vec_mul_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_vec_mul+=(double)total_time;

}

//This is a function to do ./ operation to 2 vectors.
//Let the corresponding element to do division and make a new vector.
//Yi Ding. 2014.5.21
//void vec_div(double *a, double *b, int m, double *x)
void vec_div(double * restrict a, double * restrict b, int m, double * restrict x)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	#pragma MUST_ITERATE (1);
	for(i=0;i<m;i++)
		x[i] = a[i]/b[i];

	//end_time = Timestamp_get32();
	//vec_div_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_vec_div+=(double)total_time;

}

//This function implement the m-fun bsxfun() is Matlab, but only implement the mat-vec times function.
//Not the same usage with that in Matlab, be careful to use.
//Yi Ding. 2014.5.21
//void bsxfun(double *A, double *x, int row, int col, double *B)
void bsxfun(double * restrict A, double * restrict x, int row, int col, double * restrict B)
{
	int i,j;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	#pragma MUST_ITERATE (12);
	for(i=0;i<row;i++)
		#pragma MUST_ITERATE (3);
		for(j=0;j<col;j++)
			B[i*col+j] = A[i*col+j]*x[i];

	//end_time = Timestamp_get32();
	//bsxfun_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_bsxfun += (double)total_time;
}

//Calculate -a of a
//2014.4.9 Yi Ding.
//void vec_rev(double *a, int n, double *b)
void vec_rev(double * restrict a, int n, double * restrict b)
{
	int i;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//start_time = Timestamp_get32();

	#pragma MUST_ITERATE (1);
	for(i=0;i<n;i++)
	{
		b[i] = -a[i];
	}

	//end_time = Timestamp_get32();
	//vec_rev_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_vec_rev += (double)total_time;
}

//Calculate |a| of a
//2014.5.14 Yi Ding.
//void vec_abs(double * a, int n, double *b)
void vec_abs(double * restrict a, int n, double *restrict b)
{
	int i;
	#pragma MUST_ITERATE (1);
	for(i=0;i<n;i++)
	{
		if(a[i]<0)
			b[i] = -a[i];
		else
			b[i] = a[i];
	}
}

//Make the vector into a diagonal form.
//2014.4.9 Yi Ding.
//Modified @2014.4.16
void vec_diag(double *a, int n, double *diag)
{
	int i;
	int j=0;

	//Here we should clear the diag first! @2014.4.16
	//printf("a is:\n");show_matrix(a,1,n);
	for(i=0;i<n*n;i++)
		diag[i]=0;

	//Insert the diagonal value.
	for(i=0;i<n*n;i+=n)
	{
		diag[i+j]=a[j];
		//printf("a[%d] is %f:\n",j,a[j]);
		//printf("d[%d+%d]] is %f:\n",i,j,diag[i+j]);
		//printf("diag is:\n");show_matrix(diag,n,n);
		j++;
	}
}

// Extract the diagonal elements of a diagonal matrix and form a vector.
void diag2vec(double *diag, int n , double *vec)
{
	int i;
	int j=0;
	for(i=0;i<n*n;i+=n)
	{
		vec[j] = diag[i+j];
		j++;
	}
}

//Calculate the inverse of a diagonal matrix.
//2014.4.9 Yi Ding.
void inv_diag(double *diag, int n ,double *inv_diag)
{
	int i;
	int j=0;
	//Here's a trick to access diagnose element using a single cycle but not a double cycle.
	for(i=0;i<n*n;i+=n)
	{
		inv_diag[i+j]=1/diag[i+j];
		j++;
	}
}


//Add two matrices
//2014.4.9 Yi Ding.
void mat_add(double *a, double *b, int row, int column, double *sum)
{
	int i;
	for(i=0;i<row*column;i++)
	{
		sum[i] = a[i]+b[i];
	}
}

//Do the multiplication A*x = b;
//Yi Ding. 2014.4.25
//void mat_vec(double * A, double * x, int row, int col, double *b)
void mat_vec(double * restrict A, double * restrict x, int row, int col, double *restrict b)
{
	int i,j;

	#pragma MUST_ITERATE (3);
	for(i=0;i<row;i++)
	{
		b[i]=0;
		#pragma MUST_ITERATE (3);
		for(j=0;j<col;j++)
		{
				b[i]+=A[i*col+j]*x[j];				//Attention!2014.5.5
		}
	}

}

//Do the multiplication A*x = b;
//Yi Ding. 2014.4.25
//void mat_vec(double * A, double * x, int row, int col, double *b)
void mat_vec_offline(double * restrict A, double * restrict x, int row, int col, double *restrict b)
{
	int i,j;

	for(i=0;i<row;i++)
	{
		b[i]=0;
		for(j=0;j<col;j++)
		{
				b[i]+=A[i*col+j]*x[j];				//Attention!2014.5.5
		}
	}

}


//Do the multiplication x'*A = b;
//Yi Ding. 2014.7.23
void vec_mat(double *A, double *x, int row, int col, double *b)
{
	int i,j;
	#pragma MUST_ITERATE (3);
	for(j=0;j<col;j++)
	{
		b[j]=0;
		#pragma MUST_ITERATE (3);
		for(i=0;i<row;i++)
		{
				b[j]+=A[i*col+j]*x[i];				//Attention!2014.5.5
		}
	}
}

//Do the multiplication A*diag(x) = b;
//Yi Ding. 2014
void mat_diag(double *A, double *diag, int row, int col, double *b)
{
	int i,j;

	for(i=0;i<row;i++)
	{
		for(j=0;j<col;j++)
			b[i*col+j] = A[i*col+j]*diag[j*col+j];
	}
}

//Used for LU Evaluate!
//2014.5.15
void Lower_diag_vec(double *A, double *diag_vec, int m, double *b)
{
	int i,j;
	double temp;
	for(j=0;j<m;j++)
	{
		temp = diag_vec[j];
		for(i=j;i<m;i++)
		{
			b[i*m+j] = A[i*m+j]*temp;
			//show_matrix(b,m,m);
		}
	}
}

//Do the A*X = B, here A, X, B are all diagnal matrix
//Yi Ding 2014.4.26
void diag_diag(double *A, double *X, int n, double *B)
{
	int i;

	for(i=0;i<n;i++)
		B[i*n+i] = A[i*n+i]*X[i*n+i];
}

//Do the A*x = b, here A is diag matrix, x is a vector
//Yi Ding 2014.4.26
void diag_vec(double *A, double *x, int n, double *b)
{
	int i;
	for(i=0;i<n;i++)
		b[i] = A[i*n+i]*x[i];
}

//Matrix multiplication. Just like DSP_sp_mat_mul.c in ti DSPlib.
//2014.4.8 Yi Ding.
void mat_mul(double * restrict x1, int r1, int c1, double * restrict x2, int c2, double * restrict y)
{
	int i,j,k;
	int length;

	//Uint32 start_time,end_time,total_time,dead_time;
	//start_time = Timestamp_get32();
	//end_time = Timestamp_get32();
	//dead_time = end_time - start_time;

	//length = r1*c2;

	//for(i=0;i<length;i++)
	//{
	//	y[i]=0;
	//}

	#pragma MUST_ITERATE (3);
	for(i=0;i<r1;i++)
	{
	#pragma MUST_ITERATE (3);
		for(j=0;j<c2;j++)
		{
			y[(i)*c2+(j)]=0;
			#pragma MUST_ITERATE (12);
			for(k=0;k<c1;k++)
			{
				y[(i)*c2+(j)]+=x1[(i)*c1+(k)]*x2[(k)*c2+(j)];
			}
		}
	}
	//end_time = Timestamp_get32();
	//mat_mul_times++;
	//total_time = end_time - start_time - dead_time;
	//Time_Mat_mul+=(double)total_time;
}

//Matrix multiplication. For off line use. For the accurate running time of on line mat_mul.
//2014.7.18 Yi Ding.
void mat_mul_offline(double *x1, int r1, int c1, double *x2, int c2, double *y)
{
	int i,j,k;

	for(i=0;i<r1;i++)
	{
		for(j=0;j<c2;j++)
		{
			y[(i)*c2+(j)]=0;
			for(k=0;k<c1;k++)
			{
				y[(i)*c2+(j)]+=x1[(i)*c1+(k)]*x2[(k)*c2+(j)];
			}
		}
	}
}


//Compute the power of a matrix
//2014.4.16
void mat_pow(double *mat, int row, int pow, double *MatPow)
{
	int i;
	double *temp1;

	temp1 = (double*)calloc(row*row,sizeof(double));
	//temp2 = (double*)calloc(row*row,sizeof(double));

	if(pow==0)
	{
		sca2vec(1,row,temp1);
		vec_diag(temp1,row,MatPow);
	}
	else
	{
		memcpy(temp1,mat,row*row*sizeof(double));
		memcpy(MatPow,mat,row*row*sizeof(double));
		for(i=1;i<pow;i++)
		{
			mat_mul_offline(temp1,row,row,mat,row,MatPow);
			memcpy(temp1,MatPow,row*row*sizeof(double));
		}
	}
	free(temp1);
}

// Cholesky decomposition.
// 2014.3.29 Yi Ding. Rosetta.
void cholesky(double *A, int n, double *L)
{
	int i,j,k;
	double s;
	//double *L = (double*)calloc(n * n, sizeof(double));
 //   if (L == NULL)
 //       exit(EXIT_FAILURE);

    for (i = 0; i < n; i++)
        for (j = 0; j < (i+1); j++) {
        	s = 0;
            for (k = 0; k < j; k++)
                s += L[i * n + k] * L[j * n + k];
				L[i * n + j] = (i == j) ?
            			sqrt(A[i * n + i] - s) :
						//sqrtsp(A[i * n + i] - s) :
                           (1.0 / L[j * n + j] * (A[i * n + j] - s));
        }
}

//A Cholesky Factorization routine from the book Numerical Recipes in C
// 2014.7.15 Yi Ding
//void chol_NRC(double *A, int n,double *p)
int chol_NRC(double * restrict A, int n,double * restrict p)
{
	int i,j,k;
	int flag = 1;
	double sum;
	#pragma MUST_ITERATE (3);
	for (i=1;i<=n;i++)
	{
		#pragma MUST_ITERATE (3);
		for (j=i;j<=n;j++)
		{
			for (sum=A[(i-1)*n+(j-1)],k=i-1;k>=1;k--)
			{
				sum -= A[(i-1)*n+(k-1)]*A[(j-1)*n+(k-1)];
				//printf("sum: %f\n",sum);
				//show_matrix(A,n,n);putchar('\n');
				//printf("A[%d][%d]: %f\n",i,k,A[(i-1)*n+(k-1)]);
				//printf("A[%d][%d]: %f\n",j,k,A[(j-1)*n+(k-1)]);
			}
			//show_matrix(A,n,n);putchar('\n');
			//printf("A[%d][%d]: %f\n",i,j,A[(i-1)*n+(j-1)]);
			if (i == j)
			{
				//printf("sum: %f\n",sum);
				if (sum <= 0.0)
				{
					printf("Factorization Failed.\n");
					flag = 0;
				}
				p[i-1]=sqrt(sum);
				//printf("p[%d]: %f\n",i,p[i-1]);
			}
			else A[(j-1)*n+(i-1)]=sum/p[i-1];
			//show_matrix(A,n,n);putchar('\n');
		}
	}
	return flag;
}

//Here, the diagnal D is stored in a vector.
//Modified Cholesky Factorization 2014.5.15
void mcf(double *A, int m, double *L, double *D)
{
	int i,j,k;

	double beta = 100;
	double sigma = 0.1;
	double temp;


	for(i=0;i<m;i++)
	{
		L[i] = 0;
		D[i] = 0;
	}
	for(i=0;i<m;i++)
		L[i*m+i] = 1;

	mcf_C[0] = A[0];
	D[0] = scalar_max(scalar_max(scalar_abs(mcf_C[0]),(mcf_C[0]/beta)*(mcf_C[0]/beta)),sigma);
	//printf("D[0] is: %f\n",D[0]);

	for(k=1;k<m;k++)
	{
		L[k*m] = A[k*m]/D[0];
	}
	//show_matrix(L,m,m);putchar('\n');

	for(j=1;j<m;j++)
	{
		temp = 0;
		for(k=0;k<i;k++)
			temp += D[k]*L[j*m+k]*L[j*m+k];
		mcf_C[j*m+j] = A[j*m+j]-temp;
		//printf("C is:\n");
		//show_matrix(C,m,m);putchar('\n');

		for(k=j-1;k<m;k++)
			mcf_temp_vec[k-j+1] = mcf_C[k*m+j];
		vec_abs(mcf_temp_vec,m-j+1,mcf_temp_vec);
		//show_matrix(temp_vec,1,m-j+1);putchar('\n');
		temp = vec_max(mcf_temp_vec,m-j+1);
		temp = temp*temp;
		D[j] = scalar_max(scalar_max(scalar_abs(mcf_C[j*m+j]),temp),sigma);
		//printf("D is:\n");
		//show_matrix(D,1,m);putchar('\n');

		for(i=j+1;i<m;i++)
		{
			temp = 0;
			for(k=0;k<j;k++)
			{
				temp += D[k]*L[i*m+k]*L[j*m+k];
				//printf("L is:\n");
				//show_matrix(L,m,m);
				//printf("i is:%d\n",i);
				//printf("j is:%d\n",j);
				//printf("k is:%d\n",k);
				//printf("D[k] is:%f\n",D[k]);
				//printf("L[i*m+k] is:%f\n",L[i*m+k]);
				//printf("L[j*m+k] is:%f\n",L[j*m+k]);
				//printf("temp is:%f\n",temp);
				//putchar('\n');putchar('\n');
			}

			mcf_C[i*m+j] = A[i*m+j]-temp;
			//printf("C is:\n");
			//show_matrix(C,m,m);putchar('\n');

			L[i*m+j] = mcf_C[i*m+j]/D[j];
			//printf("L is:\n");
			//show_matrix(L,m,m);putchar('\n');
		}
	}
}


// Do the dot product for two array. n is the number of elements
// 2014.3.29 Yi Ding. Rosetta.
//double dot_product(double * a, double * b, int n)
double dot_product(double * restrict a, double * restrict b, int n)
{
        double sum = 0;
        int i;

    	//Uint32 start_time,end_time,total_time,dead_time;
        //start_time = Timestamp_get32();
        //end_time = Timestamp_get32();
        //dead_time = end_time - start_time;

    	//start_time = Timestamp_get32();

		#pragma MUST_ITERATE (3);
        for (i = 0; i < n; i++) {
                sum += a[i] * b[i];
        }

        //end_time = Timestamp_get32();
        //dot_product_times++;
        //total_time = end_time - start_time - dead_time;
        //Time_dot_product += (double)total_time;

        return sum;

}

//This function implement a special searching method for alpha.
//The idea comes from the quad_wright() and symbol meanings see algorithms.
// Yi Ding. 2014.5.21
double alpha_decreas(double *vec, double *delta_vec,int length, double tau)
{
	int i;
	int neg_num = 0;
	//double *neg_index;
	double alpha = 0.9999995;

	//neg_index = (double*)calloc(length,sizeof(double));

	#pragma MUST_ITERATE (12);
	for(i=0;i<length;i++)
	{
		if(delta_vec[i]<0)
		{
			//neg_index[neg_num] = i;
			alpha_decreas_neg_delta_value[neg_num] = delta_vec[i];
			alpha_decreas_neg_vec_value[neg_num] = vec[i];
			neg_num = neg_num + 1;
		}
	}

	if(neg_num>0)
	{
		if(tau==1)
			vec_rev(alpha_decreas_neg_vec_value,neg_num,alpha_decreas_temp1);
		else
			sca_vec_mut(-tau,alpha_decreas_neg_vec_value,neg_num,alpha_decreas_temp1);
		vec_div(alpha_decreas_temp1,alpha_decreas_neg_delta_value,neg_num,alpha_decreas_temp2);
		alpha = 0.9995*vec_min(alpha_decreas_temp2,neg_num);
	}

	return alpha;

}

// Alpha calculation method (while loop)
double alpha_decreas_while(double *y, double *lambda, double *delta_y, double *delta_lambda, int m)
{

	double alpha_cond;	//Use to check whether the alpha is valid.
	double alpha = 1;
	double *temp1,*temp2,*temp3;

	temp1 = (double*)calloc(m,sizeof(double));
	temp2 = (double*)calloc(m,sizeof(double));
	temp3 = (double*)calloc(m,sizeof(double));

	sca_vec_mut(alpha,delta_y,m,temp1);
	vec_add(y,temp1,m,temp2);
	sca_vec_mut(alpha,delta_lambda,m,temp1);
	vec_add(lambda,temp1,m,temp3);
	alpha_cond = scalar_min(vec_min(temp2,m),vec_min(temp3,m));
	while(alpha_cond<0)
	{
		alpha = alpha - 0.01;		//Attention!
		//alpha_aff = alpha_aff/2;
		if(alpha<=0)
			break;
		sca_vec_mut(alpha,delta_y,m,temp1);
		vec_add(y,temp1,m,temp2);
		sca_vec_mut(alpha,delta_lambda,m,temp1);
		vec_add(lambda,temp1,m,temp3);
		alpha_cond = scalar_min(vec_min(temp2,m),vec_min(temp3,m));
	}

	free(temp1);free(temp2);free(temp3);

	return alpha;
}
// Test accomplished! 2014.4.9


// Do back substatution to solve linear equation using LU.
// 2013.3.30 Yi Ding.
//Modified @2014.4.18
void luEvaluate(double *L,double *U, double*b,int n,double *x)
{
	double *y = (double*)calloc(n,sizeof(double));
	int i,j;
	double temp = 0;
	if(x == NULL || y == NULL)
		exit(0);

	//Foward solve Ly = b;
	y[0] = b[0]/L[0];

	#pragma MUST_ITERATE (3);
	for(i=1;i<n;i++)
	{
		for(j=0;j<i;j++)
		{
			temp += L[i*n+j]*y[j];
		}
		y[i] = b[i] - temp;
		y[i] = y[i]/L[i*n+i];
		temp = 0;
	}
	//show_matrix(y,1,n);

	//Backward solve Ux = y
	x[n-1] = y[n-1]/U[n*n-1];
	temp = 0;
	#pragma MUST_ITERATE (1);
	for(i=n-2;i>=0;i--)
	{
		for(j=i+1;j<n;j++)
		{
			temp += U[i*n+j]*x[j];
		}
		x[i] = y[i] - temp;
		x[i] = x[i]/U[i*n+i];
		temp = 0;
	}
	free(y);
}


// Do back substatution to solve linear equation using LU.
// A routine from the book Numerical Recipes in C
// Yi Ding 2014.7.15
//void luEvaluate_NRC(double *L, double*b,double *p,int n,double *x)
void luEvaluate_NRC(double * restrict L, double* restrict b,double * restrict p,int n,double * restrict x)
{
int i,k;
double sum;

	#pragma MUST_ITERATE (3);
	for (i=1;i<=n;i++)
	{
		//#pragma MUST_ITERATE (3);
		for (sum=b[i-1],k=i-1;k>=1;k--)
			sum -= L[(i-1)*n+(k-1)]*x[k-1];
		x[i-1]=sum/p[i-1];
	}

	#pragma MUST_ITERATE (3);
	for (i=n;i>=1;i--)
	{
		//#pragma MUST_ITERATE (3);
		for (sum=x[i-1],k=i+1;k<=n;k++)
			sum -= L[(k-1)*n+(i-1)]*x[k-1];
		x[i-1]=sum/p[i-1];
	}
}

//Transpose the matrix m and store the matrix in w.
//2014.4.8 Yi Ding.
void transpose(double *m, int w, int h)
{
	int start, next, i;
	double tmp;

	for (start = 0; start <= w * h - 1; start++) {
		next = start;
		i = 0;
		do {	i++;
			next = (next % h) * w + next / h;
		} while (next > start);
		if (next < start || i == 1) continue;

		tmp = m[next = start];
		do {
			i = (next % h) * w + next / h;
			m[next] = (i == start) ? tmp : m[i];
			next = i;
		} while (next > start);
	}
}

//Solve linear equation using cholesky decomposition.
//Simply calling two other functions: cf() and luEvaluate()
//2014.3.30 Yi Ding.
//double *line_solve(double *A, double *b, int n)
void line_solve(double *A, double *b, int n, double *x)
{
	//double *U = (double*)calloc(n * n, sizeof(double));
	//double *D = (double*)calloc(n, sizeof(double));
	//double *LD = (double*)calloc(n * n, sizeof(double));


	////Below is Standard Cholesky Factorization
	//cholesky(A,n,L);
	//memcpy(U,L,n*n*sizeof(double));
	//transpose(U, n, n);
	//luEvaluate(L,U,b,n,x);

	//Below is the Linear solve process from the book Numerical Recipes in C
	memcpy(line_solve_L,A,n*n*sizeof(double));
	chol_NRC(line_solve_L,n,line_solve_p);
	luEvaluate_NRC(line_solve_L,b,line_solve_p,n,x);

	////Below is Modified Cholesky Factorization
	//mcf(A,n,L,D);
	//memcpy(U,L,n*n*sizeof(double));
	//transpose(U, n, n);
	//Lower_diag_vec(L,D,n,LD);
	//luEvaluate(LD,U,b,n,x);


	//free(U);
	//free(D);
	//free(LD);

}


//The feedback step in MPC, a transformation of funciton fankui.m
void feedback_v3(double *x_k_1, double *y_k, double *u_k_1, double *u_k_2, double *x_k,double *L)
{
	int i;

	//vec_sub(u_k_1,u_k_2,nu,fb_temp1);
	for(i=0;i<nu;i++)
	{
		fb_temp1[i] = u_k_1[i] - u_k_2[i];
	}
	mat_vec_offline(B_e,fb_temp1,nx+ny,nu,fb_temp2);
	mat_vec_offline(A_e,x_k_1,nx+ny,nx+ny,fb_temp1);
	vec_add(fb_temp1,fb_temp2,nx+ny,x_k);

	mat_vec_offline(C_e,x_k,ny,nx+ny,fb_temp1);
	//vec_sub(y_k,fb_temp1,ny,fb_temp2);
	for(i=0;i<ny;i++)
	{
		fb_temp2[i] = y_k[i] - fb_temp1[i];
	}
	mat_vec_offline(L,fb_temp2,nx+ny,ny,fb_temp1);
	vec_add(x_k,fb_temp1,nx+ny,fb_temp2);

	memcpy(x_k,fb_temp2,(nx+ny)*sizeof(double));

}
//Feedbcak test accomplished @2014.4.14

//The feedback step in MPC, a transformation of funciton fankui.m
//Use Kalman filter to do states estimation.
void feedback_kalman(double *x_k_1, double *y_k, double *u_k_1, double *u_k_2, double *x_k,int iter)
{

	int i;

	//if(iter>=99)
	//	iter=99;

	double K[12] = {0.4624,-0.0189,0.4175,-0.3376,-0.0056,0.6674,-0.0489,3.8177,1.2928,-0.0209,-0.0076,1.5291};

	//vec_sub(u_k_1,u_k_2,nu,fb_temp1);
	for(i=0;i<nu;i++)
	{
		fb_temp1[i] = u_k_1[i] - u_k_2[i];
	}
	mat_vec_offline(B_e,fb_temp1,nx+ny,nu,fb_temp2);
	mat_vec_offline(A_e,x_k_1,nx+ny,nx+ny,fb_temp1);
	vec_add(fb_temp1,fb_temp2,nx+ny,x_k);

	//printf("K[%d] is :\n",iter);show_matrix(Kalman_K[iter],2,6);

	mat_vec_offline(C_e,x_k,ny,nx+ny,fb_temp1);
	//vec_sub(y_k,fb_temp1,ny,fb_temp2);
	for(i=0;i<ny;i++)
	{
		fb_temp2[i] = y_k[i] - fb_temp1[i];
	}
	//mat_vec_offline(Kalman_K[iter],fb_temp2,nx+ny,ny,fb_temp1);
	mat_vec_offline(K,fb_temp2,nx+ny,ny,fb_temp1);
	vec_add(x_k,fb_temp1,nx+ny,fb_temp2);

	memcpy(x_k,fb_temp2,(nx+ny)*sizeof(double));
}

//The formation of omega_r.
//This version has y contraints.
void omega_r_form(double *aug_u_k_1, double *x_k, double *omega_r,double *F1)
{
	int i;
	double *temp1;
	temp1 = (double*)calloc(nyc*P,sizeof(double));
	//mat_vec_offline(F1,x_k,nyc*P,nx+ny,temp1);
	//printf("F1 is:\n");show_matrix(F1,nx+ny,nyc*P);putchar('\n');
	//printf("x_k is:\n");show_matrix(x_k,1,nx+ny);putchar('\n');
	mat_vec(F1,x_k,nyc*P,nx+ny,temp1);
	//printf("temp1 is:\n");show_matrix(temp1,1,nyc*P);putchar('\n');

	for(i=0;i<nu*M;i++)
		omega_r[i] = delta_U_p;
	for(i=nu*M;i<2*nu*M;i++)
		omega_r[i] = -delta_U_n;
	for(i=2*nu*M;i<3*nu*M;i++)
		omega_r[i] = U_p - aug_u_k_1[i-2*nu*M];
	for(i=3*nu*M;i<4*nu*M;i++)
		omega_r[i] = -U_n + aug_u_k_1[i-3*nu*M];
	for(i=4*nu*M;i<4*nu*M+nyc*P;i++)
		omega_r[i] = Y_p - temp1[i-4*nu*M];
	for(i=4*nu*M+nyc*P;i<4*nu*M+2*nyc*P;i++)
		omega_r[i] = -Y_n + temp1[i-4*nu*M-nyc*P];
	free(temp1);
}

//The formation of omega_r.
//This version does not have y contraints.
void omega_r_form2(double *aug_u_k_1, double *x_k, double *omega_r, double *F)
{
	int i;
	//double *temp1;
	//temp1 = (double*)calloc(ny*P,sizeof(double));
	////mat_mul(F,ny*P,nx+ny,x_k,1,temp1);
	//mat_vec(F,x_k,ny*P,nx+ny,temp1);

	for(i=0;i<nu*M;i++)
		omega_r[i] = delta_U_p;
	for(i=nu*M;i<2*nu*M;i++)
		omega_r[i] = -delta_U_n;
	for(i=2*nu*M;i<3*nu*M;i++)
		omega_r[i] = U_p - aug_u_k_1[i-2*nu*M];
	for(i=3*nu*M;i<4*nu*M;i++)
		omega_r[i] = -U_n + aug_u_k_1[i-3*nu*M];
	//for(i=4*nu*M;i<4*nu*M+ny*P;i++)
	//	omega_r[i] = Y_p - temp1[i-4*nu*M];
	//for(i=4*nu*M+ny*P;i<4*nu*M+2*ny*P;i++)
	//	omega_r[i] = -Y_n + temp1[i-4*nu*M-ny*P];
	//free(temp1);
}

//The formation of omega_r.
//This version has y contraints but do not have delta_u constraints
void omega_r_form3(double *aug_u_k_1, double *x_k, double *omega_r, double *F,double Y_p1, double Y_p2, double Y_n1, double Y_n2)
{
	int i;
	double *temp1;
	temp1 = (double*)calloc(nyc*P,sizeof(double));
	//mat_mul(F,ny*P,nx+ny,x_k,1,temp1);
	mat_vec(F,x_k,nyc*P,nx+ny,temp1);

	for(i=0;i<nu*M;i++)
		omega_r[i] = U_p - aug_u_k_1[i];
	for(i=nu*M;i<2*nu*M;i++)
		omega_r[i] = -U_n + aug_u_k_1[i-nu*M];
	for(i=2*nu*M;i<2*nu*M+nyc*P;i=i+2)
	{
		omega_r[i] = Y_p1 - temp1[i-2*nu*M];
		omega_r[i+1] = Y_p2 - temp1[i-2*nu*M+1];
	}
	for(i=2*nu*M+nyc*P;i<2*nu*M+2*nyc*P;i=i+2)
	{
		omega_r[i] = -Y_n1 + temp1[i-2*nu*M-nyc*P];
		omega_r[i+1] = -Y_n2 + temp1[i-2*nu*M-nyc*P+1];
	}
	free(temp1);
}

//Find the parameter F and Phi, not in the MPC loop. Relaxation on performance.
//Yi Ding. 2014.4.16
void fphi(double *F, double *Phi)
{
	int i,j,k,kk;
	double *temp1, *temp2;
	temp1 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));
	temp2 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));

	mat_mul_offline(C_e,ny,nx+ny,A_e,nx+ny,temp2);
	for(i=0;i<P;i++)
	{
		for(j=0;j<ny*(ny+nx);j++)
		{
			F[i*(ny*(ny+nx))+j] = temp2[j];
		}
		memcpy(temp1,temp2,ny*(nx+ny)*sizeof(double));
		mat_mul_offline(temp1,ny,nx+ny,A_e,nx+ny,temp2);
	}

	for(i=0;i<M;i++)
	{
		for(j=0;j<P;j++)
		{
			if(j<i)
			{
				sca2vec(0,nu*ny,temp1);
			}
			else
			{
				mat_pow(A_e,nx+ny,j-i,temp1);
				mat_mul_offline(C_e,ny,nx+ny,temp1,nx+ny,temp2);
				mat_mul_offline(temp2,ny,nx+ny,B_e,nu,temp1);
			}
			for(k=0;k<nu;k++)
			{
				for(kk=0;kk<ny;kk++)
				{
					Phi[j*ny*nu*M+kk*nu*M+i*nu+k] = temp1[kk*nu+k];
				}
			}
		}
	}
	free(temp1);free(temp2);
}
//fphi() Test accomplished @2014.4.16
//Attention! Find a bug in vec_diag(), modified!

//Find the parameter F and Phi, not in the MPC loop. Relaxation on performance.
//Yi Ding. 2014.4.16
void fphi1(double *F1, double *Phi1)
{
	int i,j,k,kk;
	double *temp1, *temp2;
	int ny1 = 1;
	temp1 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));
	temp2 = (double*)calloc((nx+ny)*(nx+ny),sizeof(double));

	mat_mul(C_e1,1,nx+ny,A_e,nx+ny,temp2);
	for(i=0;i<P;i++)
	{
		for(j=0;j<(ny+nx);j++)
		{
			F1[i*(ny+nx)+j] = temp2[j];
		}
		memcpy(temp1,temp2,(nx+ny)*sizeof(double));
		mat_mul(temp1,1,nx+ny,A_e,nx+ny,temp2);
	}

	for(i=0;i<M;i++)
	{
		for(j=0;j<P;j++)
		{
			if(j<i)
			{
				sca2vec(0,nu*ny1,temp1);
			}
			else
			{
				mat_pow(A_e,nx+ny,j-i,temp1);
				mat_mul(C_e1,ny1,nx+ny,temp1,nx+ny,temp2);
				mat_mul(temp2,ny1,nx+ny,B_e,nu,temp1);
			}
			//printf("temp1 is :\n");show_matrix(temp1,ny1,nu);putchar('\n');
			for(k=0;k<nu;k++)
			{
				for(kk=0;kk<ny1;kk++)
				{
					//printf("Phi1 No. %d\n:",j*ny1*nu*M+kk*nu*M+i*nu+k);
					//printf("temp No. %d\n:",kk*nu+k);
					Phi1[j*ny1*nu*M+kk*nu*M+i*nu+k] = temp1[kk*nu+k];
				}
			}
		}
	}
	free(temp1);free(temp2);
}


//Form the B used in the OMEGA_L
void form_B(double *B)
{
	int i,j,k;
	for(i=0;i<nu*M;i++)
	{
		for(j=0;j<nu*M;j++)
		{
			for(k=0;k<M;k++)
				if(i==(j+k*nu))
					B[i*nu*M+j]=1;
		}
	}
}

//Form the OMEGA_L for the constrains in the optimization problems.
//This version has y constraints.
void OMEGA_L_form(double *B, double *Phi1, double *OMEGA_L)
{
	int i,j;

	for(i=0;i<nu*M;i++)
	{
		for(j=0;j<nu*M;j++)
		{
			if(i==j)
				OMEGA_L[i*nu*M+j] = 1;
		}
	}
	for(i=nu*M;i<2*nu*M;i++)
	{
		for(j=0;j<nu*M;j++)
		{
			if((i-nu*M)==j)
				OMEGA_L[i*nu*M+j] = -1;
		}
	}
	for(i=(2*nu*M)*nu*M;i<(3*nu*M)*nu*M;i++)
		OMEGA_L[i] = B[i-(2*nu*M)*nu*M];
	for(i=(3*nu*M)*nu*M;i<(4*nu*M)*nu*M;i++)
		OMEGA_L[i] = -B[i-(3*nu*M)*nu*M];
	for(i=(4*nu*M)*nu*M;i<(4*nu*M+nyc*P)*nu*M;i++)
		OMEGA_L[i] = Phi1[i-(4*nu*M)*nu*M];
	for(i=(4*nu*M+nyc*P)*nu*M;i<(4*nu*M+2*nyc*P)*nu*M;i++)
		OMEGA_L[i] = -Phi1[i-(4*nu*M+nyc*P)*nu*M];
}

//Form the OMEGA_L for the constrains in the optimization problems.
//This version does not have y constraints.
void OMEGA_L_form2(double *B, double *Phi, double *OMEGA_L)
{
	int i,j;

	for(i=0;i<nu*M;i++)
	{
		for(j=0;j<nu*M;j++)
		{
			if(i==j)
				OMEGA_L[i*nu*M+j] = 1;
		}
	}
	for(i=nu*M;i<2*nu*M;i++)
	{
		for(j=0;j<nu*M;j++)
		{
			if((i-nu*M)==j)
				OMEGA_L[i*nu*M+j] = -1;
		}
	}
	for(i=(2*nu*M)*nu*M;i<(3*nu*M)*nu*M;i++)
		OMEGA_L[i] = B[i-(2*nu*M)*nu*M];
	for(i=(3*nu*M)*nu*M;i<(4*nu*M)*nu*M;i++)
		OMEGA_L[i] = -B[i-(3*nu*M)*nu*M];
	//for(i=(4*nu*M)*nu*M;i<(4*nu*M+ny*P)*nu*M;i++)
	//	OMEGA_L[i] = Phi[i-(4*nu*M)*nu*M];
	//for(i=(4*nu*M+ny*P)*nu*M;i<(4*nu*M+2*ny*P)*nu*M;i++)
	//	OMEGA_L[i] = -Phi[i-(4*nu*M+ny*P)*nu*M];
}

//Form the OMEGA_L for the constrains in the optimization problems.
//This version has y constraints. But DO not have delta_u constraints.
void OMEGA_L_form3(double *B, double *Phi, double *OMEGA_L)
{
	int i,j;

	for(i=0;i<nu*M*nu*M;i++)
		OMEGA_L[i] = B[i];
	for(i=nu*M*nu*M;i<(2*nu*M)*nu*M;i++)
		OMEGA_L[i] = -B[i-(nu*M)*nu*M];
	for(i=(2*nu*M)*nu*M;i<(2*nu*M+nyc*P)*nu*M;i++)
		OMEGA_L[i] = Phi[i-(2*nu*M)*nu*M];
	for(i=(2*nu*M+nyc*P)*nu*M;i<(2*nu*M+2*nyc*P)*nu*M;i++)
		OMEGA_L[i] = -Phi[i-(2*nu*M+nyc*P)*nu*M];
}

//Form the G for the opimization problem.
void form_G(double *Phi,double *G)
{
	double *G_temp1,*G_temp2,*G_temp3;
	int max_mem;

	//printf("Phi is: \n");show_matrix(Phi,nu*M,ny*P);putchar('\n');

	max_mem = scalar_max(P*ny*P*ny,M*nu*P*ny);
	G_temp1 = (double*)calloc(max_mem,sizeof(double));
	G_temp2 = (double*)calloc(max_mem,sizeof(double));
	G_temp3 = (double*)calloc(max_mem,sizeof(double));
	//sca2vec(Q,ny*P,temp1);
	//vec_aug(Q,ny,P,G_temp1);

	//printf("G_temp1 is:\n");show_matrix(G_temp1,1,ny*P);
	//vec_diag(G_temp1,ny*P,G_temp2);
	//printf("G_temp2 is:\n");show_matrix(G_temp2,ny*P,ny*P);
	memcpy(G_temp1,Phi,M*nu*P*ny*sizeof(double));
	transpose(G_temp1, M*nu, P*ny);
	//printf("G_temp1 is:\n");show_matrix(G_temp1,M*nu,ny*P);

	//Problems!Here!
	mat_mul(G_temp1,M*nu,P*ny,G_temp2,P*ny,G_temp3);
	//printf("G_temp3 is:\n");show_matrix(G_temp3,nu*M,P*ny);
	mat_mul(G_temp3,M*nu,P*ny,Phi,M*nu,G_temp1);
	//printf("G_temp1 is:\n");show_matrix(G_temp1,M*nu,M*nu);
	sca2vec(R,nu*M,G_temp2);
	//printf("G_temp2 is:\n");show_matrix(G_temp2,1,M*nu);
	vec_diag(G_temp2,nu*M,G_temp3);
	//printf("G_temp3 is:\n");show_matrix(G_temp3,nu*M,M*nu);
	mat_add(G_temp1,G_temp3,nu*M,nu*M,G);

	//printf("G is: \n");show_matrix(G,ndec,ndec);putchar('\n');

	free(G_temp1);free(G_temp2);free(G_temp3);
}

void form_G2(double *Phi,double *G)
{
	double *temp1,*temp2,*temp3;

	temp1 = (double*)calloc(M*nu*P*ny,sizeof(double));
	temp2 = (double*)calloc(M*nu*P*ny,sizeof(double));
	temp3 = (double*)calloc(M*nu*M*nu,sizeof(double));

	memcpy(temp1,Phi,M*nu*P*ny*sizeof(double));
	transpose(temp1, M*nu, P*ny);
	bsxfun(Phi,Q_temp,P*ny,M*nu,temp2);
	mat_mul(temp1,M*nu,P*ny,temp2,M*nu,temp3);

	sca2vec(R,nu*M,temp1);
	vec_diag(temp1,nu*M,temp2);
	mat_add(temp2,temp3,nu*M,nu*M,G);

	free(temp1);free(temp2);free(temp3);
}

//This is a routine to mimic the model using the disturbed state-space model.
void process_model(double* xm, double *u_k, double *y_k)
{
	double Ac[] = {1,0.1,0,0.9};
	double Bc[] = {0,0.0787};
	double Cc[] = {1,0};

	int i;

	double *temp1,*temp2;
	temp1 = (double*)calloc(nx,sizeof(double));
	temp2 = (double*)calloc(nx,sizeof(double));

	mat_vec_offline(Ac,xm,nx,nx,temp1);
	mat_vec_offline(Bc,u_k,nx,nu,temp2);
	//vec_add(temp1,temp2,nx,xm);
	for(i=0;i<nx;i++)
	{
		xm[i]=temp1[i]+temp2[i];
	}

	mat_vec_offline(Cc,xm,ny,nx,y_k);

	free(temp1);free(temp2);
}

void process_model2(double* xm, double *u_k, double *y_k)
{
	double Ac[] = {1,0.05,-0.0057,-0.000094883,0,1,-0.2308,-0.0057,0,0,1.0593,0.0510,0,0,2.3968,1.0593};
	double Bc[] = {0.0028,0.1106,-0.0091,-0.3674};
	double Cc[] = {1,0,0,0,0,0,1,0};

	int i;

	double *temp1,*temp2;
	temp1 = (double*)calloc(nx,sizeof(double));
	temp2 = (double*)calloc(nx,sizeof(double));

	for(i=0;i<nx;i++)
		xm_old[i] = xm[i];

	mat_mul(Ac,nx,nx,xm,1,temp1);
	mat_mul(Bc,nx,nu,u_k,1,temp2);
	vec_add(temp1,temp2,nx,xm);

	mat_mul(Cc,ny,nx,xm,1,y_k);

	free(temp1);free(temp2);
}

void process_model3(double* xm, double *u_k, double *y_k)
{
	double Ac[] = { 0.24,0,0.1787,0,-0.3722,1,0.2703,0,-0.9901,0,0.1389,0,-48.9354,64.1,2.3992,1};
	double Bc[] = {-1.2346,-1.4383,-4.4828,-1.7999};
	double Cc[] = {0,1,0,0,0,0,0,1,-128.2,128.2,0,0};

	double *temp1,*temp2;
	temp1 = (double*)calloc(nx,sizeof(double));
	temp2 = (double*)calloc(nx,sizeof(double));

	mat_mul(Ac,nx,nx,xm,1,temp1);
	mat_mul(Bc,nx,nu,u_k,1,temp2);
	vec_add(temp1,temp2,nx,xm);

	mat_mul(Cc,ny,nx,xm,1,y_k);

	free(temp1);free(temp2);
}

void SP_wright(double *delta_u_ini, double *y_bar, double *lambda_bar, double *G, double *c, double *A,double *A_t, double *b,double *y_ini, double *lambda_ini)
{

	int i;

	double *delta_x_aff = (double*)calloc(nu*M,sizeof(double));
	double *delta_y_aff = (double*)calloc(4*nu*M,sizeof(double));
	double *delta_lambda_aff = (double*)calloc(4*nu*M,sizeof(double));

	mat_vec(G,delta_u_ini,nu*M,nu*M,MPC_temp1);
	mat_vec(A_t,lambda_bar,nu*M,4*nu*M,MPC_temp2);
	vec_sub(MPC_temp1, MPC_temp2,nu*M,MPC_temp3);
	vec_add(MPC_temp3,c,nu*M,rd);
	//printf("rd is:\n");show_matrix(rd,1,n);putchar('\n');

	mat_vec(A,delta_u_ini,4*nu*M,nu*M,MPC_temp1);
	vec_sub(b,MPC_temp1,4*nu*M,rp_y);
	vec_sub(MPC_temp1,y_bar,4*nu*M,MPC_temp2);
	vec_sub(MPC_temp2,b,4*nu*M,rp);
	//printf("rp is:\n");show_matrix(rp,1,m);putchar('\n');

	vec_div(lambda_bar,y_bar,4*nu*M,inv_y_lambda);
	bsxfun(A,inv_y_lambda,4*nu*M,nu*M,y_lambda_A);
	//printf("inv_y_lambda is :\n");show_matrix(inv_y_lambda,1,m);putchar('\n');
	//printf("y_lambda_A is :\n");show_matrix(y_lambda_A,n,m);putchar('\n');

	mat_mul(A_t,nu*M,4*nu*M,y_lambda_A,nu*M,MPC_temp1);
	vec_add(G,MPC_temp1,nu*M*nu*M,equation_F);
	//printf("equation_F is: \n");show_matrix(equation_F,n,n);putchar('\n');

	vec_mul(inv_y_lambda,rp_y,4*nu*M,MPC_temp1);
	mat_vec(A_t,MPC_temp1,nu*M,4*nu*M,MPC_temp2);
	vec_sub(MPC_temp2,rd,nu*M,equation_b);
	//printf("equation_b is: \n");show_matrix(equation_b,1,n);putchar('\n');

	line_solve(equation_F,equation_b,nu*M,delta_x_aff);

	//printf("delta_x is: \n");show_matrix(delta_x,1,n);putchar('\n');

	mat_vec(A,delta_x_aff,4*nu*M,nu*M,MPC_temp1);
	vec_add(MPC_temp1,rp,4*nu*M,delta_y_aff);
	//printf("delta_y is: \n");show_matrix(delta_y,1,m);putchar('\n');

	vec_mul(inv_y_lambda,delta_y_aff,4*nu*M,MPC_temp2);
	vec_rev(MPC_temp2,4*nu*M,MPC_temp3);
	vec_sub(MPC_temp3,lambda_bar,4*nu*M,delta_lambda_aff);
	//printf("delta_lambda is: \n");show_matrix(delta_lambda,1,m);putchar('\n');

	vec_add(y_bar,delta_y_aff,4*nu*M,MPC_temp1);
	vec_abs(MPC_temp1,4*nu*M,y_ini);
	for(i=0;i<4*nu*M;i++)
	{
		if(y_ini[i]<1)
			y_ini[i]=1;
	}

	vec_add(lambda_bar,delta_lambda_aff,4*nu*M,MPC_temp1);
	vec_abs(MPC_temp1,4*nu*M,lambda_ini);
	for(i=0;i<4*nu*M;i++)
	{
		if(lambda_ini[i]<1)
			lambda_ini[i]=1;
	}

	//printf("lambda_ini is: \n");show_matrix(lambda_ini,1,4*nu*M);putchar('\n');
	//printf("y_ini is: \n");show_matrix(y_ini,1,4*nu*M);putchar('\n');
	//printf("delta_u_bar is: \n");show_matrix(delta_u_ini,1,nu*M);putchar('\n');


	free(delta_x_aff);free(delta_y_aff);free(delta_lambda_aff);

}

// KKT residual (Wright) Termination Criteria
double TC_KKT_Wright(double *rd, double *rp, double b_max,double c_max, double mu, double epsilon,int mc, int ndec)
{
	int flag = 0;
	double res_1, res_2, res;
	//printf("rp:\n");show_matrix(rp,1,mc);putchar('\n');
	res_1 = scalar_max(scalar_abs(vec_min(rp,mc)),scalar_abs(vec_max(rp,mc)));
	//printf("res_1: %f\n",res_1);
	res_1 = res_1/b_max;
	//printf("res_1: %f\n\n",res_1);
	//printf("rd:\n");show_matrix(rd,1,ndec);
	res_2 = scalar_max(scalar_abs(vec_min(rd,ndec)),scalar_abs(vec_max(rd,ndec)));
	//printf("res_2: %f\n",res_2);
	res_2 = res_2/c_max;
	//printf("res_2: %f\n\n",res_2);
	res = scalar_max(scalar_max(res_1,res_2),scalar_abs(mu));
	//printf("res: %f\n\n",res);
	if(res<epsilon)
		flag = 1;
	return flag;
	//printf("Res_1: %f\n",res_1);
	//printf("Res_2: %f\n",res_2);
	//printf("mu : %f\n",abs(mu));
	//printf("Residual: %f\n",res);
}

//Sigmoid function. Looking up table form.
// Yi. 2014.8.19
double Sigmoid_tb(double delta, double epsilon)
{
	double S = 250;
	//Table 1
	if(epsilon==0.00001)
	{
		if(delta<0.0001)
			S = 1;
		else if((delta>0.0001)&&(delta<0.01))
			S = 0.8;
		else if((delta>0.01)&&(delta<1))
			S = 0.6;
		else if((delta>1)&&(delta<100))
			S = -0.6;
		else if((delta>100)&&(delta<10000))
			S = -0.8;
		else if((delta>10000)&&(delta<1000000))
			S = -1;
		else
			printf("Error: Delta(%f) range out of bounds!",delta);
	}
	//Table 2
	else if(epsilon==0.0000000001)
	{
		printf("Error: Not done yet!");
	}
	else
		printf("Error: No corresponding table found!");
	return S;
}

//Sigmoid function. Quadratic Polyfit.
// Yi. 2014.10.25
double Sigmoid_pf(double delta, double epsilon)
{
	double S = 250;
	//Table 1
	if(epsilon-0.00000001 < 0.00000000001)		// epsilon = 1e-8
	{
		if(delta <= 0.00000001)
			S = 1.001;
		else if((delta>0.00000001)&&(delta<=0.0000001)) {
			delta = delta * 100000000;
			S = 0.000221 * delta * delta - 0.0070 * delta + 1.0019;
		}
		else if((delta>0.0000001)&&(delta<=0.000001)) {
			delta = delta * 10000000;
			S = 0.000307 * delta * delta - 0.0097 * delta + 0.9581;
		}
		else if((delta>0.000001)&&(delta<=0.00001)) {
			delta = delta * 1000000;
			S = 0.000420 * delta * delta - 0.0131 * delta + 0.8975;
			}
		else if((delta>0.00001)&&(delta<=0.0001)) {
			delta = delta * 100000;
			S = 0.000558 * delta * delta - 0.0173 * delta + 0.8152;
		}
		else if((delta>0.0001)&&(delta<=0.001)) {
			delta = delta * 10000;
			S = 0.000717 * delta * delta - 0.0220 * delta + 0.7067;
		}
		else if((delta>0.001)&&(delta<=0.01)) {
			delta = delta * 1000;
			S = 0.000880 * delta * delta - 0.0267 * delta + 0.5687;
		}
		else if((delta>0.01)&&(delta<=0.1)) {
			delta = delta * 100;
			S = 0.0010 * delta * delta - 0.0306 * delta + 0.4013;
		}
		else if((delta>0.1)&&(delta<=1)) {
			delta = delta * 10;
			S = 0.0011 * delta * delta - 0.0329 * delta + 0.2096;
		}
		else if((delta>1)&&(delta<=10)) {
			delta = delta * 1;
			S = 0.0011 * delta * delta - 0.0330 * delta + 0.0038;
		}
		else if((delta>10)&&(delta<=100)) {
			delta = delta * 0.1;
			S = 0.0011 * delta * delta - 0.0309 * delta - 0.2021;
		}
		else if((delta>100)&&(delta<=1000)) {
			delta = delta * 0.01;
			S = 0.000957 * delta * delta - 0.0272 * delta - 0.3944;
		}
		else if((delta>1000)&&(delta<=10000)) {
			delta = delta * 0.001;
			S = 0.000801 * delta * delta - 0.0225 * delta - 0.5625;
		}
		else if((delta>10000)&&(delta<=100000)) {
			delta = delta * 0.0001;
			S = 0.000637 * delta * delta - 0.0178 * delta - 0.7015;
		}
		else if((delta>100000)&&(delta<=1000000)) {
			delta = delta * 0.00001;
			S = 0.000488 * delta * delta - 0.0135 * delta - 0.8110;
		}
		else if((delta>1000000)&&(delta<=10000000)) {
			delta = delta * 0.000001;
			S = 0.000362 * delta * delta - 0.0100 * delta - 0.8934;
		}
		else if((delta>1000000)&&(delta<=10000000)) {
			delta = delta * 0.000001;
			S = 0.000263 * delta * delta - 0.0073 * delta - 0.9557;
		}
		else if(delta>100000000)
			S = -1.01;
		else
			printf("Error: Delta(%f) range out of bounds!",delta);
	}
	//Table 2
	else if(epsilon==0.000001)
	{
		printf("Error: Not done yet!");
	}
	else
		printf("Error: No corresponding table found!");

	return S;
}


//Convergence Depth Control
// Yi. 2014.8.19
int CDC(double cons_obey, double cons_obey_old, double Fx, double Fx_old, double *x, double *delta_x, double *G, double *c,double mu,double theta_0,double theta_1)
{
	double epsilon_0  = 0.00000001;
	int eta = 3;
	double delta_objErr,delta_feasChg,delta_objChg;
	double theta_convg, theta_prog, theta;


	int counter = 0;
	int flag = 0;


	mat_vec(G,x,ndec,ndec,CDC_temp1);
	vec_add(CDC_temp1,c,ndec,CDC_temp2);
	delta_objErr = dot_product(CDC_temp2,delta_x,ndec);

	delta_feasChg = scalar_abs(cons_obey-cons_obey_old);
	delta_objChg = scalar_abs(Fx-Fx_old);

	//printf("Sigmoid 1 parameter: %f\n",scalar_max(scalar_max(cons_obey,delta_objErr),mu));
	//printf("Sigmoid 2 parameter: %f\n\n",scalar_max(scalar_max(delta_feasChg,delta_objChg),mu));
	theta_convg = Sigmoid_pf(scalar_max(scalar_max(cons_obey,delta_objErr),mu),epsilon_0);
	theta_prog = Sigmoid_pf(scalar_max(scalar_max(delta_feasChg,delta_objChg),mu),epsilon_0);

	if(theta_convg>theta_0)
		flag = 1;
	else
	{
		if(theta_prog>theta_1)
		{
			counter++;
			if(counter>eta)
				flag = 2;
		}
	}

	return flag;



}

void get_state(double *x_k,double *y_k)
{
	int i;
	for(i=0;i<nx;i++)
		x_k[i] = xm[i] - xm_old[i];
	for(i=nx;i<nx+ny;i++)
		x_k[i] = y_k[i-nx];
}
