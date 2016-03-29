
 //*
 //* mpc_lib.h
 //*
 //*  Created on: 2014-3-29
 //*      Author: Yi
 //*
 //*  Notation:
 //*		sca: scalar;
 //*		vec: vector;
 //*		mat: matrix;
 //*		mul: multipication;
 //*
 //*



#ifndef MPC_LIB_H_
#define MPC_LIB_H_



//Find the absloute value of a double parameter.
//2014.4.10 Yi Ding.
double scalar_abs(double x);

//Find the minimun one of two scalars.
//2014.4.9 Yi Ding.
double scalar_min(double x,double y);

//Find the maximun one of two scalars.
//2014.4.11 Yi Ding.
double scalar_max(double x,double y);

//Filling a vector with a specific scalar.
//double *sca2vec(double sca, int n);
void sca2vec(double sca, int n,double *vec);

//Find the minimum one of a vector.
////2014.4.9 Yi Ding.
double vec_min(double *x, int n);

//Find the maximum one of a vector.
////2014.4.10 Yi Ding.
double vec_max(double *x, int n);

//Vector augment!
void vec_aug(double *vec,int m, int n,double *aug);

//Vector addation
//double *vec_add(double *a, double *b, int n);
void vec_add(double *a, double *b, int n,double *c);

//Vector subtraction
//double *vec_sub(double *a, double *b, int n);
void vec_sub(double *a, double *b, int n,double *c);

//Calculate |a| of a
//2014.5.14 Yi Ding.
//double *vec_abs(double *a, int n)
void vec_abs(double *a, int n, double *b);

// Return the product of a scalar and vector;
//double *sca_vec_mut(double sca, double *vec, double n);
void sca_vec_mut(double sca, double *vec, double n, double *pro);

//This is a function to make a scalar[a] into a vecotr and substract another vector vec_b to get vec_x
//Yi Ding. 2014.5.21
void sca_vec_sub(double a, double *vec_b, int length, double *vec_x);

//This is a function to do .* operation to 2 vectors.
//Not dot product, but the corresponding element to do multiplication and make a new vector.
//Yi Ding. 2014.5.21
void vec_mul(double *a, double *b, int m, double *x);

//This is a function to do ./ operation to 2 vectors.
//Let the corresponding element to do division and make a new vector.
//Yi Ding. 2014.5.21
void vec_div(double *a, double *b, int m, double *x);

//This function implement the m-fun bsxfun() is Matlab, but only implement the mat-vec times function.
//Not the same usage with that in Matlab, be careful to use.
//Yi Ding. 2014.5.21
void bsxfun(double *A, double *x, int row, int col, double *B);

//Calculate -a of a
//double *vec_rev(double *a, int n);
void vec_rev(double *a, int n, double *b);

//Make the vector into a diagonal form.
//double *vec_diag(double *a, int n);
void vec_diag(double *a, int n, double *diag);

// Extract the diagonal elements of a diagonal matrix and form a vector.
//double *diag2vec(double *diag, int n);
void diag2vec(double *diag, int n , double *vec);

//Calculate the inverse of a diagonal matrix.
//double *inv_diag(double *diag, int n);
void inv_diag(double *diag, int n ,double *inv_diag);

//Add two matrices
//double *mat_add(double *a, double *b, int row, int column);
void mat_add(double *a, double *b, int row, int column, double *sum);

//Do the multiplication A*x = b;
//Yi Ding. 2014.4.25
void mat_vec(double *A, double *x, int row, int col, double *b);

//Do the multiplication x'*A = b;
//Yi Ding. 2014.7.23
void vec_mat(double *A, double *x, int row, int col, double *b);


//Do the multiplication A*diag(x) = b;
//Yi Ding. 2014.4.26
void mat_diag(double *A, double *diag, int row, int col, double *b);

//Used for LU Evaluate!
//2014.5.15
void Lower_diag_vec(double *A, double *diag_vec, int m, double *b);

//Do the A*X = B, here A, X, B are all diagnal matrix
//Yi Ding 2014.4.26
void diag_diag(double *A, double *X, int n, double *B);

//Do the A*x = b, here A is diag matrix, x is a vector
//Yi Ding 2014.4.26
void diag_vec(double *A, double *x, int n, double *b);

//Matrix multiplication. Just like DSP_sp_mat_mul.c in ti DSPlib.
//double *mat_mul(double *x1, int r1, int c1, double *x2, int c2);
void mat_mul(double *x1, int r1, int c1, double *x2, int c2, double *y);

//Matrix multiplication. For off line use. For the accurate running time of on line mat_mul.
//2014.7.18 Yi Ding.
void mat_mul_offline(double *x1, int r1, int c1, double *x2, int c2, double *y);

//Compute the power of a matrix
void mat_pow(double *mat, int row, int pow, double *MatPow);

//Cholesky decomposition.
//double *cholesky(double *A, int n);
void cholesky(double *A, int n, double *L);

//A Cholesky Factorization routine from the book Numerical Recipes in C
// 2014.7.15 Yi Ding
int chol_NRC(double *A, int n,double *p);

//Here, the diagnal D is stored in a vector.
//Modified Cholesky Factorization 2014.5.15
void mcf(double *A, int m, double *L, double *D);

//Dot product
double dot_product(double *a, double *b, int n);

//Call two functions to solve linear function.
//double *line_solve(double *A, double *b, int n);
void line_solve(double *A, double *b, int n, double *x);

//This function implement a special searching method for alpha.
//The idea comes from the quad_wright() and symbol meanings see algorithms.
// Yi Ding. 2014.5.21
double alpha_decreas(double *vec, double *delta_vec, int length,  double tau);

// Alpha calculation method (while loop)
double alpha_decreas_while(double *y, double *lambda, double *delta_y, double *delta_lambda, int m);

//Do back substatution to solve linear equation using LU.
//double *luEvaluate(double *L,double *U,double *b,int n);
void luEvaluate(double *L,double *U, double*b,int n,double *x);

// Do back substatution to solve linear equation using LU.
// A routine from the book Numerical Recipes in C
// Yi Ding 2014.7.15
void luEvaluate_NRC(double *L, double*b,double *p,int n,double *x);

//Print the matrix to the consolo. Matrix stored in a array. w column, h rows.
void show_matrix(double *m, int w, int h);

//Transpose the matrix m and store the matrix in w.
void transpose(double *m, int w, int h);

//Primal Dual Interior Point Method to solve QP.
void priduip(double *G, double *c, double *A, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u);

//The feedback step in MPC, a transformation of funciton fankui.m
void feedback_v3(double *x_k_1, double *y_k, double *u_k_1, double *u_k_2, double *x_k,double *L);

//Kalman filter
void feedback_kalman(double *x_k_1, double *y_k, double *u_k_1, double *u_k_2, double *x_k,int iter);

//The formation of omega_r.
//This version has y contraints.
void omega_r_form(double *aug_u_k_1, double *x_k, double *omega_r, double*F);

//The formation of omega_r.
//This version does not have y contraints.
void omega_r_form2(double *aug_u_k_1, double *x_k, double *omega_r, double *F);

//The formation of omega_r.
//This version has y contraints but do not have delta_u constraints
void omega_r_form3(double *aug_u_k_1, double *x_k, double *omega_r, double *F,double Y_p1, double Y_p2, double Y_n1, double Y_n2);

//Find the parameter F and Phi, not in the MPC loop. Relaxation on performance.
//Yi Ding. 2014.4.16
void fphi(double *F, double *Phi);

//Form the B used in the OMEGA_L
void form_B(double *B);

//Form the OMEGA_L for the constrains in the optimization problems.
//This version does not have y constraints.
void OMEGA_L_form(double *B, double *Phi, double *OMEGA_L);

//Form the OMEGA_L for the constrains in the optimization problems.
//This version does not have y constraints.
void OMEGA_L_form2(double *B, double *Phi, double *OMEGA_L);

//Form the OMEGA_L for the constrains in the optimization problems.
//This version has y constraints. But DO not have delta_u constraints.
void OMEGA_L_form3(double *B, double *Phi, double *OMEGA_L);

//Form the G for the opimization problem.
void form_G(double *Phi,double *G);

void form_G2(double *Phi,double *G);

//This is a routine to mimic the model using the disturbed state-space model.
void process_model(double* xm, double *u_k, double *y_k);

void process_model2(double* xm, double *u_k, double *y_k);

void process_model3(double* xm, double *u_k, double *y_k);

//MPC online computing!
//Yi Ding @2014.4.18
void mpc(double *delta_u, double *u_k, double *y_k,double *x_k,double *r_k, double *u_k_1, double *u_k_2,double *F,double *Phi,double *OMEGA_L,double *N_OMEGA_L,double *L,double *G);

void SP_wright(double *delta_u_ini, double *y_bar, double *lambda_bar, double *G, double *c, double *A,double *A_t, double *b,double *y_ini, double *lambda_ini);

double TC_KKT_Wright(double *rd, double *rp, double b_max,double c_max, double mu, double epsilon,int mc, int ndec);

//Sigmoid function. Looking up table form.
// Yi. 2014.8.19
double Sigmoid_tb(double delta, double epsilon);

//Sigmoid function. Quadratic Polyfit.
// Yi. 2014.10.25
double Sigmoid_pf(double delta, double epsilon);

//Convergence Depth Control
// Yi. 2014.8.19
int CDC(double cons_obey, double cons_obey_old, double Fx, double Fx_old, double *x, double *delta_x, double *G, double *c,double mu,double theta_0,double theta_1);

void get_state(double *x_k,double *y_k);

#endif /* MPC_LIB_H_ */
