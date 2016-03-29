//
 #include <stdio.h>
 #include <string.h>
 #include <stdlib.h>
 #include <math.h>
 #include <limits.h>
 #include <stdint.h>
 //#include <time.h>

 #include <xdc/std.h>
 #include <xdc/runtime/Error.h>
 #include <xdc/runtime/System.h>
 #include <xdc/runtime/Types.h>
 #include <xdc/runtime/Memory.h>
 #include <xdc/runtime/Timestamp.h>
 #include <ti/sysbios/BIOS.h>
 #include <ti/sysbios/knl/Task.h>

 #include "ti\platform\platform.h"

 #include "mpc_lib_DP.h"
 #include "PDIPM_V3_DP_Opt.h"
 #include "IPM_V3_DP.h"


 //Here are some global variables.
 int P = 25;		//Prediction horizon
 int M = 5;		//Control horizon
 double Q[] = {1.5,0};	//
 double R = 2;

 int nu = 1;		//The number of inputs
 int ny = 2;		//The number of outputs
 int nx = 4;		//The number of original states.(Augment model has one more state)
 int nyc = 2;	//The number of constrained outputs

 int ndec;		//Number of decision variables in QP solving.
 int mc;			//?Number of inequality constraints in QP solving.

 // Model parameters
 double A_e[] = {1,0.05,-0.0057,-0.000094883,0,0,0,1,-0.2308,-0.0057,0,0,0,0,1.0593,0.0510,0,0,0,0,2.3968,1.0593,0,0,1,0.05,-0.0057,-0.000094883,1,0,0,0,1.0593,0.0510,0,1};
 double B_e[] = {0.0028,0.1106,-0.0091,-0.3674,0.0028,-0.0091};
 double C_e[] = {0,0,0,0,1,0,0,0,0,0,0,1};
 double C_e1[] = {0,0,0,0,1,0,0};

 //Constraint initialization
 double delta_U_p = 5;
 double delta_U_n = -5;
 double U_p = 5;
 double U_n = -5;
 double Y_p;
 double Y_n;
 double Y_p1 = 2, Y_p2 = 0.79, Y_n1 = -2, Y_n2 = -0.79;

 //VARIABLES DECLARARION
 //Variable declaration for main
 //Constant in online computing.
 double *L,*B,*G,*GL,*GL_T,*F,*Phi,*OMEGA_L,*N_OMEGA_L,*N_OMEGA_L_T,*F1,*Phi1,*QQ,*Q_temp;
 //Constant in online computing.
 double *delta_u_ini,*y_ini,*lambda_ini;

 //Variable declaration for mpc_DP
 double *aug_u_k_1;
 double *omega_r,*n_omega_r;
 double *c,*c_rev;
 double *h_x;
 double *delta_u_M; //delta_u in the M control horizon
 double *MPC_temp1,*MPC_temp2,*MPC_temp3,*MPC_temp4;
 double b_max,c_max,G_max,A_max,p_max;

 //Make xm a global variable
 double xm[] = {0,0,0,0};
 double xm_old[] = {0,0,0,0};

 //Variable declaration for IPM_v3_DP
 int MaxIteration = 30;
 double *x_old;
 double *rd, *rp, *rp_y;
 double *inv_y_lambda, *y_lambda_A;
 double *equation_F, *equation_b;
 double *delta_x, *delta_y, *delta_lambda;
 double *temp_cond;	//To check condition in the choice of alpha_tau_pri & alpha_tau_dual.
 double *x_try;
 double *sig_mu, *center_part;
 double *QP_temp1, *QP_temp2, *QP_temp3;

 //Additional variable delcartion for priduip_v4
 //double *delta_x_aff, *delta_y_aff, *delta_lambda_aff;
 //double *sig_mu_lam_y, *aff_part;

 //Variable declaration for feedback_v3
 double *fb_temp1, *fb_temp2;

 //Variable declaration for mcf
 double *mcf_C;
 double *mcf_temp_vec;

 //Variable declaration for alpha_decreas
 double *alpha_decreas_neg_delta_value, *alpha_decreas_neg_vec_value;
 double *alpha_decreas_temp1,*alpha_decreas_temp2;

 //Variable declaration for line_solve
 double *line_solve_L,*line_solve_p;
 double *equation_F_scaling,*equation_b_scaling;
 double *line_solve_U,*line_solve_D,*line_solve_LD;

 //Variable declaration for CDC
 double *CDC_temp1,*CDC_temp2;

 //Time consuming record
 int QP_Iteration = 0;
 double Time_MPC = 0,Time_QP = 0,Time_QP_non_Iter = 0,Time_QP_Iter = 0;
 double Time_Section1 = 0,Time_Section2 = 0,Time_Section3 = 0,Time_Section4 = 0,Time_Section5 = 0;
 //For priduip.
 double Time_Section6 = 0,Time_Section7= 0,Time_Section8 = 0;
 double Time_LineSolve1 = 0,Time_LineSolve2 = 0;

 //Time consuming for subfunctions
 int mat_mul_times = 0;		double Time_Mat_mul = 0;
 int Mat_vec_times = 0;
 int vec_add_times = 0;		double Time_vec_add = 0;
 int vec_sub_times = 0;		double Time_vec_sub = 0;
 int vec_div_times = 0;		double Time_vec_div = 0;
 int vec_mul_times = 0;		double Time_vec_mul = 0;
 int bsxfun_times = 0; 		double Time_bsxfun = 0;
 int sca2vec_times = 0;		double Time_sca2vec = 0;
 int vec_rev_times = 0;		double Time_vec_rev = 0;
 int sca_vec_mut_times = 0;	double Time_sca_vec_mut = 0;
 int dot_product_times = 0;	double Time_dot_product = 0;
 int alpha_compute_times = 0;

 int total_QP_iterations = 0;		double Time_MPC_Total = 0;
 int worst_QP_iterations = 0;		double Time_MPC_Worst = 0;

 //// Control performance record
 double IAE = 0;
 double ITAE = 0;
 double ISE = 0;

 FILE *fp_Fx, *fp_cons_obey;


 //============================mpc_lib_DP.c==START======================================

 //============================mpc_lib_DP.c==END========================================


 //
 //Read seria port：
 int post_read_uart(uint8_t* msg, uint32_t msg_len,uint32_t delay)
 {
     uint32_t i;

	 for(i=0;i<msg_len;i++)
	 {
		 if(platform_uart_read(&(msg[i]),delay) != Platform_EOK)
		 {
			 return -1;
		 }
	 }
	 return 0;
 }

 //Write seria port：
 int post_write_uart( char* msg,uint32_t msg_len)
 {
     uint32_t i;

     for (i = 0; i < msg_len; i++)
     {
         if (platform_uart_write(msg[i]) != Platform_EOK)
         {
             return -1;
         }
     }
     return 0;
 }


 //MPC online computing!
 void mpc_V3_DP(double *delta_u, double *delta_u_ini,double *y_ini,double *lambda_ini, double *u_k, double *y_k,double *x_k,double *r_k, double *u_k_1, double *u_k_2,int iter)
 {
	 int i;
	 double WarmVal = 0.5;

	 UInt32 start_time,end_time,dead_time;
	 UInt32 QP_start_time,QP_end_time,QP_total_time;

     start_time = Timestamp_get32();
     end_time = Timestamp_get32();
     dead_time = end_time - start_time;

	 // State Updating
	 //feedback_v3(x_k,y_k,u_k_1,u_k_2,x_k,L);
	 //show_matrix(x_k,1,nx+ny);

	 //Use kalman filter to update states
	 feedback_kalman(x_k,y_k,u_k_1,u_k_2,x_k,iter);
	 //printf("x_k is:\n");show_matrix(x_k,1,nx+ny);putchar('\n');

	 ////Or get real state
	 //get_state(x_k,y_k);
	 //printf("x_k is:\n");show_matrix(x_k,1,nx+ny);putchar('\n');

	 vec_aug(u_k_1,nu,M,aug_u_k_1);

	 //Below form the omega_r with y constraints.
	 omega_r_form3(aug_u_k_1, x_k, n_omega_r,F, Y_p1, Y_p2, Y_n1, Y_n2);
	 vec_rev(n_omega_r,mc,n_omega_r);
	 //printf("omega_r is: \n");show_matrix(omega_r,1,mc);putchar('\n');

	 ////Below form the omega_r without y constraints.
	 //omega_r_form2(aug_u_k_1, x_k, omega_r,F);
	 //vec_rev(omega_r,4*nu*M,n_omega_r);
	 ////show_matrix(omega_r,1,4*nu*M+2*ny*P);

	 //Below form the c in QP.
	 mat_vec(F,x_k,ny*P,nx+ny,MPC_temp1);
	 //printf("MPC_temp1 is:\n");show_matrix(MPC_temp1,1,ny*P);putchar('\n');
	 vec_sub(MPC_temp1,r_k,ny*P,MPC_temp2);
	 //printf("MPC_temp2 is:\n");show_matrix(MPC_temp2,1,ny*P);putchar('\n');
	 //vec_mat(QQ,MPC_temp2,ny*P,ny*P,MPC_temp3);
	 //vec_mat(Phi,MPC_temp3,ny*P,nu*M,c);
	 bsxfun(Phi,Q_temp,P*ny,M*nu,MPC_temp3);
	 //printf("MPC_temp3 is:\n");show_matrix(MPC_temp3,M*nu,ny*P);putchar('\n');
	 vec_mat(MPC_temp3,MPC_temp2,P*ny,M*nu,c);
	 ////c scaling
	 //for(i=0;i<ndec;i++)
	 //	c[i] = c[i]*0.00001;
	 vec_rev(c,nu*M,c_rev);
	 //printf("c is:\n");show_matrix(c,1,nu*M);putchar('\n');
	 //printf("c_rev is:\n");show_matrix(c_rev,1,nu*M);putchar('\n');

	 //QP Starting Point initializing (delta_u)
	 for(i=0;i<ndec;i++)
		 delta_u_ini[i] = delta_u_ini[i+1];
	 delta_u_ini[ndec-1] = 0;

	 // Below is a parameter used in an initializing method from LOQO
	 //mat_vec(OMEGA_L,delta_u_ini,4*nu*M+2*ny*P,nu*M,MPC_temp4);
	 //vec_sub(omega_r,temp4,4*nu*M+2*ny*P,h_x);

	 //QP Starting Point initializing (A heuristic method for y and lambda)
	 for(i=0;i<mc;i++)
	 {
		 y_ini[i] = scalar_max(y_ini[i],0.5);
		 lambda_ini[i]  = scalar_max(lambda_ini[i],0.5);
	 }

	 //QP Starting Point initializing (A simple method for y and lambda)
	 //for(i=0;i<4*nu*M+2*ny*P;i++)
	 //{
	 //	y_ini[i] = 0.5;
	 //	lambda_ini[i]  = 0.5;
	 //}


	 ////QP Starting Point initializing (AUT's method)
	 // Scale preparation
	 vec_abs(n_omega_r,mc,MPC_temp1);
	 b_max = vec_max(MPC_temp1,mc);
	 vec_abs(c,ndec,MPC_temp1);
	 c_max = vec_max(MPC_temp1,ndec);
	 vec_abs(G,ndec*ndec,MPC_temp1);
	 G_max = vec_max(MPC_temp1,ndec*ndec);
	 vec_abs(N_OMEGA_L,mc*ndec,MPC_temp1);
	 A_max = vec_max(MPC_temp1,mc*ndec);
	 p_max = scalar_max(scalar_max(b_max,c_max),scalar_max(G_max,A_max));
	 //printf("b_max, c_max is :%f,%f\n",b_max,c_max);
	 //printf("p_max: %f\n",p_max);
	 //if(p_max>1)
	 //	WarmVal = sqrt(p_max);
	 //for(i=0;i<4*nu*M+2*ny*P;i++)
	 //{
	 //	y_ini[i] = WarmVal;
	 //	lambda_ini[i]  = WarmVal;
	 //}

	 ////QP Starting Point initializing (Wright's method)
	 //SP_wright(delta_u_ini, y_ini, lambda_ini, G, c, N_OMEGA_L,N_OMEGA_L_T, n_omega_r,y_ini, lambda_ini);
	 //printf("delta_u_ini is: \n");show_matrix(delta_u_ini,1,nu*M);putchar('\n');
	 //printf("y_ini is: \n");show_matrix(y_ini,1,4*nu*M);putchar('\n');
	 //printf("lambda_ini is: \n");show_matrix(lambda_ini,1,4*nu*M);putchar('\n');

	 ////The below is a way to initialize lambda with a random value between 0 and 1.
	 ////In the process of debug with Matlab, the lambda_ini is define just like which in Matlab
	 //lambda_ini = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
	 //for(i=0;i<4*nu*M+2*ny*P;i++)
	 //	lambda_ini[i] = ((double)rand()/RAND_MAX);
	 //Show initial values.
	 //putchar('\n');show_matrix(delta_u_ini,1,nu*M);
	 //putchar('\n');show_matrix(y_ini,1,4*nu*M);
	 //putchar('\n');show_matrix(lambda_ini,1,4*nu*M);
	 // This random method can be used to compare with other method.


	 ////QP solve
	 QP_start_time = Timestamp_get32();
	 IPM_V3_DP(G,GL,GL_T,c, N_OMEGA_L,N_OMEGA_L_T, n_omega_r, delta_u_ini, y_ini, lambda_ini, mc, ndec, delta_u_M);
	 //priduip_v4_DP(G, GL,GL_T, c, N_OMEGA_L, N_OMEGA_L_T, n_omega_r, delta_u_ini, y_ini, lambda_ini, mc, ndec, delta_u_M);
	 //printf("delta_u is: \n");show_matrix(delta_u_M,1,nu*M);
	 memcpy(delta_u,delta_u_M,nu*sizeof(double));
	 QP_end_time = Timestamp_get32();
	 QP_total_time = QP_end_time - QP_start_time - dead_time;
	 Time_QP = (double)(QP_total_time);

	 //vec_add(u_k,delta_u,nu,MPC_temp1);
	 for(i=0;i<nu;i++)
		 MPC_temp1[i] = u_k[i]+delta_u[i];
	 memcpy(u_k,MPC_temp1,nu*sizeof(double));
 }

 Void taskFxn(UArg a0, UArg a1)
 {
     printf("enter taskFxn()\n");
     int i,j;
     int maxmem;

	 //double *u_k,*delta_u,*x_k,*y_k,*r_k;
	 //double *u_k,*u_k_1,*u_k_2,*delta_u,*x_k,*y_k,*xm;

	 //The variables below are used to test time consuming.
	 double u_k[] = {0};
	 double u_k_1[] = {0};
	 double u_k_2[] = {0};
	 double delta_u[] = {0.1};
	 double x_k[] = {0,0,0,0,0,0,0};
	 double y_k[] = {0,0,0};

	 //Reference curve
	 double r_k1[] = {1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0,1,0};
	 double r_k2[] = {-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0,-1,0};
	 //double r_k2[] = {1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1};

	 UInt32 start_time,end_time,dead_time;
	 UInt32 MPC_start_time,MPC_end_time,MPC_total_time;
	 Types_FreqHz freq;

	 FILE *fp_u, *fp_y1, *fp_y2, *fp_y3, *fp_t, *fp_i,*fp_mpct;
	 FILE *fp_qpt_nonIter,*fp_qpt_Iter;
	 FILE *fp_s1t, *fp_s2t, *fp_s3t, *fp_s4t, *fp_s5t;

	 //Communication Settings
	 char token_read = 255;
	 char token_write = 240;
	 float u_float = 0;
	 float y1_float = 0;
	 float y2_float = 0;
	 double y1_double = 0;
	 double y2_double = 0;
	 float count_write = 1;
	 float count_read = 1;
	 float count_r_verify;
	 int t_f1, t_f2;
	 //uint32_t delay = 1700000;	//摸索出来的，不要随便乱改
	 uint32_t delay = 2000000;	//摸索出来的，不要随便乱改

	 uint8_t Buffer_read[12];
	 char Buffer_write[8];


	 Timestamp_getFreq(&freq);
	 //printf("freq = {0x%x : 0x%x}\n", freq.hi, freq.lo);
	 printf("Frequency is: %u\n",freq.lo);

	 //Determine the number of inequality constraints and number of decision values
	 ndec = nu*M;
	 //mc = 4*nu*M;	//Outputs constraints not included.
	 mc = 2*nu*M+2*nyc*P;	//Outputs constraints included.

	 //Memory allocation!============================================================
	 // Memory allocation for main function
	 L = (double*)calloc((nx+ny)*ny,sizeof(double));
	 F = (double*)calloc(P*ny*(nx+ny),sizeof(double));
	 Phi = (double*)calloc(M*nu*P*ny,sizeof(double));
	 B = (double*)calloc(P*nu*M*nu*M,sizeof(double));
	 //OMEGA_L = (double*)calloc(mc*ndec,sizeof(double));
	 N_OMEGA_L = (double*)calloc(mc*ndec,sizeof(double));
	 N_OMEGA_L_T = (double*)calloc(mc*ndec,sizeof(double));
	 G = (double*)calloc(ndec*ndec,sizeof(double));
	 GL = (double*)calloc(ndec*ndec,sizeof(double));
	 GL_T = (double*)calloc(ndec*ndec,sizeof(double));
	 delta_u_ini = (double*)calloc(ndec,sizeof(double));
	 y_ini = (double*)calloc(mc,sizeof(double));
	 lambda_ini = (double*)calloc(mc,sizeof(double));
	 //F1 = (double*)calloc(P*(nx+ny),sizeof(double));
	 //Phi1 = (double*)calloc(M*nu*P,sizeof(double));
	 //QQ = (double*)calloc(P*ny*P*ny,sizeof(double));
	 Q_temp = (double*)calloc(P*ny,sizeof(double));

	 // Memory allocation for mpc_DP
	 aug_u_k_1 = (double*)calloc(ndec,sizeof(double));
	 //omega_r = (double*)calloc(mc,sizeof(double));
	 n_omega_r = (double*)calloc(mc,sizeof(double));
	 c = (double*)calloc(ndec,sizeof(double));
	 c_rev = (double*)calloc(ndec,sizeof(double));
	 //h_x = (double*)calloc(4*nu*M+2*ny*P,sizeof(double));
	 delta_u_M = (double*)calloc(ndec,sizeof(double));
	 maxmem = scalar_max(scalar_max(ny*P*(nx+ny),mc),scalar_max(ndec*ndec,mc*ndec));
	 MPC_temp1 = (double*)calloc(maxmem,sizeof(double));
	 MPC_temp2 = (double*)calloc(ny*P,sizeof(double));
	 MPC_temp3 = (double*)calloc(ny*P*M*nu,sizeof(double));
	 //MPC_temp4 = (double*)calloc(mc*mc,sizeof(double));

	 //Memory allocation for IPM_V3_DP
	 //Error may emerge if the temp variable exceed the allocated area!
	 x_old = (double*)calloc(ndec,sizeof(double));
	 rd = (double*)calloc(ndec,sizeof(double));
	 rp = (double*)calloc(mc,sizeof(double));
	 rp_y = (double*)calloc(mc,sizeof(double));
	 inv_y_lambda = (double*)calloc(mc,sizeof(double));
	 y_lambda_A = (double*)calloc(mc*ndec,sizeof(double));
	 equation_F = (double*)calloc(ndec*ndec,sizeof(double));
	 equation_b = (double*)calloc(ndec,sizeof(double));
	 delta_x = (double*)calloc(ndec,sizeof(double));
	 delta_y = (double*)calloc(mc,sizeof(double));
	 delta_lambda = (double*)calloc(mc,sizeof(double));
	 //temp_cond = (double*)calloc(mc,sizeof(double));
	 x_try = (double*)calloc(ndec,sizeof(double));
	 sig_mu = (double*)calloc(mc,sizeof(double));
	 center_part = (double*)calloc(mc,sizeof(double));
	 maxmem = scalar_max(mc,ndec*ndec);
	 QP_temp1 = (double*)calloc(maxmem,sizeof(double));
	 QP_temp2 = (double*)calloc(mc,sizeof(double));
	 QP_temp3 = (double*)calloc(mc,sizeof(double));

	 //Additional Memory allocation for priduip_v4
	 //Error may emerge if the temp variable exceed the allocated area!
	 //delta_x_aff = (double*)calloc(ndec,sizeof(double));
	 //delta_y_aff = (double*)calloc(mc,sizeof(double));
	 //delta_lambda_aff = (double*)calloc(mc,sizeof(double));
	 //diag = (double*)calloc(m*m,sizeof(double));
	 //sig_mu_lam_y = (double*)calloc(mc,sizeof(double));
	 //aff_part = (double*)calloc(mc,sizeof(double));

	 //Memory allocation for feedback_v3
	 fb_temp1 = (double*)calloc(nx+ny,sizeof(double));
	 fb_temp2 = (double*)calloc(nx+ny,sizeof(double));

	 //Memory allocation for mcf
	 mcf_C = (double*)calloc(ndec*ndec+1,sizeof(double));
	 mcf_temp_vec = (double*)calloc(ndec+1,sizeof(double));

	 //Memory allocation for alpha_decreas
	 alpha_decreas_neg_delta_value = (double*)calloc(mc,sizeof(double));
	 alpha_decreas_neg_vec_value = (double*)calloc(mc,sizeof(double));
	 alpha_decreas_temp1 = (double*)calloc(mc,sizeof(double));
	 alpha_decreas_temp2 = (double*)calloc(mc,sizeof(double));

	 //Memory allocation for line_solve
	 line_solve_L = (double*)calloc(ndec*ndec, sizeof(double));
	 line_solve_p = (double*)calloc(ndec, sizeof(double));
	 //equation_F_scaling = (double*)calloc(ndec*ndec, sizeof(double));
	 //equation_b_scaling = (double*)calloc(ndec, sizeof(double));
     line_solve_U = (double*)calloc(ndec*ndec, sizeof(double));
	 line_solve_D = (double*)calloc(ndec, sizeof(double));
	 line_solve_LD = (double*)calloc(ndec*ndec, sizeof(double));

	 //Memory allocation for line_solve
	 CDC_temp1 = (double*)calloc(ndec,sizeof(double));
	 CDC_temp2 = (double*)calloc(ndec,sizeof(double));

	 //Memory allocation Done!============================================================

	 start_time = Timestamp_get32();
	 end_time = Timestamp_get32();
	 dead_time = end_time - start_time;

	 //Initialization of Feedback gain coefficient L.
	 j = ny+1;
	 for(i=nx*ny;i<(nx+ny)*ny;i+=j){
		 L[i] = 1;
	 }
	 printf("L is: \n");show_matrix(L,ny,(nx+ny));putchar('\n');

	 //Below for the QQ
	 vec_aug(Q,ny,P,Q_temp);
	 //vec_diag(Q_temp,ny*P,QQ);
	 //printf("QQ is: \n");show_matrix(QQ,ny*P,ny*P);putchar('\n');

	 //Below form the MPC parameters F, Phi and B.
	 fphi(F,Phi);
	 //fphi1(F1,Phi1);
	 form_B(B);
	 printf("F is: \n");show_matrix(F,nx+ny,ny*P);putchar('\n');
	 printf("Phi is: \n");show_matrix(Phi,nu*M,ny*P);putchar('\n');
	 //printf("F1 is: \n");show_matrix(F1,nx+ny,P);putchar('\n');
	 //printf("Phi1 is: \n");show_matrix(Phi1,nu*M,P);putchar('\n');
	 //show_matrix(B,nu*M,nu*M);

	 ////Below form the OMEGA_L with y constraints
	 OMEGA_L_form3(B,Phi,N_OMEGA_L);
	 vec_rev(N_OMEGA_L,mc*ndec,N_OMEGA_L);
	 memcpy(N_OMEGA_L_T,N_OMEGA_L,mc*ndec*sizeof(double));
	 transpose(N_OMEGA_L_T, ndec, mc);
	 //printf("OMEGA_L is: \n");show_matrix(OMEGA_L,ndec,mc);putchar('\n');
	 printf("N_OMEGA_L is: \n");show_matrix(N_OMEGA_L,ndec,mc);putchar('\n');
	 //printf("N_OMEGA_L_T is: \n");show_matrix(N_OMEGA_L_T,mc,ndec);putchar('\n');

	 ////Below form the OMEGA_L without y constraints
	 //OMEGA_L_form2(B,Phi,OMEGA_L);
	 //vec_rev(OMEGA_L,(4*nu*M)*nu*M,N_OMEGA_L);
	 //memcpy(N_OMEGA_L_T,N_OMEGA_L,(4*nu*M)*nu*M*sizeof(double));
	 //transpose(N_OMEGA_L_T, nu*M, (4*nu*M));
	 ////show_matrix(OMEGA_L,nu*M,(4*nu*M+2*ny*P));

	 //Below form the QP parameters G, GL and GL_T.
	 //form_G(Phi,G);
	 form_G2(Phi,G);
	 //printf("G is: \n");show_matrix(G,ndec,ndec);putchar('\n');
	 //G scaling
	 //for(i=0;i<ndec*ndec;i++)
	 //	G[i] = G[i]*0.00001;
	 cholesky(G,ndec,GL);
	 memcpy(GL_T,GL,ndec*ndec*sizeof(double));
	 transpose(GL_T,ndec,ndec);
	 printf("G is: \n");show_matrix(G,ndec,ndec);putchar('\n');
	 //printf("GL_T is: \n");show_matrix(GL_T,ndec,ndec);putchar('\n');

	 //Initialization of QP parameters in main function.
	 for(i=0;i<ndec;i++)
		 delta_u_ini[i] = 0.1;
	 for(i=0;i<mc;i++)
	 {
		 y_ini[i] = 0.5;
		 lambda_ini[i] = 0.5;
	 }

	 // Runtime data record preparation.
	 fp_u=fopen("DSP_u_k.txt","w");
	 fp_y1=fopen("DSP_y_k1.txt","w");fp_y2=fopen("DSP_y_k2.txt","w");fp_y3=fopen("DSP_y_k3.txt","w");
	 fp_t=fopen("DSP_QPtime.txt","w");fp_i=fopen("DSP_QPIter.txt","w");
	 fp_mpct=fopen("DSP_MPCtime.txt","w");
	 fp_qpt_nonIter = fopen("Time_QP_non_Iter.txt","w");fp_qpt_Iter = fopen("Time_QP_Iter.txt","w");
	 fp_s1t = fopen("Time_QP_Section1.txt","w");
	 fp_s2t = fopen("Time_QP_Section2.txt","w");fp_s3t = fopen("Time_QP_Section3.txt","w");
	 fp_s4t = fopen("Time_QP_Section4.txt","w");fp_s5t = fopen("Time_QP_Section5.txt","w");
	 fp_Fx = fopen("Fx.txt","w");fp_cons_obey = fopen("cons_obey.txt","w");

	 for(i=0;i<300;i++)
	 {

		 //==============================DSP读串口===============================
		 //发送读数据指令
		 t_f2 = post_write_uart(&token_read,1);
		 //读8字节数据
		 t_f1 = post_read_uart(&Buffer_read[0],16,delay);

		 //数据解封装，化为float型
		 memcpy(&count_r_verify,&Buffer_read[0],sizeof(float));
		 memcpy(&y1_float,&Buffer_read[4],sizeof(float));
		 memcpy(&y2_float,&Buffer_read[8],sizeof(float));
		 //校验
		 if(count_r_verify!=count_read)
		 {
		    printf("Count_read Verify failed! \n\r");
		    printf("Is read successful? %d \n\r",t_f1);
		    printf("What's in the buffer? %c \n\r",Buffer_read);
		    break;
	     }

	    count_read = count_read + 1;

	    //精度转换
	    y1_double = (double)y1_float;
	    y2_double = (double)y2_float;
	    //printf("y1 is %f\n",y1_double);
	    //printf("y2 is %f\n",y2_double);

	    y_k[0] = y1_double;y_k[1] = y2_double;

	    //==============================DSP读串口结束===============================

		 MPC_start_time = Timestamp_get32();

		 if((i<125)||(i>250))
			 mpc_V3_DP(delta_u,delta_u_ini,y_ini,lambda_ini,u_k,y_k,x_k,r_k1,u_k_1,u_k_2,i);
		 else
			 mpc_V3_DP(delta_u,delta_u_ini,y_ini,lambda_ini,u_k,y_k,x_k,r_k2,u_k_1,u_k_2,i);


		 MPC_end_time = Timestamp_get32();
		 MPC_total_time = MPC_end_time - MPC_start_time - dead_time;
		 Time_MPC = (double)(MPC_total_time);
		 printf("%f\n",Time_MPC/freq.lo*1000);
		 fprintf(fp_mpct,"%f\n",Time_MPC/freq.lo*1000);Time_MPC_Total += Time_MPC/freq.lo*1000;
		 if (Time_MPC/freq.lo*1000 > Time_MPC_Worst) Time_MPC_Worst = Time_MPC/freq.lo*1000;

		 //Here the DSP communicate with Matlab to transmite u(k) and get y(k)

		 //==============================DSP写串口==================================

		 //精度转换
		 u_float = (float)u_k[0];
		 //printf("u is %f\n",u_float);
		 //封装u
		 memcpy(&Buffer_write[0],&count_write,sizeof(float));
		 memcpy(&Buffer_write[4],&u_float,sizeof(float));
		 //发送写数据指令
		 t_f2 = post_write_uart(&token_write,1);
		 //写8位数据
		 t_f2 = post_write_uart(&Buffer_write[0],8);
		 count_write = count_write + 1;
		 //==============================DSP写串口结束===============================

		 ////Here we use a simple model instead.
		 //process_model2(xm,u_k,y_k);

		 memcpy(u_k_2,u_k_1,nu*sizeof(double));
		 memcpy(u_k_1,u_k,nu*sizeof(double));

		 //printf("Manipulated Variables	(u(k)):\n");
		 //show_matrix(u_k,1,nu);putchar('\n');
		 fprintf(fp_u,"%f\n",u_k[0]);
		 //printf("Manipulated Variables	(u(k-1)):\n");
		 //show_matrix(u_k_1,1,nu);putchar('\n');
		 //printf("Manipulated Variables	(u(k-2)):\n");
		 //show_matrix(u_k_2,1,nu);putchar('\n');
		 //printf("Aug State Variables	(x(k)):\n");
		 //show_matrix(x_k,1,nx+ny);putchar('\n');
		 //printf("State Variables		(xm(k)):\n");
		 //show_matrix(xm,1,nx);putchar('\n');
		 //printf("Controlled Variables	(y(k)):\n");
		 //show_matrix(y_k,1,ny);putchar('\n');putchar('\n');
		 fprintf(fp_y1,"%f\n",y_k[0]);fprintf(fp_y2,"%f\n",y_k[1]);fprintf(fp_y3,"%f\n",y_k[2]);
		 if(i<100) {
			 IAE += scalar_abs(y_k[0] - r_k1[0])+scalar_abs( y_k[1] - r_k1[1]);
			 ITAE += i * 0.5 * (scalar_abs(y_k[0] - r_k1[0])+scalar_abs( y_k[1] - r_k1[1]));
			 ISE += (y_k[0] - r_k1[0])*(y_k[0] - r_k1[0])+(y_k[1] - r_k1[1])*(y_k[1] - r_k1[1]);
		 }
		 else {
			 IAE += scalar_abs(y_k[0] - r_k1[0])+scalar_abs( y_k[1] - r_k1[1]);
			 ITAE += i * 0.5 * (scalar_abs(y_k[0] - r_k1[0])+scalar_abs( y_k[1] - r_k1[1]));
			 ISE += (y_k[0] - r_k1[0])*(y_k[0] - r_k1[0])+(y_k[1] - r_k1[1])*(y_k[1] - r_k1[1]);
		 }

		 fprintf(fp_t,"%f\n",Time_QP/freq.lo*1000);
		 //printf("QP_Iter:%d\n",QP_Iteration);

		 total_QP_iterations = QP_Iteration + total_QP_iterations;
		 if ( QP_Iteration > worst_QP_iterations) worst_QP_iterations =  QP_Iteration;

		 fprintf(fp_i,"%d\n",QP_Iteration);
		 total_QP_iterations = QP_Iteration + total_QP_iterations;
		 //fprintf(fp_qpt_nonIter,"%f\n",Time_QP_non_Iter/freq.lo*1000);
		 //fprintf(fp_qpt_Iter,"%f\n",Time_QP_Iter/freq.lo*1000);
		 //fprintf(fp_s1t,"%f\n",Time_Section1/freq.lo*1000);fprintf(fp_s2t,"%f\n",Time_Section2/freq.lo*1000);
		 //fprintf(fp_s3t,"%f\n",Time_Section3/freq.lo*1000);fprintf(fp_s4t,"%f\n",Time_Section4/freq.lo*1000);
		 //fprintf(fp_s5t,"%f\n",Time_Section5/freq.lo*1000);
		 //fprintf(fp_Fx,"Control Step: %d\n\n",i);
		 //fprintf(fp_cons_obey,"Control Step: %d\n\n",i);

		 Time_QP = 0;QP_Iteration = 0;Time_QP_non_Iter = 0;Time_QP_Iter = 0;
		 Time_Section1 = 0;Time_Section2 = 0;Time_Section3 = 0;Time_Section4 = 0;Time_Section5 = 0;
	 }

	 fclose(fp_u);fclose(fp_y1);fclose(fp_y2);fclose(fp_y3);fclose(fp_t);fclose(fp_i);
	 fclose(fp_qpt_nonIter);fclose(fp_qpt_Iter);
	 fclose(fp_s1t);fclose(fp_s2t);fclose(fp_s3t);fclose(fp_s4t);fclose(fp_s5t);
	 fclose(fp_Fx);fclose(fp_cons_obey);

	 //end_time = clock();
	 //total_time = end_time - start_time - dead_time;
	 //Time_MPC = (double)total_time;
	 //printf("%d times duration: %f\n",i,(double)total_time);

	 //printf("MPC time consuming: %f\n",Time_MPC);
	 //printf("QP solving time consuming(Clocks): %f\n",Time_QP/freq.lo);
	 //printf("QP solving time consuming(ms Per Single time): %f\n",((Time_QP/freq.lo)/i)*1000);
	 //printf("QP solving step 1 time consuming: %f\n",Time_QP_step1);
	 //printf("QP solving step 2 time consuming: %f\n",Time_QP_step2);
	 //printf("QP solving step 3 time consuming: %f\n",Time_QP_step3);
	 //printf("QP solving step 4 time consuming: %f\n",Time_QP_step4);
	 //printf("Linear system1 solving time consuming(ms): %f\n",(Time_LineSolve1/freq.lo)*1000);
	 //printf("Linear system1 solving times: %d\n",LS_solve_times);
	 //printf("Linear system1 average solving time(ms): %f\n",(Time_LineSolve1/freq.lo)*1000/LS_solve_times);
	 //printf("Linear system2 solving time cnsuming: %f\n\n",Time_LineSolve2);

	 //printf("QP Section 1 time consuming: %f\n",Time_Section1);
	 //printf("QP Section 2 time consuming: %f\n",Time_Section2);
	 //printf("QP Section 3 time consuming(ms): %f\n",(Time_Section3/freq.lo)*1000);
	 //printf("alpha_compute_times:%d\n",alpha_compute_times);
	 //printf("QP Section 4 time consuming: %f\n",Time_Section4);
	 //printf("QP Section 5 time consuming: %f\n",Time_Section5);
	 //printf("QP Section 6 time consuming: %f\n",Time_Section6);
	 //printf("QP Section 7 time consuming: %f\n",Time_Section7);
	 //printf("QP Section 8 time consuming: %f\n\n",(Time_Section8/freq.lo)*1000);
	 //printf("Mat_vec solving times: %d\n",Mat_vec_times);
	 //printf("Mat_vec_times average solving time(ms): %f\n",(Time_Section8/freq.lo)*1000/Mat_vec_times);

	 //printf("dot_product time consuming (ms): %f\n",(Time_bsxfun/freq.lo)*1000);
	 //printf("dot_product solving times: %d\n",dot_product_times);
	 //printf("dot_product average solving time(ms): %f\n",(Time_dot_product/freq.lo)*1000/dot_product_times);

	 //printf("vec_add time consuming (ms): %f\n",(Time_vec_add/freq.lo)*1000);
	 //printf("vec_add solving times: %d\n",vec_add_times);
	 //printf("vec_add average solving time(ms): %f\n\n",(Time_vec_add/freq.lo)*1000/vec_add_times);
	 //printf("vec_sub time consuming (ms): %f\n",(Time_vec_sub/freq.lo)*1000);
	 //printf("vec_sub solving times: %d\n",vec_sub_times);
	 //printf("vec_sub average solving time(ms): %f\n\n",(Time_vec_sub/freq.lo)*1000/vec_sub_times);

	 //printf("vec_div time consuming (ms): %f\n",(Time_vec_div/freq.lo)*1000);
	 //printf("vec_div solving times: %d\n",vec_div_times);
	 //printf("vec_div average solving time(ms): %f\n\n",(Time_vec_div/freq.lo)*1000/vec_div_times);
	 //printf("vec_mul time consuming (ms): %f\n",(Time_vec_mul/freq.lo)*1000);
	 //printf("vec_mul solving times: %d\n",vec_mul_times);
	 //printf("vec_mul average solving time(ms): %f\n\n",(Time_vec_mul/freq.lo)*1000/vec_mul_times);

	 //printf("bsxfun time consuming (ms): %f\n",(Time_bsxfun/freq.lo)*1000);
	 //printf("vec_div solving times: %d\n",bsxfun_times);
	 //printf("vec_div average solving time(ms): %f\n\n",(Time_bsxfun/freq.lo)*1000/bsxfun_times);

	 //printf("Mat_mul time consuming (ms): %f\n",(Time_Mat_mul/freq.lo)*1000);
	 //printf("Mat_mul solving times: %d\n",mat_mul_times);
	 //printf("Mat_mul average solving time(ms): %f\n\n",(Time_Mat_mul/freq.lo)*1000/mat_mul_times);

	 //printf("sca2vec time consuming (ms): %f\n",(Time_sca2vec/freq.lo)*1000);
	 //printf("sca2vec solving times: %d\n",sca2vec_times);
	 //printf("sca2vec average solving time(ms): %f\n\n",(Time_sca2vec/freq.lo)*1000/sca2vec_times);

	 //printf("vec_rev time consuming (ms): %f\n",(Time_vec_rev/freq.lo)*1000);
	 //printf("vec_rev solving times: %d\n",vec_rev_times);
	 //printf("vec_rev average solving time(ms): %f\n\n",(Time_vec_rev/freq.lo)*1000/vec_rev_times);

	 //printf("sca_vec_mut time consuming (ms): %f\n",(Time_sca_vec_mut/freq.lo)*1000);
	 //printf("sca_vec_mut solving times: %d\n",sca_vec_mut_times);
	 //printf("sca_vec_mut average solving time(ms): %f\n\n",(Time_sca_vec_mut/freq.lo)*1000/sca_vec_mut_times);

	 // Memory free
	 // main function memory free
	 free(L);free(F);free(Phi);free(B);free(OMEGA_L);free(N_OMEGA_L);free(N_OMEGA_L_T);
	 free(G);free(GL);free(GL_T);free(delta_u_ini);free(y_ini);free(lambda_ini);
	 free(F1);free(Phi1);free(QQ);free(Q_temp);
	 //free(F1);free(Phi1);

	 // MPC function memory free
	 free(aug_u_k_1);free(omega_r);free(n_omega_r);free(c);free(c_rev);//free(h_x);
	 free(MPC_temp1);free(MPC_temp2);free(MPC_temp3);free(MPC_temp4);

	 // QP function memory free
	 free(x_old);free(rd);free(rp);free(rp_y);
	 //free(y_diag);free(inv_y_diag);free(lambda_diag);free(inv_lambda_diag);
	 free(inv_y_lambda);free(y_lambda_A);
	 free(equation_F);free(equation_b);//free(delta_x_aff);free(delta_y_aff);free(delta_lambda_aff);free(diag);
	 free(delta_x);free(delta_y);free(delta_lambda);free(temp_cond);
	 free(x_try);free(sig_mu);free(center_part);
	 free(QP_temp1);free(QP_temp2);free(QP_temp3);

	 // Priduip function memory free
	 //free(delta_x_aff);free(delta_y_aff);free(delta_lambda_aff);
	 //free(sig_mu_lam_y);free(aff_part);

	 // feedback function memory free
	 free(fb_temp1);free(fb_temp2);

	 // mcf function memory free
	 free(mcf_C);free(mcf_temp_vec);

	 // alpha_decreas function memory free
	 free(alpha_decreas_neg_delta_value);free(alpha_decreas_neg_vec_value);
	 free(alpha_decreas_temp1);free(alpha_decreas_temp2);

	 // line_solve function memory free
	 free(line_solve_L);free(line_solve_p);
	 free(equation_F_scaling);free(equation_b_scaling);
	 free(line_solve_U);free(line_solve_D);free(line_solve_LD);

	 //CDC function memory free
	 free(CDC_temp1);free(CDC_temp2);

	 //printf("Time_MPC_Total: %f\n",Time_MPC_Total);Time_MPC_Total = 0;
	 //printf("total_QP_iterations: %d\n",total_QP_iterations);total_QP_iterations = 0;

     printf("Time_MPC_Average: %f\n",Time_MPC_Total/i);Time_MPC_Total = 0;
     printf("Time_MPC_Worst: %f\n",Time_MPC_Worst);Time_MPC_Worst = 0;
     printf("average_QP_iterations: %f\n",(double)total_QP_iterations/i);total_QP_iterations = 0;
     printf("worst_QP_iterations: %d\n",worst_QP_iterations);worst_QP_iterations = 0;
     printf("IAE: %f\n",IAE);	IAE = 0;
     printf("ITAE: %f\n",ITAE);	ITAE = 0;
     printf("ISE: %f\n",ISE);	ISE = 0;

	 printf("exit taskFxn()\n");
 }

 Void main()
 {
 //    Task_Handle task;
     Error_Block eb;

     System_printf("enter main()\n");

     Error_init(&eb);
 //    task = Task_create(taskFxn, NULL, &eb);
 //    if (task == NULL) {
 //        System_printf("Task_create() failed!\n");
 //        BIOS_exit(0);
 //    }

     BIOS_start();     /* enable interrupts and start SYS/BIOS */
 }
