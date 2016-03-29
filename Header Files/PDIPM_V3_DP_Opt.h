/*
 * PDIPM_V3_DP_Opt.h
 *
 *  Created on: 2014-8-18
 *      Author: WORK
 */

#ifndef PDIPM_V3_DP_OPT_H_
#define PDIPM_V3_DP_OPT_H_


//V4 modification:
//	1 Himmelblau Termination Criteria
//	2 Initial Value Reset
void priduip_v4(double *G, double *GL,double *GL_T, double *c, double *A, double *A_t, double *b, double *x, double *y, double *lambda, const int m, const int n, double *delta_u);



#endif /* PDIPM_V3_DP_OPT_H_ */
