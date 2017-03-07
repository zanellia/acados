/*
 *    This file is part of ACADOS.
 *
 *    ACADOS is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    ACADOS is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with ACADOS; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#if defined(__APPLE__)
#include <mach/mach_time.h>
#elsev
#include <sys/stat.h>
#endif

#define _GNU_SOURCE

// system headers
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

// HPMPC headers
#include "hpmpc/include/aux_d.h"

// ACADOS headers
#include "acados/utils/types.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_hpmpc.h"
#include "acados/utils/tools.h"

#include "acados_kitty_interface.h"
#include "../models/aircraft_integrator.h"
#include "../ocps/residual_x_eval_wrapper.h"
#include "../ocps/residual_u_eval_wrapper.h"


#include <fenv.h>
#include "eispack.h"
#include <assert.h>

#include <plotting.h>

// extern int feenableexcept (unsigned int excepts);
// flush denormals to zero
#if defined(TARGET_X64_AVX2) || defined(TARGET_X64_AVX) ||  \
    defined(TARGET_X64_SSE3) || defined(TARGET_X86_ATOM) || \
    defined(TARGET_AMD_SSE3)
#include <xmmintrin.h>  // needed to flush to zero sub-normals with _MM_SET_FLUSH_ZERO_MODE (_MM_FLUSH_ZERO_ON); in the main()
#endif

// define IP solver arguments && number of repetitions
#define NREP 1000
#define MAX_IP_ITER 50
#define TOL 1e-8
#define MINSTEP 1e-8

#define NP_dummy 43

#define STEP_SIZE 0.05   // shooting interval size
#define nn 10
#define RKSTEPS 10 // intermediate rk4 steps
#define NSIM 100


static void print_states_controls(real_t *w, real_t T, int_t N, int_t NX, int_t NU) {
    FILE *fp = fopen("acados_log.txt","w+");
    printf("node\tx\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\tu\t\t\t\t\n");
    for (int_t i = 0; i < N; i++) {
        printf("%4d\t%+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e%+e %+e\n", i, w[i*(NX+NU)], w[i*(NX+NU)+1],
        w[i*(NX+NU)+2], w[i*(NX+NU)+3], w[i*(NX+NU)+4], w[i*(NX+NU)+5] , w[i*(NX+NU)+6] , w[i*(NX+NU)+7],
        w[i*(NX+NU)+8], w[i*(NX+NU)+9], w[i*(NX+NU)+10], w[i*(NX+NU)+11], w[i*(NX+NU)+12], w[i*(NX+NU)+13], w[i*(NX+NU)+14]);

        // print to log file too
        fprintf(fp,"%e,\t%+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e,%+e, %+e\n", i*T/N, w[i*(NX+NU)], w[i*(NX+NU)+1],
        w[i*(NX+NU)+2], w[i*(NX+NU)+3], w[i*(NX+NU)+4], w[i*(NX+NU)+5] , w[i*(NX+NU)+6] , w[i*(NX+NU)+7],
        w[i*(NX+NU)+8], w[i*(NX+NU)+9], w[i*(NX+NU)+10], w[i*(NX+NU)+11], w[i*(NX+NU)+12], w[i*(NX+NU)+13], w[i*(NX+NU)+14]);
    }

    printf("%4d\t%+e %+e %+e %+e %+e %+e %+e %+e %+e %+e %+e\n", N, w[N*(NX+NU)], w[N*(NX+NU)+1],
    w[N*(NX+NU)+2], w[N*(NX+NU)+3], w[N*(NX+NU)+4], w[N*(NX+NU)+5] , w[N*(NX+NU)+6] , w[N*(NX+NU)+7],
    w[N*(NX+NU)+8], w[N*(NX+NU)+9], w[N*(NX+NU)+10]);

    // print to log file too
    fprintf(fp, "%e,\t%+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e, %+e,%+e, %+e\n", N*T/N, w[N*(NX+NU)], w[N*(NX+NU)+1],
    w[N*(NX+NU)+2], w[N*(NX+NU)+3], w[N*(NX+NU)+4], w[N*(NX+NU)+5] , w[N*(NX+NU)+6] , w[N*(NX+NU)+7],
    w[N*(NX+NU)+8], w[N*(NX+NU)+9], w[N*(NX+NU)]+10, 0.0, 0.0, 0.0, 0.0);
    fclose(fp);
}


// static void shift_states(real_t *w, real_t *x_end, int_t N) {
//     for (int_t i = 0; i < N; i++) {
//         for (int_t j = 0; j < NX; j++) w[i*(NX+NU)+j] = w[(i+1)*(NX+NU)+j];
//     }
//     for (int_t j = 0; j < NX; j++) w[N*(NX+NU)+j] = x_end[j];
// }

// static void shift_controls(real_t *w, real_t *u_end, int_t N) {
//     for (int_t i = 0; i < N-1; i++) {
//         for (int_t j = 0; j < NU; j++) w[i*(NX+NU)+NX+j] = w[(i+1)*(NX+NU)+NX+j];
//     }
//     for (int_t j = 0; j < NU; j++) w[(N-1)*(NX+NU)+NX+j] = u_end[j];
// }

// TODO(Andrea): need to fix this stuff below...
extern int d_ip2_res_mpc_hard_work_space_size_bytes_libstr(int N, int *nx,
  int *nu, int *nb, int *ng);

extern  int d_size_strmat(int m, int n);
extern  int d_size_strvec(int m);

void init_acados(nmpc_data* nmpc_data, rk4_int* rk4_int, init_nmpc_data* init_data, acados_options *acados_options){

    // get dimensions
    const int NN = init_data->NN;   // number of stages
    const int NX = init_data->NX;   // number of states
    const int NU = init_data->NU;   // number of inputs
    const int NP = init_data->NP;   // numbers of parameters
    const int NR = init_data->NR;   // number of residuals (least-square)

    const int NB0 = init_data->NB0; // number of input bounds on stage 0
    const int NB = init_data->NB;   // number of input bounds on stage 1 to N-1
    const int NBN = init_data->NBN; // number of input bounds on stage N

    // and copy them into nmpc_data (apart from lb0, lb and lbN, we'll need some mallocs...)
    nmpc_data->NN = init_data->NN;   // number of stages
    nmpc_data->NX = init_data->NX;   // number of states
    nmpc_data->NU = init_data->NU;   // number of inputs
    nmpc_data->NP = init_data->NP;   // numbers of parameters
    nmpc_data->NR = init_data->NR;   // number of residuals (least-square)

    nmpc_data->NB0 = init_data->NB0; // number of input bounds on stage 0
    nmpc_data->NB = init_data->NB;   // number of input bounds on stage 1 to N-1
    nmpc_data->NBN = init_data->NBN; // number of input bounds on stage N

    // allocate memory for rk4_int
    real_t h  = 0.1;
    real_t t0 = 0.0;

    // rk4_int = malloc(sizeof *rk4_int);

    rk4_int->x_in = malloc(sizeof(*rk4_int->x_in) * NX);
    rk4_int->u_in = malloc(sizeof(*rk4_int->u_in) * NU);
    rk4_int->p_in = malloc(sizeof(*rk4_int->p_in) * NP);

    rk4_int->x_out = malloc(sizeof(*rk4_int->x_out) * (NX));
    rk4_int->Sx_out = malloc(sizeof(*rk4_int->Sx_out) * (NX*NX));
    rk4_int->Su_out = malloc(sizeof(*rk4_int->Su_out) * (NX*NU));

    rk4_int->h_in = h;
    rk4_int->t0_in = t0;
    rk4_int->nx = NX;
    rk4_int->nu = NU;
    rk4_int->np = NP;
    rk4_int->n_steps = 1;

    // allocate internal buffers for rk4 integrator
    rk4_int->intermediate_Sx_in = malloc(sizeof(rk4_int->intermediate_Sx_in)*NX*NX);
    rk4_int->intermediate_Sx_out = malloc(sizeof(rk4_int->intermediate_Sx_out)*NX*NX);
    rk4_int->intermediate_x_in = malloc(sizeof(rk4_int->intermediate_x_in)*NX);
    rk4_int->intermediate_x_out = malloc(sizeof(rk4_int->intermediate_x_out)*NX);
    rk4_int->intermediate_Su = malloc(sizeof(rk4_int->intermediate_Su)*NX*NU);
    rk4_int->temp_Sx_out = malloc(sizeof(rk4_int->temp_Sx_out)*NX*NX);
    rk4_int->temp_Su_out = malloc(sizeof(rk4_int->temp_Su_out)*NX*NX);

    int  sz_arg;
    int  sz_res;
    int  sz_iw;
    int  sz_w;
    rk4_int->eval_dynamics_work(&sz_arg, &sz_res, &sz_iw, &sz_w);
    rk4_int->internal_mem = malloc(sizeof(real_t)*(sz_w));

    // allocate memory for qp_in
    int_t *nx;
    int_zeros(&nx, NN+1, 1);
    int_t *nu;
    int_zeros(&nu, NN+1, 1);
    int_t *nb;
    int_zeros(&nb, NN+1, 1);
    int_t *nc;
    int_zeros(&nc, NN+1, 1);
    for (int_t i = 0; i < NN; i++) {
        nx[i] = NX;
        nu[i] = NU;
    }
    nx[0]  = 11;  // x0 is NOT eliminated
    nx[NN] = NX;
    nu[NN] = 0;

    real_t **pA;
    pA = malloc(NN*sizeof(real_t*));
    real_t **pB;
    pB = malloc(NN*sizeof(real_t*));
    real_t **pb;
    pb = malloc(NN*sizeof(real_t*));
    real_t **pQ;
    pQ = malloc((NN+1)*sizeof(real_t*));
    real_t **pS;
    pS = malloc(NN*sizeof(real_t*));
    real_t **pR;
    pR = malloc(NN*sizeof(real_t*));
    real_t **pq;
    pq = malloc((NN+1)*sizeof(real_t*));
    real_t **pr;
    pr = malloc(NN*sizeof(real_t*));
    real_t **px;
    px = malloc((NN+1)*sizeof(real_t*));
    real_t **pu;
    pu = malloc(NN*sizeof(real_t*));
    // real_t **px0;
    // px0 = malloc(sizeof(real_t*));
    int    **hidxb;
    hidxb = malloc((NN+1)*sizeof(int*));
    real_t **pub;
    pub = malloc((NN+1)*sizeof(real_t*));
    real_t **plb;
    plb = malloc((NN+1)*sizeof(real_t*));
    real_t **pC;
    pC = malloc((NN+1)*sizeof(real_t*));
    real_t **pD;
    pD = malloc(NN*sizeof(real_t*));
    real_t **plg;
    plg = malloc((NN+1)*sizeof(real_t*));
    real_t **pug;
    pug = malloc((NN+1)*sizeof(real_t*));

    // box constraints
    int ii, jj;
    nb[0] = NB0;
    for (ii = 1; ii < NN; ii++ ) nb[ii] = NB;
    nb[NN] = NBN;

    int *idxb0;  // TODO(Andrea): need to swap these guys in the interface...
    int_zeros(&idxb0, NB0, 1);
    for (jj = 0; jj < NB0; jj++) idxb0[jj] = init_data->idxb0[jj];

    int *idxb;
    int_zeros(&idxb, NB, 1);
    for (jj = 0; jj < NB; jj++ ) idxb[jj] = init_data->idxb[jj];

    int *idxbN;
    int_zeros(&idxbN, NBN, 1);
    for ( jj = 0; jj < NBN; jj++ ) idxbN[jj] = init_data->idxbN[jj];

    hidxb[0] = idxb0;
    for (ii = 1; ii < NN; ii++ ) hidxb[ii] = idxb;
    hidxb[NN] = idxbN;

    d_zeros(&plb[0], NB0, 1);
    d_zeros(&pub[0], NB0, 1);
    d_zeros(&pA[0], nx[0+1], nx[0]);
    d_zeros(&pB[0], nx[0+1], nu[0]);
    d_zeros(&pb[0], nx[0+1], 1);
    d_zeros(&pQ[0], nx[0], nx[0]);
    d_zeros(&pR[0], nu[0], nu[0]);
    d_zeros(&pS[0], nu[0], nx[0]);
    d_zeros(&pq[0], nx[0], 1);
    d_zeros(&pr[0], nu[0], 1);
    d_zeros(&px[0], nx[0], 1);
    d_zeros(&pu[0], nu[0], 1);
    for (int_t i = 1; i < NN; i++) {
        d_zeros(&pA[i], nx[i+1], nx[i]);
        d_zeros(&pB[i], nx[i+1], nu[i]);
        d_zeros(&pb[i], nx[i+1], 1);
        d_zeros(&pQ[i], nx[i], nx[i]);
        d_zeros(&pR[i], nu[i], nu[i]);
        d_zeros(&pS[i], nu[i], nx[i]);
        d_zeros(&pq[i], nx[i], 1);
        d_zeros(&pr[i], nu[i], 1);
        d_zeros(&px[i], nx[i], 1);
        d_zeros(&pu[i], nu[i], 1);
        d_zeros(&plb[i], NB, 1);
        d_zeros(&pub[i], NB, 1);
    }
    d_zeros(&plb[NN], NBN, 1);
    d_zeros(&pub[NN], NBN, 1);
    d_zeros(&pq[NN], nx[NN], 1);
    d_zeros(&px[NN], nx[NN], 1);
    d_zeros(&pQ[NN], nx[NN], nx[NN]);

    int ng = 0;  // number of general constraints
    int ngN = 0; // number of general constraints at the last stage

    int ngg[NN + 1];
    for (ii = 0; ii < NN; ii++) ngg[ii] = ng;
    ngg[NN] = ngN;

    // general constraints
    double *C0;
    d_zeros(&C0, 0, NX);
    double *D0;
    d_zeros(&D0, 0, NU);
    double *lg0;
    d_zeros(&lg0, 0, 1);
    double *ug0;
    d_zeros(&ug0, 0, 1);

    double *C;
    d_zeros(&C, 0, NX);
    double *D;
    d_zeros(&D, 0, NU);
    double *lg;
    d_zeros(&lg, 0, 1);
    double *ug;
    d_zeros(&ug, 0, 1);

    double *CN;
    d_zeros(&CN, ngN, NX);
    double *lgN;
    d_zeros(&lgN, ngN, 1);  // force all states to 0 at the last stage
    double *ugN;
    d_zeros(&ugN, ngN, 1);  // force all states to 0 at the last stage

    pC[0] = C0;
    pD[0] = D0;
    plg[0] = lg0;
    pug[0] = ug0;
    real_t **ppi;
    ppi = malloc(sizeof(real_t *)*NN);
    real_t **plam;
    plam = malloc(sizeof(real_t *)*(NN+1));

    ii = 0;
    d_zeros(&ppi[ii], nx[ii+1], 1);
    d_zeros(&plam[ii], 2*nb[ii]+2*nb[ii], 1);
    for (ii = 1; ii < NN; ii++) {
        pC[ii] = C;
        pD[ii] = D;
        plg[ii] = lg;
        pug[ii] = ug;
        d_zeros(&ppi[ii], nx[ii+1], 1);
        d_zeros(&plam[ii], 2*nb[ii]+2*nb[ii], 1);
    }

    d_zeros(&plam[NN], 2*nb[NN]+2*nb[NN], 1);

    pC[NN] = CN;
    plg[NN] = lgN;
    pug[NN] = ugN;

    // solver arguments
    ocp_qp_hpmpc_args *hpmpc_args;
    hpmpc_args = malloc(sizeof *hpmpc_args);
    hpmpc_args->tol = TOL;
    hpmpc_args->max_iter = MAX_IP_ITER;
//  hpmpc_args.min_step = MINSTEP;
    hpmpc_args->mu0 = 10;
//  hpmpc_args.sigma_min = 1e-3;
    hpmpc_args->warm_start = 0;
    hpmpc_args->N2 = NN;

    // work space

    int work_space_size = d_ip2_res_mpc_hard_work_space_size_bytes_libstr(NN,
      nx, nu, nb, ngg);

    // Adding memory for data
    for ( int ii=0; ii <NN; ii++ ) {
        work_space_size+= d_size_strmat(nu[ii]+nx[ii]+1, nx[ii+1]);
        work_space_size+= d_size_strvec(nx[ii+1]);
        work_space_size+= d_size_strmat(nu[ii]+nx[ii]+1, nu[ii]+nx[ii]);
        work_space_size+= d_size_strvec(nu[ii]+nx[ii]);
        work_space_size+= d_size_strmat(nu[ii]+nx[ii]+1, ngg[ii]);
        work_space_size+= d_size_strvec(2*nb[ii]+2*ngg[ii]);
        work_space_size+= d_size_strvec(nu[ii]+nx[ii]);
        work_space_size+= d_size_strvec(nx[ii+1]);
        work_space_size+= d_size_strvec(2*nb[ii]+2*ngg[ii]);
        work_space_size+= d_size_strvec(2*nb[ii]+2*ngg[ii]);
    }


    work_space_size+= d_size_strmat(nu[NN]+nx[NN]+1, nu[NN]+nx[NN]);
    work_space_size+= d_size_strvec(2*nb[NN]+2*ngg[NN]);
    work_space_size+= d_size_strvec(nu[NN]+nx[NN]);
    work_space_size+= d_size_strvec(nx[NN]);
    work_space_size+= d_size_strvec(2*nb[NN]+2*ngg[NN]);
    work_space_size+= d_size_strvec(2*nb[NN]+2*ngg[NN]);

    work_space_size += sizeof(int_t)*MAX_IP_ITER*5;

    // work_space_size += 1*sizeof(double)*(NN+1); // TODO(Andrea): ??

    void *workspace;

    v_zeros_align(&workspace, work_space_size);

    // allocate OCP QP variables
    ocp_qp_in *qp_in;
    qp_in = malloc(sizeof *qp_in);
    qp_in->N = NN;
    ocp_qp_out *qp_out;
    qp_out = malloc(sizeof *qp_out);
    qp_in->nx = nx;
    qp_in->nu = nu;
    qp_in->nb = nb;
    qp_in->nc = nc;

    qp_in->Q = (const real_t **) pQ;
    qp_in->S = (const real_t **) pS;
    qp_in->R = (const real_t **) pR;
    qp_in->q = (const real_t **) pq;
    qp_in->r = (const real_t **) pr;
    qp_in->A = (const real_t **) pA;
    qp_in->B = (const real_t **) pB;
    qp_in->b = (const real_t **) pb;
    qp_in->lb = (const real_t **) plb;
    qp_in->ub = (const real_t **) pub;
    qp_in->idxb = (const int_t **) hidxb;
    qp_in->Cx = (const real_t **) pC;
    qp_in->Cu = (const real_t **) pD;
    qp_in->lc = (const real_t **) plg;
    qp_in->uc = (const real_t **) pug;

    qp_out->x = px;
    qp_out->u = pu;
    qp_out->pi = ppi;
    qp_out->lam = plam;

    // allocate memory for nmpc_data
    // Problem data
    real_t  *w; // States and controls stacked
    d_zeros(&w, NN*(NX+NU)+NX, 1);

    real_t *lb0 = malloc(sizeof(real_t)*NB0);
    real_t *lb = malloc(sizeof(real_t)*NB);
    real_t *lbN = malloc(sizeof(real_t)*NBN);

    real_t *ub0 = malloc(sizeof(real_t)*NB0);
    real_t *ub = malloc(sizeof(real_t)*NB);
    real_t *ubN = malloc(sizeof(real_t)*NBN);

    nmpc_data->lb0 = lb0;
    nmpc_data->lb = lb;
    nmpc_data->lbN = lbN;

    nmpc_data->ub0 = ub0;
    nmpc_data->ub = ub;
    nmpc_data->ubN = ubN;

    nmpc_data->idxb0 = idxb0;
    nmpc_data->idxb = idxb;
    nmpc_data->idxbN = idxbN;

    int_t   max_sqp_iters = 1;
    int_t   max_iters     = 1;

    nmpc_data->w = w;
    nmpc_data->max_sqp_iters = max_sqp_iters;
    nmpc_data->max_iters = max_iters;

    nmpc_data->qp_in = qp_in;
    nmpc_data->qp_out = qp_out;
    nmpc_data->hpmpc_args = hpmpc_args;
    nmpc_data->workspace = workspace;

    // allocate memory for residual evalutation
    // get internal memory dimension
    nmpc_data->res_x_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

    d_zeros((double **)&nmpc_data->residual_x_eval_mem, sz_w, 1);
    d_zeros(&nmpc_data->residual_x_out, NR + NR*NX + NR, 1 );
    d_zeros(&nmpc_data->residual_x_in, NR + NP + 1 +1, 1);
    d_zeros(&nmpc_data->drdx_tran, NR*NX, 1 );

    nmpc_data->res_u_work(&sz_arg, &sz_res, &sz_iw, &sz_w);

    d_zeros((double **)&nmpc_data->residual_u_eval_mem, sz_w, 1);
    d_zeros(&nmpc_data->residual_u_out, NU + NU*NU + NU, 1 );
    d_zeros(&nmpc_data->residual_u_in, NU + NP + 1 +1, 1 );
    d_zeros(&nmpc_data->drdu_tran, NU*NU, 1 );

    acados_options->non_symmetry_tol = 1e-4;
    acados_options->print_level = 1;

    if (acados_options->use_gnuplot){
      nmpc_data->gnuplotPipe = popen("gnuplot -persist", "w");
    }

    return;
  }

  int_t de_init_acados(nmpc_data* nmpc_data, rk4_int* rk4_int){

      const int NN = nmpc_data->NN;   // number of stages

      free(rk4_int->x_in);
      free(rk4_int->u_in);
      free(rk4_int->p_in);

      free(rk4_int->x_out);
      free(rk4_int->Sx_out);
      free(rk4_int->Su_out);

      free(rk4_int->intermediate_Sx_in);
      free(rk4_int->intermediate_Sx_out);
      free(rk4_int->intermediate_x_in);
      free(rk4_int->intermediate_x_out);
      free(rk4_int->intermediate_Su);
      free(rk4_int->temp_Sx_out);
      free(rk4_int->temp_Su_out);

      free(rk4_int->internal_mem);

      // free memory for qp_in
      free((void *)nmpc_data->qp_in->nx);
      free((void *)nmpc_data->qp_in->nu);
      free((void *)nmpc_data->qp_in->nb);
      free((void *)nmpc_data->qp_in->nc);
      // free(nmpc_data->qp_in->idxb[0]);
      // free(nmpc_data->qp_in->idxb[1]);
      // free(nmpc_data->qp_in->idxb[NN]);

      free((void *)nmpc_data->qp_in->lb[0]);
      free((void *)nmpc_data->qp_in->ub[0]);
      free((void *)nmpc_data->qp_in->A[0]);
      free((void *)nmpc_data->qp_in->B[0]);
      free((void *)nmpc_data->qp_in->b[0]);
      free((void *)nmpc_data->qp_in->Q[0]);
      free((void *)nmpc_data->qp_in->R[0]);
      free((void *)nmpc_data->qp_in->S[0]);
      free((void *)nmpc_data->qp_in->q[0]);
      free((void *)nmpc_data->qp_in->r[0]);

      free(nmpc_data->qp_out->x[0]);
      free(nmpc_data->qp_out->u[0]);

      for (int_t i = 1; i < NN; i++) {
          free((void *)nmpc_data->qp_in->A[i]);
          free((void *)nmpc_data->qp_in->B[i]);
          free((void *)nmpc_data->qp_in->b[i]);
          free((void *)nmpc_data->qp_in->Q[i]);
          free((void *)nmpc_data->qp_in->R[i]);
          free((void *)nmpc_data->qp_in->S[i]);
          free((void *)nmpc_data->qp_in->q[i]);
          free((void *)nmpc_data->qp_in->r[i]);
          free((void *)nmpc_data->qp_in->lb[i]);
          free((void *)nmpc_data->qp_in->ub[i]);

          free(nmpc_data->qp_out->x[i]);
          free(nmpc_data->qp_out->u[i]);
      }

      free((void *)nmpc_data->qp_in->lb[NN]);
      free((void *)nmpc_data->qp_in->ub[NN]);
      free((void *)nmpc_data->qp_in->q[NN]);
      free((void *)nmpc_data->qp_in->Q[NN]);

      free((void *)nmpc_data->qp_out->x[NN]);

      // general constraints
      free((void *)nmpc_data->qp_in->Cx[0]);
      free((void *)nmpc_data->qp_in->Cu[0]);
      free((void *)nmpc_data->qp_in->lc[0]);
      free((void *)nmpc_data->qp_in->uc[0]);

      free((void *)nmpc_data->qp_in->Cx[1]);
      free((void *)nmpc_data->qp_in->Cu[1]);
      free((void *)nmpc_data->qp_in->lc[1]);
      free((void *)nmpc_data->qp_in->uc[1]);

      free((void *)nmpc_data->qp_in->Cx[NN]);
      free((void *)nmpc_data->qp_in->lc[NN]);
      free((void *)nmpc_data->qp_in->uc[NN]);


      free(nmpc_data->qp_out->pi[0]);
      free(nmpc_data->qp_out->lam[0]);
      for (int_t ii = 1; ii < NN; ii++) {
          free(nmpc_data->qp_out->pi[ii]);
          free(nmpc_data->qp_out->lam[ii]);
      }

      free(nmpc_data->qp_out->lam[NN]);

      // free double pointers
      free(nmpc_data->qp_in->A);
      free(nmpc_data->qp_in->B);
      free(nmpc_data->qp_in->b);
      free(nmpc_data->qp_in->Q);
      free(nmpc_data->qp_in->S);
      free(nmpc_data->qp_in->R);
      free(nmpc_data->qp_in->q);
      free(nmpc_data->qp_in->r);
      // px0 = malloc(sizeof(real_t*));

      free(nmpc_data->qp_in->ub);
      free(nmpc_data->qp_in->lb);
      free(nmpc_data->qp_in->Cx);
      free(nmpc_data->qp_in->Cu);
      free(nmpc_data->qp_in->lc);
      free(nmpc_data->qp_in->uc);

      free(nmpc_data->qp_in->idxb);


      free(nmpc_data->qp_out->x);
      free(nmpc_data->qp_out->u);
      free(nmpc_data->qp_out->pi);
      free(nmpc_data->qp_out->lam);

      // free solver arguments
      free(nmpc_data->hpmpc_args);

      // free work space
      free(nmpc_data->workspace);


      // free memory for nmpc_data
      free(nmpc_data->w);

      free(nmpc_data->lb0);
      free(nmpc_data->lb);
      free(nmpc_data->lbN);

      free(nmpc_data->ub0);
      free(nmpc_data->ub);
      free(nmpc_data->ubN);

      free(nmpc_data->idxb0);
      free(nmpc_data->idxb);
      free(nmpc_data->idxbN);

      // free OCP QP struct
      free(nmpc_data->qp_in);
      free(nmpc_data->qp_out);

      // allocate memory for residual evalutation
      // get internal memory dimension

      free(nmpc_data->residual_x_eval_mem);
      free(nmpc_data->residual_x_out);
      free(nmpc_data->residual_x_in);
      free(nmpc_data->drdx_tran);

      free(nmpc_data->residual_u_eval_mem);
      free(nmpc_data->residual_u_out);
      free(nmpc_data->residual_u_in );
      free(nmpc_data->drdu_tran);

      return 0;

  }

  int_t run_acados(nmpc_data *nmpc_data, rk4_int* rk4_int, acados_options *acados_options){

    int status = -1;
    int sanity_checks = 1;
    real_t qp_step_size = 0;

    // get nmpc data
    const int NN = nmpc_data->NN;   // number of stages
    const int NX = nmpc_data->NX;   // number of states
    const int NU = nmpc_data->NU;   // number of inputs
    const int NB0 = nmpc_data->NB0; // number of input bounds on stage 0
    const int NB = nmpc_data->NB;   // number of input bounds on stage 1 to N-1
    const int NBN = nmpc_data->NBN; // number of input bounds on stage N

    const int NP = nmpc_data->NP;   // numbers of parameters
    const int NR = nmpc_data->NR;   // number of residuals (least-square)

    real_t *lb0 = nmpc_data->lb0;
    real_t *lb = nmpc_data->lb;
    real_t *lbN = nmpc_data->lbN;

    real_t *ub0 = nmpc_data->ub0;
    real_t *ub = nmpc_data->ub;
    real_t *ubN = nmpc_data->ubN;

    real_t *w = nmpc_data->w;
    // int_t max_sqp_iters=  nmpc_data->max_sqp_iters;
    // int_t max_iters= nmpc_data->max_iters;

    const int nls = acados_options->nls;

    // get qp data
    ocp_qp_in *qp_in= nmpc_data->qp_in;
    ocp_qp_out *qp_out= nmpc_data->qp_out;

    ocp_qp_hpmpc_args *hpmpc_args= nmpc_data->hpmpc_args;
    void *workspace = nmpc_data->workspace;

    int_t *nx = (int_t *)qp_in->nx;
    // int_t *nu = (int_t *)qp_in->nu;
    // int_t *nb = (int_t *)qp_in->nb;
    // int_t *nc = (int_t *)qp_in->nc;
    real_t **pA = (real_t **)qp_in->A;
    real_t **pB = (real_t **)qp_in->B;
    real_t **pb = (real_t **)qp_in->b;
    real_t **pQ = (real_t **)qp_in->Q;

    // real_t **pS = qp_in->S;
    real_t **pR = (real_t **)qp_in->R;
    real_t **pq = (real_t **)qp_in->q;
    real_t **pr = (real_t **)qp_in->r;
    int    **hidxb = (int_t **)qp_in->idxb;
    real_t **pub = (real_t **)qp_in->ub;
    real_t **plb = (real_t **)qp_in->lb;

    // get pointers to memory for residual evaluation
    real_t *residual_x_eval_mem = nmpc_data->residual_x_eval_mem;
    real_t *residual_x_out = nmpc_data->residual_x_out;
    real_t *residual_x_in = nmpc_data->residual_x_in;
    real_t *drdx_tran = nmpc_data->drdx_tran;

    real_t *residual_u_eval_mem = nmpc_data->residual_u_eval_mem;
    real_t *residual_u_out = nmpc_data->residual_u_out;
    real_t *residual_u_in = nmpc_data->residual_u_in;
    real_t *drdu_tran = nmpc_data->drdu_tran;

    real_t regQ = nmpc_data->regQ;
    real_t regR = nmpc_data->regR;

    acado_timer timer;

    real_t timings = 0;
    printf("running acados\n");
    // printf("\n------ ITERATION %d ------\n", iter);
    acado_tic(&timer);
    // for ( int_t ii = 0; ii < NX; ii++ ) w[ii] = lb0[ii];
    for (int_t sqp_iter = 0; sqp_iter < acados_options->sqp_steps; sqp_iter++) {

    // Pass state and control to integrator
    for (int_t j = 0; j < NX; j++) rk4_int->x_in[j] = w[0*(NX+NU)+j];
    for (int_t j = 0; j < NU; j++) rk4_int->u_in[j] = w[0*(NX+NU)+NX+j];

    integrate_aircraft_ode(rk4_int);

    // Construct QP matrices
    real_t *r = residual_x_out;
    real_t *drdx = &residual_x_out[NR];
    real_t *drdu = &residual_u_out[NU];
    real_t *rref = &residual_x_out[NR+NR*NX];

    if (nls){
      // build Gauss-Newton Hessian and gradient
      // evaluate residuals

      // states
      for (int_t j = 0; j < NX; j++) residual_x_in[j] = w[0*(NX+NU)+j];
      for (int_t j = 0; j < NP; j++) residual_x_in[j+NX] = rk4_int->p_in[j];

      residual_x_in[NX+NP] = rk4_int->t0_in;
      residual_x_in[NX+NP+1] = rk4_int->h_in;

      residual_x_eval_wrapper(nmpc_data, residual_x_in, residual_x_out, residual_x_eval_mem);

      for (int_t j = 0; j < NX; j++)
        for (int_t k = 0; k < NR; k++) {
          drdx_tran[k*NX + j] = drdx[j*NR + k];
        }

      // init pQ to zeros
      for (int_t j = 0; j < NX*NX; j++) pQ[0][j] = 0.0;
      dgemm_nn_3l(NX, NX, NR, drdx_tran, NX, drdx, NR, pQ[0], NX);

      // check symmetry
      for (int_t j = 0; j < NX; j++) {
        for (int_t k = 0; k < NX; k++) {
          if (fabs(pQ[0][j*NX + k] - pQ[0][k*NX + j]) > acados_options->non_symmetry_tol) {
            printf("nonsymmetric Hessian!");
            sanity_checks = 0;
          }
        }
      }

      // compute eta_tilde
      if ( acados_options->print_level > 2) {
        for (int_t j = 0; j < NR; j++) printf("rref[%i]=%f\n", j, rref[j] );
        for (int_t j = 0; j < NR; j++) rref[j] = -rref[j] + r[j];
        printf("\n\n");
        for (int_t j = 0; j < NR; j++) printf("r[%i]=%f\n", j, r[j] );
      }

      for (int_t j = 0; j < NX; j++) pq[0][j] = 0.0;

      dgemv_n_3l(NX, NR, drdx_tran, NX, rref, pq[0]);

      // inputs
      for (int_t j = 0; j < NU; j++) residual_u_in[j] = w[0*(NX+NU)+NX+j];
      for (int_t j = 0; j < NP; j++) residual_u_in[j+NU] = rk4_int->p_in[j];

      residual_u_in[NU+NP] = rk4_int->t0_in;
      residual_u_in[NU+NP+1] = rk4_int->h_in;

      // !! code for u_residual is different! we have NU residuals instead of NR
      residual_u_eval_wrapper(nmpc_data, residual_u_in, residual_u_out, residual_u_eval_mem);

      r = residual_u_out;
      drdu = &residual_u_out[NU];
      rref = &residual_u_out[NU+NU*NU];

      for (int_t j = 0; j < NU; j++)
        for (int_t k = 0; k < NU; k++) {
          drdu_tran[k*NU + j] = drdu[j*NU + k];
        }

      // init pR to zeros
      for (int_t j = 0; j < NU*NU; j++) pR[0][j] = 0.0;
      dgemm_nn_3l(NU, NU, NU, drdu_tran, NU, drdu, NU, pR[0], NU);

      // compute eta_tilde
      for (int_t j = 0; j < NU; j++) rref[j] = -rref[j] + r[j];
      for (int_t j = 0; j < NU; j++) pr[0][j] = 0.0;
      dgemv_n_3l(NU, NU, drdu_tran, NU, rref, pr[0]);

      // regularize
      for (int_t j = 0; j < NX; j++) pQ[0][j*(NX+1)] += regQ;
      for (int_t j = 0; j < NU; j++) pR[0][j*(NU+1)] += regR;

      if ( acados_options->print_level > 2 ){
        // compute eigenvalues of Q and R (this is for debugging greg don't hate me)
        int_t matz = 0;
        real_t eigenValQ[11];
        real_t eigenVecQ[11*11];
        real_t Q_copy[11*11] = {0.0};
        for (int_t j = 0; j < 11*11; j++) Q_copy[j] = pQ[0][j];
        rs (NX, Q_copy, eigenValQ, matz, eigenVecQ);

        printf("\neigenvalus of Gauss-Newton Q:\n");
        for(int_t j = 0; j < NX; j++) printf("eigenValQ[%i] = %f\n",j,eigenValQ[j]);

        // compute eigenvalues of Q and R
        // int_t matz = 0;
        real_t eigenValR[4];
        real_t eigenVecR[4*4];
        real_t R_copy[4*4] = {0.0};
        for (int_t j = 0; j < 4*4; j++) R_copy[j] = pR[0][j];
        rs (NU, R_copy, eigenValR, matz, eigenVecR);

        printf("\neigenvalus of Gauss-Newton R:\n");
        for(int_t j = 0; j < NU; j++) printf("eigenValR[%i] = %f\n", j, eigenValR[j]);
      }
    } else {
      // constant matrices
      const real_t xref[11] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0};
      const real_t Q[11] = {10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01,};
      for (int_t j = 0; j < NX; j++) {
        for (int_t k = 0; k < NX; k++) {
          pQ[0][j*NX + k] = 0.0;
          if (j==k) pQ[0][j*NX + k] = Q[j];
        }
      }

      const real_t uref[4] = {0.0, 0.0, 0.0, 0.0};
      const real_t R[4] = {0.1, 0.1, 0.1, 0.1};
      for (int_t j = 0; j < NU; j++) {
        for (int_t k = 0; k < NU; k++) {
          pR[0][j*NU + k] = 0.0;
          if (j==k) pR[0][j*NU + k] = R[j];
        }
      }

      for (int_t j = 0; j < NX; j++) {
          pq[0][j] = pQ[0][j*(NX+1)]*(w[0*(NX+NU)+j]-xref[j]);
      }

      for (int_t j = 0; j < NU; j++) {
          pr[0][j] = pR[0][j*(NU+1)]*(w[0*(NX+NU)+NX+j]-uref[j]);
      }
    }

    for (int_t j = 0; j < NX; j++) {
        pb[0][j] = rk4_int->x_out[j] - w[(0+1)*(NX+NU)+j];
        for (int_t k = 0; k < nx[0]; k++) pA[0][j*NX+k] = rk4_int->Sx_out[j*(NX)+k];
    }

    for (int_t j = 0; j < NU; j++)
        for (int_t k = 0; k < NX; k++) pB[0][j*NX+k] = rk4_int->Su_out[NX*j+k];

    for (int_t i = 1; i < NN; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < NX; j++) rk4_int->x_in[j] = w[i*(NX+NU)+j];
        for (int_t j = 0; j < NU; j++) rk4_int->u_in[j] = w[i*(NX+NU)+NX+j];

        integrate_aircraft_ode(rk4_int);

        // Construct QP matrices

        if (nls){
          // build Gauss-Newton Hessian and gradient
          // evaluate residuals

          for (int_t j = 0; j < NX; j++) residual_x_in[j] = w[i*(NX+NU)+j];
          for (int_t j = 0; j < NP; j++) residual_x_in[j+NX] = rk4_int->p_in[j];

          residual_x_in[NX+NP] = rk4_int->t0_in;
          residual_x_in[NX+NP+1] = rk4_int->h_in;

          residual_x_eval_wrapper(nmpc_data, residual_x_in, residual_x_out, residual_x_eval_mem);

          r = residual_x_out;
          drdx = &residual_x_out[NR];
          rref = &residual_x_out[NR+NR*NX];

          for (int_t j = 0; j < NX; j++) {
            for (int_t k = 0; k < NR; k++) {
              drdx_tran[k*NX + j] = drdx[j*NR + k];
            }
          }
          // init pQ to zeros
          for (int_t j = 0; j < NX*NX; j++) pQ[i][j] = 0.0;
          dgemm_nn_3l(NX, NX, NR, drdx_tran, NX, drdx, NR, pQ[i], NX);


          // check symmetry
          for (int_t j = 0; j < NX; j++) {
            for (int_t k = 0; k < NX; k++) {
              if (fabs(pQ[i][j*NX + k] - pQ[i][k*NX + j]) > acados_options->non_symmetry_tol) {
                printf("nonsymmetric Hessian!");
                sanity_checks = 0;
              }
            }
          }

          // compute eta_tilde
          for (int_t j = 0; j < NR; j++) rref[j]= -rref[j] + r[j];
          for (int_t j = 0; j < NX; j++) pq[i][j] = 0.0;

          dgemv_n_3l(NX, NR, drdx_tran, NX, rref, pq[i]);

          // inputs
          for (int_t j = 0; j < NU; j++) residual_u_in[j] = w[i*(NX+NU)+NX+j];
          for (int_t j = 0; j < NP; j++) residual_u_in[j+NU] = rk4_int->p_in[j];

          residual_u_in[NU+NP] = rk4_int->t0_in;
          residual_u_in[NU+NP+1] = rk4_int->h_in;

          // code for u residuals is different! NU residuals instead of NR
          residual_u_eval_wrapper(nmpc_data, residual_u_in, residual_u_out, residual_u_eval_mem);

          r = residual_u_out;
          drdu = &residual_u_out[NU];
          rref = &residual_u_out[NU+NU*NU];

          // real_t drdu_tran[NU*NU] = {0.0};
          for (int_t j = 0; j < NU; j++) {
            for (int_t k = 0; k < NU; k++) {
              drdu_tran[k*NU + j] = drdu[j*NU + k];
            }
          }

          // init pR to zeros
          for (int_t j = 0; j < NU*NU; j++) pR[i][j] = 0.0;

          dgemm_nn_3l(NU, NU, NU, drdu_tran, NU, drdu, NU, pR[i], NU);

          // compute eta_tilde
          for (int_t j = 0; j < NU; j++) rref[j] = -rref[j] + r[j];
          for (int_t j = 0; j < NU; j++) pr[i][j] = 0.0;
          dgemv_n_3l(NU, NU, drdu_tran, NU, rref, pr[i]);

          // regularize
          for (int_t j = 0; j < NX; j++) pQ[i][j*(NX+1)] += regQ;
          for (int_t j = 0; j < NU; j++) pR[i][j*(NU +1)] += regR;

          if ( acados_options->print_level > 2 ){
            // compute eigenvalues of Q and R
            int_t matz = 0;
            real_t eigenValQ[11];
            real_t eigenVecQ[11*11];
            real_t Q_copy[11*11] = {0.0};
            for (int_t j = 0; j < 11*11; j++) Q_copy[j] = pQ[i][j];
            rs (NX, Q_copy, eigenValQ, matz, eigenVecQ);

            printf("\neigenvalus of Gauss-Newton Q:\n");
            for(int_t j = 0; j < NX; j++) printf("eigenValQ[%i] = %f\n",j,eigenValQ[j]);

            // compute eigenvalues of Q and R
            // int_t matz = 0;
            real_t eigenValR[4];
            real_t eigenVecR[4*4];
            real_t R_copy[4*4] = {0.0};
            for (int_t j = 0; j < 4*4; j++) R_copy[j] = pR[i][j];
            rs (NU, R_copy, eigenValR, matz, eigenVecR);

            printf("\neigenvalus of Gauss-Newton R:\n");
            for(int_t j = 0; j < NU; j++) printf("eigenValR[%i] = %f\n", j, eigenValR[j]);
          }
        } else {
          // constant matrices
          const real_t xref[11] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0};
          const real_t Q[11] = {10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01,};
          for (int_t j = 0; j < NX; j++) {
            for (int_t k = 0; k < NX; k++) {
              pQ[i][j*NX + k] = 0.0;
              if (j==k) pQ[i][j*NX + k] = Q[j];
            }
          }

          const real_t uref[4] = {0.0, 0.0, 0.0, 0.0};
          const real_t R[4] = {0.1, 0.1, 0.1, 0.1};
          for (int_t j = 0; j < NU; j++) {
            for (int_t k = 0; k < NU; k++) {
              pR[i][j*NU + k] = 0.0;
              if (j==k) pR[i][j*NU + k] = R[j];
            }
          }

          for (int_t j = 0; j < NX; j++) {
              pq[i][j] = pQ[i][j*(NX+1)]*(w[i*(NX+NU)+j]-xref[j]);
          }
          for (int_t j = 0; j < NU; j++) {
              pr[i][j] = pR[i][j*(NU+1)]*(w[i*(NX+NU)+NX+j]-uref[j]);
          }
        }

        for (int_t j = 0; j < NX; j++) {
            pb[i][j] = rk4_int->x_out[j] - w[(i+1)*(NX+NU)+j];
            for (int_t k = 0; k < nx[i]; k++) pA[i][j*NX+k] = rk4_int->Sx_out[j*(NX)+k];
        }
        for (int_t j = 0; j < NU; j++)
            for (int_t k = 0; k < NX; k++) pB[i][j*NX+k] = rk4_int->Su_out[NX*j+k];

        // assign bounds taking into account "u-x" order (hpmpc) vs "x-u" order
        real_t relax_bounds = 0.000001;
        for (int_t j = 0; j < NB; j++) {
          if (hidxb[i][j] < NU) {
            plb[i][j] = lb[j] - w[i*(NX+NU)+hidxb[i][j] + NX] - relax_bounds;
            pub[i][j] = ub[j] - w[i*(NX+NU)+hidxb[i][j] + NX] + relax_bounds;
          } else {
            plb[i][j] = lb[j] - w[i*(NX+NU)+hidxb[i][j] - NU] - relax_bounds;
            pub[i][j] = ub[j] - w[i*(NX+NU)+hidxb[i][j] - NU] + relax_bounds;
          }
        }
    }

    // assign bounds taking into account "u-x" order (hpmpc) vs "x-u" order
    real_t relax_bounds = 0.000001;
    for (int_t j = 0; j < NBN; j++) {
      if (hidxb[NN][j] < NU) {
        plb[NN][j] = lbN[j] - w[NN*(NX+NU)+hidxb[NN][j] + NX] - relax_bounds;
        pub[NN][j] = ubN[j] - w[NN*(NX+NU)+hidxb[NN][j] + NX] + relax_bounds;
      } else {
        plb[NN][j] = lbN[j] - w[NN*(NX+NU)+hidxb[NN][j] - NU] - relax_bounds;
        pub[NN][j] = ubN[j] - w[NN*(NX+NU)+hidxb[NN][j] - NU] + relax_bounds;
      }
    }

    // assign bounds taking into account "u-x" order (hpmpc) vs "x-u" order
    for (int_t j = 0; j < NB0; j++) {
      if (hidxb[0][j] < NU) {
        plb[0][j] = lb0[j] - w[0*(NX+NU)+hidxb[0][j] + NX] - relax_bounds;
        pub[0][j] = ub0[j] - w[0*(NX+NU)+hidxb[0][j] + NX] + relax_bounds;
      } else {
        plb[0][j] = lb0[j] - w[0*(NX+NU)+hidxb[0][j] - NU] - relax_bounds;
        pub[0][j] = ub0[j] - w[0*(NX+NU)+hidxb[0][j] - NU] + relax_bounds;
      }
    }


    if (acados_options->print_level > 1) {
      for (int_t j = 0; j < NB0; j++) {
        printf("lb0[%i]=%f  ",j,plb[0][j] );
      }

      printf("\n");

      for (int_t j = 0; j < NB0; j++) {
        printf("ub0[%i]=%f  ",j,pub[0][j] );
      }

      printf("\n\n");

      for (int_t j = 0; j < NB; j++) {
        printf("lb[%i]=%f  ",j,plb[1][j] );
      }
      printf("\n\n");

      for (int_t j = 0; j < NB; j++) {
        printf("ub[%i]=%f  ",j,pub[1][j] );
      }

      printf("\n\n");

      for (int_t j = 0; j < NBN; j++) {
        printf("lbN[%i]=%f  ",j,plb[NN][j] );
      }

      printf("\n");

      for (int_t j = 0; j < NBN; j++) {
        printf("ubN[%i]=%f  ",j,pub[NN][j] );
      }
    }

    // for (int_t j = 0; j < NB0; j++) if (plb[0] > pub[0])  {
    //   printf("infeasbile QP!");
    //   sanity_checks = 0;
    // }
    //
    //
    // for (int_t j = 0; j < NB; j++) if (plb[1] > pub[1]) {
    //   printf("infeasbile QP!");
    //   sanity_checks = 0;
    //
    // }
    //
    //
    // for (int_t j = 0; j < NB; j++) if (plb[NN] > pub[NN]) {
    //   printf("infeasbile QP!");
    //   sanity_checks = 0;
    // }

    if (nls){
      // build Gauss-Newton Hessian and gradient
      // evaluate residuals

      for (int_t j = 0; j < NX; j++) residual_x_in[j] = w[NN*(NX+NU)+j];
      for (int_t j = 0; j < NP; j++) residual_x_in[j+NX] = rk4_int->p_in[j];

      residual_x_in[NX+NP] = rk4_int->t0_in;
      residual_x_in[NX+NP+1] = rk4_int->h_in;

      residual_x_eval_wrapper(nmpc_data, residual_x_in, residual_x_out, residual_x_eval_mem);

      r = residual_x_out;
      drdx = &residual_x_out[NR];
      rref = &residual_x_out[NR+NR*NX];

      for (int_t j = 0; j < NX; j++) {
        for (int_t k = 0; k < NR; k++)
          drdx_tran[k*NX + j] = drdx[j*NR + k];
      }

      // init pQ to zeros
      for (int_t j = 0; j < NX*NX; j++) pQ[NN][j] = 0.0;

      dgemm_nn_3l(NX, NX, NR, drdx_tran, NX, drdx, NR, pQ[NN], NX);

      // check symmetry
      for (int_t j = 0; j < NX; j++) {
        for (int_t k = 0; k < NX; k++) {
          if (fabs(pQ[NN][j*NX + k] - pQ[NN][k*NX + j]) > acados_options->non_symmetry_tol) {
            printf("nonsymmetric Hessian!");
            sanity_checks = 0;
          }
        }
      }

      // compute eta_tilde
      for (int_t j = 0; j < NR; j++) rref[j]= -rref[j] + r[j];
      for (int_t j = 0; j < NX; j++) pq[NN][j] = 0.0;

      dgemv_n_3l(NX, NR, drdx_tran, NX, rref, pq[NN]);

      // regularize
      for (int_t j = 0; j < NX; j++) pQ[NN][j*(NX+1)] += regQ;

      // compute eigenvalues of Q
      if ( acados_options->print_level > 2 ){
        int_t matz = 0;
        real_t eigenValQ[11];
        real_t eigenVecQ[11*11];
        real_t Q_copy[11*11] = {0.0};
        for (int_t j = 0; j < 11*11; j++) Q_copy[j] = pQ[NN][j];
        rs (NX, Q_copy, eigenValQ, matz, eigenVecQ);

        printf("\neigenvalus of Gauss-Newton Q:\n");
        for(int_t j = 0; j < NX; j++) printf("eigenValQ[%i] = %f\n", j, eigenValQ[j]);
    }

    } else {
      // constant matrices
      const real_t xref[11] = {1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ,0.0, 0.0, 0.0, 0.0};
      const real_t Q[11] = {10.0, 10.0, 10.0, 10.0, 1.0, 1.0, 1.0, 0.01, 0.01, 0.01, 0.01};
      for (int_t j = 0; j < NX; j++) {
        for (int_t k = 0; k < NX; k++) {
          pQ[NN][j*NX + k] = 0.0;
          if (j==k) pQ[NN][j*NX + k] = Q[j];
        }
      }

      for (int_t j = 0; j < NX; j++) {
          pq[NN][j] = pQ[NN][j*(NX+1)]*(w[NN*(NX+NU)+j]-xref[j]);
      }
    }

    if (sanity_checks) {
      status = ocp_qp_hpmpc_libstr(qp_in, qp_out, hpmpc_args, workspace);
    } else {
      status = 3; // sanity checks failed
    } // int status = 0;

    if (status == 1) printf("status = ACADOS_MAXITER\n");

    if (status == 2) printf("status = ACADOS_MINSTEP\n");

    if (status == 3) printf("sanity checks failed\n");

    // int this formulation we keep x0 as well
    real_t lambda = acados_options->newton_step_size;
    for (int_t j = 0; j < NX; j++) {
      w[0*(NX+NU)+j] += qp_out->x[0][j];
      qp_step_size+=lambda*lambda*qp_out->x[0][j]*qp_out->x[0][j];
    }
    for (int_t j = 0; j < NU; j++) {
      w[0*(NX+NU)+NX+j] += lambda*qp_out->u[0][j];
      qp_step_size+=lambda*lambda*qp_out->u[0][j]*qp_out->u[0][j];
    }

    for (int_t i = 1; i < NN; i++) {
        for (int_t j = 0; j < NX; j++) {
          w[i*(NX+NU)+j] += lambda*qp_out->x[i][j];
          qp_step_size+=lambda*lambda*qp_out->x[i][j]*qp_out->x[i][j];
        }
        for (int_t j = 0; j < NU; j++) {
          w[i*(NX+NU)+NX+j] += lambda*qp_out->u[i][j];
          qp_step_size+=lambda*lambda*qp_out->u[i][j]*qp_out->u[i][j];
        }
    }

    for (int_t j = 0; j < NX; j++) {
      w[NN*(NX+NU)+j] += lambda*qp_out->x[NN][j];
      qp_step_size+=lambda*lambda*qp_out->x[NN][j]*qp_out->x[NN][j];
    }


    // noralize vector
    qp_step_size = sqrt(qp_step_size);
    qp_step_size/=((NX+NU)*NN+NX);
    // for (int_t i = 0; i < NX; i++) x0[i] = w[NX+NU+i];
    // shift_states(w, x_end, NN);
    // shift_controls(w, u_end, NN);
    timings += acado_toc(&timer);
  }


  real_t Th = (real_t)NN*rk4_int->h_in;
  if (acados_options->print_level > 1)  print_states_controls(&w[0], Th, NN, NX, NU);
  if (acados_options->plot_open_loop)   plot_states_controls(&w[0],  Th, NN, NX, NU, nmpc_data->gnuplotPipe);
  if (acados_options->print_level > 0)  printf("\n%.3f ms per iteration. QP step-size = %f\n\n", 1e3*timings, qp_step_size);

  return status;
}

// uncomment for testing
#if 0
int main() {

  // feenableexcept(FE_DIVBYZERO | FE_INVALID | FE_OVERFLOW);

  const int NN  = nn;   // number of stages
  const int NX  = 11;   // number of states
  const int NU  = 4;    // number of inputs
  const int NB0 = 11;    // number of bounds on stage 0
  const int NB  = 4;    // number of bounds on stage 1 to N-1
  const int NBN = 4;    // number of bounds on stage N
  const int NP  = 43;   // numbers of parameters
  const int NR  = 14;   // number of residuals (least-square)

  const real_t REGQ = 0.0;
  const real_t REGR = 0.0;

  init_nmpc_data *init_data;
  init_data = malloc(sizeof *init_data);

  nmpc_data *nmpc_data;
  nmpc_data = malloc(sizeof *nmpc_data);

  rk4_int *rk4_int;
  rk4_int = malloc(sizeof *rk4_int);

  // get dimensions
  init_data->NN  = NN;
  init_data->NX  = NX;
  init_data->NU  = NU;
  init_data->NB0 = NB0;
  init_data->NB  = NB;
  init_data->NBN = NBN;
  init_data->NP  = NP;
  init_data->NR  = NR;

  // define bound indexes
  int_t *idxb0 = malloc(sizeof(int_t)*NB0);
  int_t *idxb = malloc(sizeof(int_t)*NB);
  int_t *idxbN = malloc(sizeof(int_t)*NBN);

  // bound-indexes for stage 0
  idxb0[0] = 4;
  idxb0[1] = 5;
  idxb0[2] = 6;
  idxb0[3] = 7;
  idxb0[4] = 8;
  idxb0[5] = 9;
  idxb0[6] = 10;

  idxb0[7] = 11;
  idxb0[8] = 12;
  idxb0[9] = 13;
  idxb0[10] = 14;

  // // bound-indexes for stage 1 to N-1
  idxb[0] = 11;
  idxb[1] = 12;
  idxb[2] = 13;
  idxb[3] = 14;

  // bound-indexes for stage N
  idxbN[0] = 11;
  idxbN[1] = 12;
  idxbN[2] = 13;
  idxbN[3] = 14;

  init_data->idxb0 = idxb0;
  init_data->idxb = idxb;
  init_data->idxbN = idxbN;

  nmpc_data->res_x = &ocp_xtracking_banana1;
  nmpc_data->res_x_work = &ocp_xtracking_banana1_work;

  nmpc_data->res_u = &ocp_utracking_banana1;
  nmpc_data->res_u_work = &ocp_utracking_banana1_work;

  nmpc_data->nls = 1;

  rk4_int->eval_dynamics = &ocp_integrate_ode_banana1;
  rk4_int->eval_dynamics_work = &ocp_integrate_ode_banana1_work;

  init_acados(nmpc_data, rk4_int, init_data);

  real_t init_angle = 1;

  rk4_int->n_steps = RKSTEPS;
  rk4_int->h_in = STEP_SIZE;

  // initial condition
  nmpc_data->lb0[0] = cos(init_angle/2);
  nmpc_data->lb0[1] = sin(init_angle/2);
  nmpc_data->lb0[2] = 0.0;
  nmpc_data->lb0[3] = 0.0;
  nmpc_data->lb0[4] = 0.0;
  nmpc_data->lb0[5] = 0.0;
  nmpc_data->lb0[6] = 0.0;

  nmpc_data->lb0[7] = -1.0;
  nmpc_data->lb0[8] = -1.0;
  nmpc_data->lb0[9] = -1.0;
  nmpc_data->lb0[10] = -1.0;

  nmpc_data->ub0[0] = cos(init_angle/2);
  nmpc_data->ub0[1] = sin(init_angle/2);
  nmpc_data->ub0[2] = 0.0;
  nmpc_data->ub0[3] = 0.0;
  nmpc_data->ub0[4] = 0.0;
  nmpc_data->ub0[5] = 0.0;
  nmpc_data->ub0[6] = 0.0;

  nmpc_data->ub0[7] = 1.0;
  nmpc_data->ub0[8] = 1.0;
  nmpc_data->ub0[9] = 1.0;
  nmpc_data->ub0[10] = 1.0;

  // bounds on stage 1 to N-1
  nmpc_data->lb[0] = -1.0;
  nmpc_data->lb[1] = -1.0;
  nmpc_data->lb[2] = -1.0;
  nmpc_data->lb[3] = -1.0;

  nmpc_data->ub[0] = 1.0;
  nmpc_data->ub[1] = 1.0;
  nmpc_data->ub[2] = 1.0;
  nmpc_data->ub[3] = 1.0;

  // bounds on stage N
  nmpc_data->lbN[0] = -1.0;
  nmpc_data->lbN[1] = -1.0;
  nmpc_data->lbN[2] = -1.0;
  nmpc_data->lbN[3] = -1.0;

  nmpc_data->ubN[0] = 1.0;
  nmpc_data->ubN[1] = 1.0;
  nmpc_data->ubN[2] = 1.0;
  nmpc_data->ubN[3] = 1.0;

  real_t p[NP_dummy] = {9.8, 0, 0, 0, 9.8, 0, 0, 0, 18.8, -0.75, 0, 0, 0, 217,
    1, 1, 1, 1, 9, 1, -1, 1, -1, 1, 0, 40, 40, 40, 0, 0, 0, -1, 1, 1, 1, 15, 15,
    1, 1000, 1000, 1000, 1000, 1000};

  for (int_t i = 0; i < NP; i++) rk4_int->p_in[i] = p[i];


  for (int i = 0; i <= NN; i++)
    for (int_t j = 0; j < NB0; j++)
      nmpc_data->w[i*(NX+NU)+j] = nmpc_data->ub0[j];

  for (int i = 0; i < NU; i++) nmpc_data->w[NX-NU+i] = 0.0;
  for (int i = 0; i < NB0; i++) nmpc_data->w[i] = nmpc_data->lb0[i];

  // for (int i = 0; i <= NN; i++) nmpc_data->w[i*(NX+NU)+1] =  sqrt(1-0.99*0.99);

  nmpc_data->regQ = REGQ;
  nmpc_data->regR = REGR;

  for (int_t i = 0; i < NSIM; i++){
    // normalize quaternion to avoid blow-up
    for (int_t j = 0; j <= NN; j++){
      real_t norm = 0.0;
      // if (nmpc_data->w[j*(NX+NU)]>0.99) nmpc_data->w[j*(NX+NU)] = 0.97;
      for (int_t k = 0; k < 4; k++) norm+=nmpc_data->w[j*(NX+NU) +k]*nmpc_data->w[j*(NX+NU) +k];
      for (int_t k = 0; k < 4; k++) nmpc_data->w[j*(NX+NU) +k] = nmpc_data->w[j*(NX+NU) +k]/sqrt(norm);
    }

    run_acados(nmpc_data, rk4_int);

    for(int_t j = 0; j < NB0; j++) {
      nmpc_data->lb0[j] = nmpc_data->w[NX+NU +j];
      nmpc_data->ub0[j] = nmpc_data->w[NX+NU +j];
    }
  }

}
#   endif
