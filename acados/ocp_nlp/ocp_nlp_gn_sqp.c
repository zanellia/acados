/*
 *    This file is part of acados.
 *
 *    acados is free software; you can redistribute it and/or
 *    modify it under the terms of the GNU Lesser General Public
 *    License as published by the Free Software Foundation; either
 *    version 3 of the License, or (at your option) any later version.
 *
 *    acados is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 *    Lesser General Public License for more details.
 *
 *    You should have received a copy of the GNU Lesser General Public
 *    License along with acados; if not, write to the Free Software Foundation,
 *    Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */

#include "acados/ocp_nlp/ocp_nlp_gn_sqp.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/sim/sim_common.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados/utils/math.h"

int_t ocp_nlp_gn_sqp_calculate_workspace_size(const ocp_nlp_in *in, void *args_,
    ocp_nlp_gn_sqp_memory* mem) {
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args*) args_;

    int_t size;

    size = sizeof(ocp_nlp_gn_sqp_work);

    // allocate mem for least-squares cost
    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *)in->cost;

    size += ocp_nlp_calculate_workspace_size(in, args->common);

    int_t raw_workspace_size = 0;
    if (!args->lin_res) {

        int_t *nr =  cost->nr;
        int_t nr_ = 0;
        for (int_t i = 0; i < in->N; i++ ) nr_+=nr[i];

        size += sizeof(int_t)*nr_;
        int_t N = in->N;
        const int_t *nx = in->nx;
        const int_t *nu = in->nu;

        int_t max_ls_res_in_size = 0;       // ls_res_in
        int_t max_ls_res_out_size = 0;      // ls_res_out
        int_t max_ls_drdw_tran_size = 0;    // drdw_tran
        int_t max_Hess_gn_size = 0;         // Hess_gn
        int_t max_grad_gn_size = 0;         // grad_gn
        int_t max_rref_size = 0;         // grad_gn

        for (int_t i = 0; i <= N; i++) {
            if (max_ls_res_in_size < nx[i] + nu[i])
                max_ls_res_in_size = nx[i] + nu[i];

            if (max_ls_res_out_size < nr[i] + nr[i]*(nx[i] + nu[i]))
                max_ls_res_out_size = nr[i] + nr[i]*(nx[i] + nu[i]);

            if (max_ls_drdw_tran_size < nr[i]*(nx[i] + nu[i]))
                max_ls_drdw_tran_size = nr[i]*(nx[i] + nu[i]);

            if (max_Hess_gn_size < (nx[i] + nu[i])*(nx[i] + nu[i]))
                max_Hess_gn_size = (nx[i] + nu[i])*(nx[i] + nu[i]);

            if (max_grad_gn_size < (nx[i] + nu[i]))
                max_grad_gn_size = (nx[i] + nu[i]);

            if (max_rref_size < nr[i])
                max_rref_size = nr[i];

        }

        raw_workspace_size +=sizeof(real_t)*max_ls_res_out_size;
        raw_workspace_size +=sizeof(real_t)*max_ls_res_in_size;
        raw_workspace_size +=sizeof(real_t)*max_ls_drdw_tran_size;
        raw_workspace_size +=sizeof(real_t)*max_Hess_gn_size;
        raw_workspace_size +=sizeof(real_t)*max_grad_gn_size;
        raw_workspace_size +=sizeof(real_t)*max_rref_size;

        mem->raw_workspace_size = raw_workspace_size;

        size+=raw_workspace_size;

    }

    return size;
}

static void ocp_nlp_gn_sqp_cast_workspace(ocp_nlp_gn_sqp_work *work,
                                          ocp_nlp_gn_sqp_memory *mem) {
    char *ptr = (char *)work;

    ptr += sizeof(ocp_nlp_gn_sqp_work);
    work->raw = (real_t *)ptr;
    ptr += mem->raw_workspace_size;
    work->common = (ocp_nlp_work *)ptr;
    ocp_nlp_cast_workspace(work->common, mem->common);
}


static void initialize_objective(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_args *args,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost*) nlp_in->cost;

    const int_t *nr = cost->nr;

    real_t **qp_Q = (real_t **) gn_sqp_mem->qp_solver->qp_in->Q;
    real_t **qp_S = (real_t **) gn_sqp_mem->qp_solver->qp_in->S;
    real_t **qp_R = (real_t **) gn_sqp_mem->qp_solver->qp_in->R;

    real_t **qp_q = (real_t **) gn_sqp_mem->qp_solver->qp_in->q;
    real_t **qp_r = (real_t **) gn_sqp_mem->qp_solver->qp_in->r;

    if (args->lin_res) {
        // TODO(rien): only for least squares cost with state and control reference atm
        for (int_t i = 0; i <= N; i++) {
            for (int_t j = 0; j < nx[i]; j++) {
                for (int_t k = 0; k < nx[i]; k++) {
                    qp_Q[i][j * nx[i] + k] = cost->W[i][j * (nx[i] + nu[i]) + k];
                }
                for (int_t k = 0; k < nu[i]; k++) {
                    qp_S[i][j * nu[i] + k] =
                        cost->W[i][j * (nx[i] + nu[i]) + nx[i] + k];
                }
            }
            for (int_t j = 0; j < nu[i]; j++) {
                for (int_t k = 0; k < nu[i]; k++) {
                    qp_R[i][j * nu[i] + k] =
                        cost->W[i][(nx[i] + j) * (nx[i] + nu[i]) + nx[i] + k];
                }
            }
        }
    } else {

        for (int_t i = 0; i <= N; i++) {

            char *ptr = (char *)work->raw;
            real_t *ls_res_out = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nr[i] + (nx[i]+nu[i])*nr[i]);

            real_t *r = ls_res_out;  // TODO(Andrea): allocate this
            real_t *drdw = &ls_res_out[nr[i]];

            real_t *ls_res_in = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i]);

            real_t *drdw_tran = (real_t *)ptr;
            ptr+=sizeof(real_t)*nr[i]*(nx[i]+nu[i]);

            real_t *Hess_gn = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i])*(nx[i]+nu[i]);

            real_t *grad_gn = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i]);

            real_t *rref = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nr[i]);

            // (dense) Gauss-Newton Hessian
            for (int_t j = 0; j < nx[i]; j++)
                ls_res_in[j] = gn_sqp_mem->common->x[i][j];

            for (int_t j = 0; j < nu[i]; j++)
                ls_res_in[nx[i]+j] = gn_sqp_mem->common->u[i][j];

            cost->ls_res_eval[i](ls_res_in, ls_res_out, cost->ls_res[i]);

            for (int_t j = 0; j < nx[i] + nu[i]; j++)
              for (int_t k = 0; k < nr[i]; k++) {
                drdw_tran[k*(nx[i] + nu[i]) + j] = drdw[j*nr[i] + k];
              }

            // init Hess_gn to zeros
            for (int_t j = 0; j < (nx[i]+nu[i])*(nx[i]+nu[i]); j++) Hess_gn[j] = 0.0;
            dgemm_nn_3l(nx[i]+nu[i], nx[i]+nu[i], nr[i], drdw_tran, nx[i]+nu[i],
                drdw, nr[i], Hess_gn, nx[i]+nu[i]);

            // compute gradient
            for (int_t j = 0; j < nr[i]; j++) rref[j] = -cost->y_ref[i][j] + r[j];
            // for (int_t j = 0; j < nr[i]; j++) printf("rref[%i]=%f\n", j, rref[j] );
            // printf("\n\n");
            // for (int_t j = 0; j < nr[i]; j++) printf("r[%i]=%f\n", j, r[j] );

            for (int_t j = 0; j < nx[i] + nu[i]; j++) grad_gn[j] = 0.0;

            dgemv_n_3l(nx[i] + nu[i], nr[i], drdw_tran, nx[i] + nu[i], rref, grad_gn);

            // copy dense Hessian and gradient into qp struct
            for (int_t j = 0; j < nx[i]; j++)
                for (int_t k = 0; k < nx[i]; k++)
                    qp_Q[i][k + j*nx[i]] = Hess_gn[k + j*(nx[i] + nu[i])];

            for (int_t j = 0; j < nu[i]; j++)
                for (int_t k = 0; k < nu[i]; k++)
                    qp_R[i][k + j*nu[i]] =
                        Hess_gn[nx[i]*(nx[i]+nu[i]) + k + nx[i] + j*(nx[i] + nu[i])];
            // for (int_t j = 0; j < nx[i]*nu[i]; j++) qp_S[i][j] =
            // Hess_gn[j + nx[i]*nx[i]];  // TODO(Andrea: untested)

            for (int_t j = 0; j < nx[i]; j++) qp_q[i][j] = grad_gn[j];
            for (int_t j = 0; j < nu[i]; j++) qp_r[i][j] = grad_gn[nx[i] + j];

        }
    }
}

static void initialize_trajectories(
    const ocp_nlp_in *nlp_in,
    ocp_nlp_gn_sqp_memory *gn_sqp_mem,
    ocp_nlp_gn_sqp_work *work) {

    const int_t N = nlp_in->N;
    const int_t *nx = nlp_in->nx;
    const int_t *nu = nlp_in->nu;
    real_t *w = work->common->w;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx + j] = gn_sqp_mem->common->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            w[w_idx + nx[i] + j] = gn_sqp_mem->common->u[i][j];
        }
        w_idx += nx[i] + nu[i];
    }
}


static void multiple_shooting(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_args *args,
    ocp_nlp_gn_sqp_work *work,
    ocp_nlp_gn_sqp_memory *mem, real_t *w) {

    ocp_nlp_ls_cost *cost = (ocp_nlp_ls_cost *)nlp->cost;
    const int_t *nr = cost->nr;
    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    sim_solver *sim = nlp->sim;
    real_t **y_ref = cost->y_ref;

    real_t **qp_A = (real_t **) mem->qp_solver->qp_in->A;
    real_t **qp_B = (real_t **) mem->qp_solver->qp_in->B;
    real_t **qp_b = (real_t **) mem->qp_solver->qp_in->b;
    real_t **qp_q = (real_t **) mem->qp_solver->qp_in->q;
    real_t **qp_r = (real_t **) mem->qp_solver->qp_in->r;
    real_t **qp_R = (real_t **) mem->qp_solver->qp_in->R;
    real_t **qp_Q = (real_t **) mem->qp_solver->qp_in->Q;
    real_t **qp_lb = (real_t **) mem->qp_solver->qp_in->lb;
    real_t **qp_ub = (real_t **) mem->qp_solver->qp_in->ub;

    int_t w_idx = 0;

    for (int_t i = 0; i < N; i++) {
        // Pass state and control to integrator
        for (int_t j = 0; j < nx[i]; j++) sim[i].in->x[j] = w[w_idx+j];
        for (int_t j = 0; j < nu[i]; j++) sim[i].in->u[j] = w[w_idx+nx[i]+j];
        sim[i].fun(sim[i].in, sim[i].out, sim[i].args, sim[i].mem, sim[i].work);

        // TODO(rien): transition functions for changing dimensions not yet implemented!
        for (int_t j = 0; j < nx[i]; j++) {
            qp_b[i][j] = sim[i].out->xn[j] - w[w_idx+nx[i]+nu[i]+j];
            for (int_t k = 0; k < nx[i]; k++)
                qp_A[i][j*nx[i]+k] = sim[i].out->S_forw[j*nx[i]+k];
        }
        for (int_t j = 0; j < nu[i]; j++)
            for (int_t k = 0; k < nx[i]; k++)
                qp_B[i][j*nx[i]+k] = sim[i].out->S_forw[(nx[i]+j)*nx[i]+k];

        // Update bounds:
        for (int_t j = 0; j < nlp->nb[i]; j++) {
#ifdef FLIP_BOUNDS
            if (nlp->idxb[i][j] < nu[i]) {
                qp_lb[i][j] = nlp->lb[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
                qp_ub[i][j] = nlp->ub[i][j] - w[w_idx + nx[i] + nlp->idxb[i][j]];
            } else {
                qp_lb[i][j] = nlp->lb[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
                qp_ub[i][j] = nlp->ub[i][j] - w[w_idx - nu[i] + nlp->idxb[i][j]];
            }
#else
            qp_lb[i][j] = nlp->lb[i][j] - w[w_idx+nlp->idxb[i][j]];
            qp_ub[i][j] = nlp->ub[i][j] - w[w_idx+nlp->idxb[i][j]];
#endif
        }

        if (args->lin_res) {
            // Update gradients
            // TODO(rien): only for diagonal Q, R matrices atm
            // TODO(rien): only for least squares cost with state and control reference atm
            sim_RK_opts *opts = (sim_RK_opts*) sim[i].args;
            for (int_t j = 0; j < nx[i]; j++) {
                qp_q[i][j] = cost->W[i][j*(nx[i]+nu[i]+1)]*(w[w_idx+j]-y_ref[i][j]);
                // adjoint-based gradient correction:
                if (opts->scheme.type != exact) qp_q[i][j] += sim[i].out->grad[j];
            }
            for (int_t j = 0; j < nu[i]; j++) {
                qp_r[i][j] = cost->W[i][(nx[i]+j)*(nx[i]+nu[i]+1)]*(w[w_idx+nx[i]+j]-y_ref[i][nx[i]+j]);
                // adjoint-based gradient correction:
                if (opts->scheme.type != exact) qp_r[i][j] += sim[i].out->grad[nx[i]+j];
            }
        } else {
            char *ptr = (char *)work->raw;
            real_t *ls_res_out = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nr[i] + (nx[i]+nu[i])*nr[i]);

            real_t *r = ls_res_out;  // TODO(Andrea): allocate this
            real_t *drdw = &ls_res_out[nr[i]];

            real_t *ls_res_in = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i]);

            real_t *drdw_tran = (real_t *)ptr;
            ptr+=sizeof(real_t)*nr[i]*(nx[i]+nu[i]);

            real_t *Hess_gn = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i])*(nx[i]+nu[i]);

            real_t *grad_gn = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nx[i]+nu[i]);

            real_t *rref = (real_t *)ptr;
            ptr+=sizeof(real_t)*(nr[i]);

            // (dense) Gauss-Newton Hessian
            for (int_t j = 0; j < nx[i]; j++)
                ls_res_in[j] = w[w_idx+j];

            for (int_t j = 0; j < nu[i]; j++)
                ls_res_in[nx[i]+j] = w[w_idx+nx[i]+j];

            cost->ls_res_eval[i](ls_res_in, ls_res_out, cost->ls_res[i]);

            for (int_t j = 0; j < nx[i] + nu[i]; j++)
              for (int_t k = 0; k < nr[i]; k++) {
                drdw_tran[k*(nx[i] + nu[i]) + j] = drdw[j*nr[i] + k];
              }

            // init Hess_gn to zeros
            for (int_t j = 0; j < (nx[i]+nu[i])*(nx[i]+nu[i]); j++) Hess_gn[j] = 0.0;
            dgemm_nn_3l(nx[i]+nu[i], nx[i]+nu[i], nr[i], drdw_tran, nx[i]+nu[i],
                drdw, nr[i], Hess_gn, nx[i]+nu[i]);

            // compute gradient
            for (int_t j = 0; j < nr[i]; j++) rref[j] = -cost->y_ref[i][j] + r[j];
            // for (int_t j = 0; j < nr[i]; j++) printf("rref[%i]=%f\n", j, rref[j] );
            // printf("\n\n");
            // for (int_t j = 0; j < nr[i]; j++) printf("r[%i]=%f\n", j, r[j] );

            for (int_t j = 0; j < nx[i] + nu[i]; j++) grad_gn[j] = 0.0;

            dgemv_n_3l(nx[i] + nu[i], nr[i], drdw_tran, nx[i] + nu[i], rref, grad_gn);

            // copy dense Hessian and gradient into qp struct

            for (int_t j = 0; j < nx[i]; j++)
                for (int_t k = 0; k < nx[i]; k++)
                    qp_Q[i][k + j*nx[i]] = Hess_gn[k + j*(nx[i] + nu[i])];

            for (int_t j = 0; j < nu[i]; j++)
                for (int_t k = 0; k < nu[i]; k++)
                    qp_R[i][k + j*nu[i]] =
                        Hess_gn[nx[i]*(nx[i]+nu[i]) + k + nx[i] + j*(nx[i] + nu[i])];

            // for (int_t j = 0; j < nx[i]*nu[i]; j++) qp_S[i][j] =
            // Hess_gn[j + nx[i]*nx[i]];  // TODO(Andrea: untested)

            for (int_t j = 0; j < nx[i]; j++) qp_q[i][j] = grad_gn[j];
            for (int_t j = 0; j < nu[i]; j++) qp_r[i][j] = grad_gn[nx[i] + j];

        }

        w_idx += nx[i]+nu[i];
    }

    for (int_t j = 0; j < nlp->nb[N]; j++) {
#ifdef FLIP_BOUNDS
        if (nlp->idxb[N][j] < nu[N]) {
            qp_lb[N][j] = nlp->lb[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
            qp_ub[N][j] = nlp->ub[N][j] - w[w_idx + nx[N] + nlp->idxb[N][j]];
        } else {
            qp_lb[N][j] = nlp->lb[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
            qp_ub[N][j] = nlp->ub[N][j] - w[w_idx - nu[N] + nlp->idxb[N][j]];
        }
#else
        qp_lb[N][j] = nlp->lb[N][j] - w[w_idx+nlp->idxb[N][j]];
        qp_ub[N][j] = nlp->ub[N][j] - w[w_idx+nlp->idxb[N][j]];
#endif
    }

    for (int_t j = 0; j < nx[N]; j++)
        qp_q[N][j] = cost->W[N][j*(nx[N]+nu[N]+1)]*(w[w_idx+j]-y_ref[N][j]);
}


static void update_variables(const ocp_nlp_in *nlp, ocp_nlp_gn_sqp_memory *mem, real_t *w) {
    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;
    sim_solver *sim = nlp->sim;

    for (int_t i = 0; i < N; i++)
        for (int_t j = 0; j < nx[i+1]; j++)
            sim[i].in->S_adj[j] = -mem->qp_solver->qp_out->pi[i][j];

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            w[w_idx+j] += mem->qp_solver->qp_out->x[i][j];
        }
        for (int_t j = 0; j < nu[i]; j++)
            w[w_idx+nx[i]+j] += mem->qp_solver->qp_out->u[i][j];
        w_idx += nx[i]+nu[i];
    }
}


static void store_trajectories(const ocp_nlp_in *nlp, ocp_nlp_memory *memory, ocp_nlp_out *out,
    real_t *w) {

    const int_t N = nlp->N;
    const int_t *nx = nlp->nx;
    const int_t *nu = nlp->nu;

    int_t w_idx = 0;
    for (int_t i = 0; i <= N; i++) {
        for (int_t j = 0; j < nx[i]; j++) {
            memory->x[i][j] = w[w_idx+j];
            out->x[i][j] = w[w_idx+j];
        }
        for (int_t j = 0; j < nu[i]; j++) {
            memory->u[i][j] = w[w_idx+nx[i]+j];
            out->u[i][j] = w[w_idx+nx[i]+j];
        }
        w_idx += nx[i] + nu[i];
    }
}


// Simple fixed-step Gauss-Newton based SQP routine
int_t ocp_nlp_gn_sqp(const ocp_nlp_in *nlp_in, ocp_nlp_out *nlp_out, void *nlp_args_,
    void *nlp_mem_, void *nlp_work_) {

    ocp_nlp_gn_sqp_memory *gn_sqp_mem = (ocp_nlp_gn_sqp_memory *) nlp_mem_;
    ocp_nlp_gn_sqp_work *work = (ocp_nlp_gn_sqp_work*) nlp_work_;
    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args*) nlp_args_;
    ocp_nlp_gn_sqp_cast_workspace(work, gn_sqp_mem);

    initialize_objective(nlp_in, args, gn_sqp_mem, work);
    initialize_trajectories(nlp_in, gn_sqp_mem, work);

    // TODO(roversch): Do we need this here?
    int_t **qp_idxb = (int_t **) gn_sqp_mem->qp_solver->qp_in->idxb;
    for (int_t i = 0; i <= nlp_in->N; i++) {
        for (int_t j = 0; j < nlp_in->nb[i]; j++) {
            qp_idxb[i][j] = nlp_in->idxb[i][j];
        }
    }

    ocp_nlp_gn_sqp_args * nlp_args = (ocp_nlp_gn_sqp_args *) nlp_args_;
    int_t max_sqp_iterations = nlp_args->common->maxIter;

    acados_timer timer;
    real_t total_time = 0;
    acados_tic(&timer);
    for (int_t sqp_iter = 0; sqp_iter < max_sqp_iterations; sqp_iter++) {

        multiple_shooting(nlp_in, nlp_args, work, gn_sqp_mem, work->common->w);

        int_t qp_status = gn_sqp_mem->qp_solver->fun(
            gn_sqp_mem->qp_solver->qp_in,
            gn_sqp_mem->qp_solver->qp_out,
            gn_sqp_mem->qp_solver->args,
            gn_sqp_mem->qp_solver->mem,
            gn_sqp_mem->qp_solver->work);

        if (qp_status != 0) {
            printf("QP solver returned error status %d\n", qp_status);
            return -1;
        }

        update_variables(nlp_in, gn_sqp_mem, work->common->w);

        for (int_t i = 0; i < nlp_in->N; i++) {
            sim_RK_opts *opts = (sim_RK_opts*) nlp_in->sim[i].args;
            nlp_in->sim[i].in->sens_adj = (opts->scheme.type != exact);
            if (nlp_in->freezeSens) {  // freeze inexact sensitivities after first SQP iteration !!
                opts->scheme.freeze = true;
            }
        }
    }

    total_time += acados_toc(&timer);
    store_trajectories(nlp_in, gn_sqp_mem->common, nlp_out, work->common->w);
    return 0;
}

void ocp_nlp_gn_sqp_create_memory(const ocp_nlp_in *in, void *args_, void *memory_) {

    ocp_nlp_gn_sqp_args *args = (ocp_nlp_gn_sqp_args *)args_;
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)memory_;

    ocp_qp_in *dummy_qp = create_ocp_qp_in(in->N, in->nx, in->nu, in->nb, in->nc);
    int_t **idxb = (int_t **) dummy_qp->idxb;
    for (int_t i = 0; i < in->N; i++)
        for (int_t j = 0; j < in->nb[i]; j++)
            idxb[i][j] = in->idxb[i][j];
    mem->qp_solver = create_ocp_qp_solver(dummy_qp, args->qp_solver_name, NULL);

    ocp_nlp_create_memory(in, mem->common);
}

void ocp_nlp_gn_sqp_free_memory(void *mem_) {
    ocp_nlp_gn_sqp_memory *mem = (ocp_nlp_gn_sqp_memory *)mem_;

    int_t N = mem->qp_solver->qp_in->N;
    ocp_nlp_free_memory(N, mem->common);

    mem->qp_solver->destroy(mem->qp_solver->mem, mem->qp_solver->work);

    free(mem->qp_solver->qp_in);
    free(mem->qp_solver->qp_out);
    free(mem->qp_solver->args);
    free(mem->qp_solver);
    // TODO(dimitris): where do we free the integrators?
}
