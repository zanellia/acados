/*
 * Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
 * Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
 * Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
 * Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
 *
 * This file is part of acados.
 *
 * The 2-Clause BSD License
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 * this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.;
 */


#include "acados/ocp_nlp/ocp_nlp_zo_sqp.h"

// external
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#if defined(ACADOS_WITH_OPENMP)
#include <omp.h>
#endif

// blasfeo
#include "blasfeo/include/blasfeo_d_aux.h"
#include "blasfeo/include/blasfeo_d_aux_ext_dep.h"
#include "blasfeo/include/blasfeo_d_blas.h"
// acados
#include "acados/ocp_nlp/ocp_nlp_common.h"
#include "acados/ocp_nlp/ocp_nlp_dynamics_cont.h"
#include "acados/ocp_nlp/ocp_nlp_reg_common.h"
#include "acados/ocp_qp/ocp_qp_common.h"
#include "acados/ocp_qp/ocp_qp_full_condensing.h"
#include "acados/dense_qp/dense_qp_qpoases.h"
#include "acados/utils/mem.h"
#include "acados/utils/print.h"
#include "acados/utils/timing.h"
#include "acados/utils/types.h"
#include "acados_c/ocp_qp_interface.h"

#define HPIPM_ZO 1


/************************************************
 * options
 ************************************************/

acados_size_t ocp_nlp_zo_sqp_opts_calculate_size(void *config_, void *dims_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_zo_sqp_opts);

    size += ocp_nlp_opts_calculate_size(config, dims);

    return size;
}



void *ocp_nlp_zo_sqp_opts_assign(void *config_, void *dims_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;

    char *c_ptr = (char *) raw_memory;

    ocp_nlp_zo_sqp_opts *opts = (ocp_nlp_zo_sqp_opts *) c_ptr;
    c_ptr += sizeof(ocp_nlp_zo_sqp_opts);

    opts->nlp_opts = ocp_nlp_opts_assign(config, dims, c_ptr);
    c_ptr += ocp_nlp_opts_calculate_size(config, dims);

    assert((char *) raw_memory + ocp_nlp_zo_sqp_opts_calculate_size(config, dims) >= c_ptr);

    return opts;
}



void ocp_nlp_zo_sqp_opts_initialize_default(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    // int ii;

    // this first !!!
    ocp_nlp_opts_initialize_default(config, dims, nlp_opts);

    // SQP opts
    opts->max_outer_iter = 100;
    opts->max_inner_iter = 10;
    opts->tol_stat = 1e-6;
    opts->tol_eq   = 1e-6;
    opts->tol_ineq = 1e-6;
    opts->tol_comp = 1e-6;

    opts->ext_qp_res = 0;

    opts->qp_warm_start = 0;
    opts->warm_start_first_qp = false;
    opts->rti_phase = 0;
    opts->print_level = 0;
    opts->initialize_t_slacks = 0;

    // overwrite default submodules opts

    // qp tolerance
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", &opts->tol_stat);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", &opts->tol_eq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", &opts->tol_ineq);
    qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", &opts->tol_comp);

    return;
}



void ocp_nlp_zo_sqp_opts_update(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_update(config, dims, nlp_opts);

    return;
}



void ocp_nlp_zo_sqp_opts_set(void *config_, void *opts_, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = (ocp_nlp_zo_sqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    int ii;

    char module[MAX_STR_LEN];
    char *ptr_module = NULL;
    int module_length = 0;

    // extract module name
    char *char_ = strchr(field, '_');
    if (char_!=NULL)
    {
        module_length = char_-field;
        for (ii=0; ii<module_length; ii++)
            module[ii] = field[ii];
        module[module_length] = '\0'; // add end of string
        ptr_module = module;
    }

    // pass options to QP module
    if ( ptr_module!=NULL && (!strcmp(ptr_module, "qp")) )
    {
        ocp_nlp_opts_set(config, nlp_opts, field, value);

        if (!strcmp(field, "qp_warm_start"))
        {
            int* i_ptr = (int *) value;
            opts->qp_warm_start = *i_ptr;
        }
    }
    else // nlp opts
    {
        if (!strcmp(field, "max_iter"))
        {
            int* max_outer_iter = (int *) value;
            opts->max_outer_iter = *max_outer_iter;
        }
        else if (!strcmp(field, "tol_stat"))
        {
            double* tol_stat = (double *) value;
            opts->tol_stat = *tol_stat;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_stat", value);
        }
        else if (!strcmp(field, "tol_eq"))
        {
            double* tol_eq = (double *) value;
            opts->tol_eq = *tol_eq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_eq", value);
        }
        else if (!strcmp(field, "tol_ineq"))
        {
            double* tol_ineq = (double *) value;
            opts->tol_ineq = *tol_ineq;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_ineq", value);
        }
        else if (!strcmp(field, "tol_comp"))
        {
            double* tol_comp = (double *) value;
            opts->tol_comp = *tol_comp;
            // TODO: set accuracy of the qp_solver to the minimum of current QP accuracy and the one specified.
            config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts, "tol_comp", value);
        }
        else if (!strcmp(field, "ext_qp_res"))
        {
            int* ext_qp_res = (int *) value;
            opts->ext_qp_res = *ext_qp_res;
        }
        else if (!strcmp(field, "warm_start_first_qp"))
        {
            bool* warm_start_first_qp = (bool *) value;
            opts->warm_start_first_qp = *warm_start_first_qp;
        }
        else if (!strcmp(field, "rti_phase"))
        {
            int* rti_phase = (int *) value;
            if (*rti_phase < 0 || *rti_phase > 0) {
                printf("\nerror: ocp_nlp_zo_sqp_opts_set: invalid value for rti_phase field."); 
                printf("possible values are: 0\n");
                exit(1);
            } else opts->rti_phase = *rti_phase;
        }
        else if (!strcmp(field, "print_level"))
        {
            int* print_level = (int *) value;
            if (*print_level < 0)
            {
                printf("\nerror: ocp_nlp_zo_sqp_opts_set: invalid value for print_level field, need int >=0, got %d.", *print_level);
                exit(1);
            }
            opts->print_level = *print_level;
        }
        else if (!strcmp(field, "initialize_t_slacks"))
        {
            int* initialize_t_slacks = (int *) value;
            if (*initialize_t_slacks != 0 && *initialize_t_slacks != 1)
            {
                printf("\nerror: ocp_nlp_zo_sqp_opts_set: invalid value for initialize_t_slacks field, need int 0 or 1, got %d.", *initialize_t_slacks);
                exit(1);
            }
            opts->initialize_t_slacks = *initialize_t_slacks;
        }
        else
        {
            ocp_nlp_opts_set(config, nlp_opts, field, value);
        }
    }

    return;

}



void ocp_nlp_zo_sqp_opts_set_at_stage(void *config_, void *opts_, size_t stage, const char *field, void* value)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = (ocp_nlp_zo_sqp_opts *) opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    ocp_nlp_opts_set_at_stage(config, nlp_opts, stage, field, value);

    return;

}



/************************************************
 * memory
 ************************************************/

acados_size_t ocp_nlp_zo_sqp_memory_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    acados_size_t size = 0;

    size += sizeof(ocp_nlp_zo_sqp_memory);

    // nlp mem
    size += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    int stat_m = opts->max_outer_iter+1;
    int stat_n = 6;
    if (opts->ext_qp_res)
        stat_n += 4;
    size += stat_n*stat_m*sizeof(double);

    size += 3*8;  // align

    make_int_multiple_of(8, &size);

    return size;
}



void *ocp_nlp_zo_sqp_memory_assign(void *config_, void *dims_, void *opts_, void *raw_memory)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    // ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;
    // ocp_nlp_dynamics_config **dynamics = config->dynamics;
    // ocp_nlp_cost_config **cost = config->cost;
    // ocp_nlp_constraints_config **constraints = config->constraints;

    char *c_ptr = (char *) raw_memory;

    // int N = dims->N;
    // int *nx = dims->nx;
    // int *nu = dims->nu;
    // int *nz = dims->nz;

    // initial align
    align_char_to(8, &c_ptr);

    ocp_nlp_zo_sqp_memory *mem = (ocp_nlp_zo_sqp_memory *) c_ptr;
    c_ptr += sizeof(ocp_nlp_zo_sqp_memory);

    align_char_to(8, &c_ptr);

    // nlp mem
    mem->nlp_mem = ocp_nlp_memory_assign(config, dims, nlp_opts, c_ptr);
    c_ptr += ocp_nlp_memory_calculate_size(config, dims, nlp_opts);

    // stat
    mem->stat = (double *) c_ptr;
    mem->stat_m = opts->max_outer_iter+1;
    mem->stat_n = 6;
    if (opts->ext_qp_res)
        mem->stat_n += 4;
    c_ptr += mem->stat_m*mem->stat_n*sizeof(double);

    mem->status = ACADOS_READY;

    align_char_to(8, &c_ptr);

    assert((char *) raw_memory + ocp_nlp_zo_sqp_memory_calculate_size(config, dims, opts) >= c_ptr);

    return mem;
}



/************************************************
 * workspace
 ************************************************/

acados_size_t ocp_nlp_zo_sqp_workspace_calculate_size(void *config_, void *dims_, void *opts_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;

    acados_size_t size = 0;

    // sqp
    size += sizeof(ocp_nlp_zo_sqp_workspace);

    // nlp
    size += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // tmp qp in
    size += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp qp out
    size += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    if (opts->ext_qp_res)
    {
        // qp res
        size += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        size += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    return size;
}



static void ocp_nlp_zo_sqp_cast_workspace(ocp_nlp_config *config, ocp_nlp_dims *dims,
         ocp_nlp_zo_sqp_opts *opts, ocp_nlp_zo_sqp_memory *mem, ocp_nlp_zo_sqp_workspace *work)
{
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    // sqp
    char *c_ptr = (char *) work;
    c_ptr += sizeof(ocp_nlp_zo_sqp_workspace);

    // nlp
    work->nlp_work = ocp_nlp_workspace_assign(config, dims, nlp_opts, nlp_mem, c_ptr);
    c_ptr += ocp_nlp_workspace_calculate_size(config, dims, nlp_opts);

    // tmp qp in
    work->tmp_qp_in = ocp_qp_in_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_in_calculate_size(dims->qp_solver->orig_dims);

    // tmp qp out
    work->tmp_qp_out = ocp_qp_out_assign(dims->qp_solver->orig_dims, c_ptr);
    c_ptr += ocp_qp_out_calculate_size(dims->qp_solver->orig_dims);

    if (opts->ext_qp_res)
    {
        // qp res
        work->qp_res = ocp_qp_res_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_calculate_size(dims->qp_solver->orig_dims);

        // qp res ws
        work->qp_res_ws = ocp_qp_res_workspace_assign(dims->qp_solver->orig_dims, c_ptr);
        c_ptr += ocp_qp_res_workspace_calculate_size(dims->qp_solver->orig_dims);
    }

    assert((char *) work + ocp_nlp_zo_sqp_workspace_calculate_size(config, dims, opts) >= c_ptr);

    return;
}



/************************************************
 * functions
 ************************************************/

int ocp_nlp_zo_sqp(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{

    acados_timer timer0, timer1;

    acados_tic(&timer0);

    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_opts *nlp_opts = opts->nlp_opts;
    ocp_nlp_zo_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_qp_xcond_solver_config *qp_solver = config->qp_solver;

    ocp_nlp_zo_sqp_workspace *work = work_;
    ocp_nlp_zo_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;
    ocp_nlp_out *inner_nlp_out = nlp_work->inner_nlp_out;

    // zero timers
    double total_time = 0.0;
    double tmp_time;
    mem->time_qp_sol = 0.0;
    mem->time_qp_solver_call = 0.0;
    mem->time_qp_xcond = 0.0;
    mem->time_lin = 0.0;
    mem->time_reg = 0.0;
    mem->time_tot = 0.0;
    mem->time_glob = 0.0;

    int N = dims->N;

    int ii;

    int qp_iter = 0;
    int qp_status = 0;

#if defined(ACADOS_WITH_OPENMP)
    // backup number of threads
    int num_threads_bkp = omp_get_num_threads();
    // set number of threads
    omp_set_num_threads(opts->nlp_opts->num_threads);
    #pragma omp parallel
    { // beginning of parallel region
#endif

    // alias to dynamics_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_ux1_ptr(nlp_out->ux+ii+1, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_ux1_ptr(nlp_work->tmp_nlp_out->ux+ii+1, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_pi_ptr(nlp_out->pi+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_tmp_pi_ptr(nlp_work->tmp_nlp_out->pi+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_BAbt_ptr(nlp_mem->qp_in->BAbt+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_dzduxt_ptr(nlp_mem->dzduxt+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_sim_guess_ptr(nlp_mem->sim_guess+ii, nlp_mem->set_sim_guess+ii, nlp_mem->dynamics[ii]);
        config->dynamics[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->dynamics[ii]);
    }

    // alias to cost_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->cost[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->cost[ii]);
        config->cost[ii]->memory_set_Z_ptr(nlp_mem->qp_in->Z+ii, nlp_mem->cost[ii]);
    }
    // alias to constraints_memory
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif
    for (ii = 0; ii <= N; ii++)
    {
        config->constraints[ii]->memory_set_ux_ptr(nlp_out->ux+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_tmp_ux_ptr(nlp_work->tmp_nlp_out->ux+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_lam_ptr(nlp_out->lam+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_tmp_lam_ptr(nlp_work->tmp_nlp_out->lam+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_z_alg_ptr(nlp_mem->z_alg+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_dzdux_tran_ptr(nlp_mem->dzduxt+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_DCt_ptr(nlp_mem->qp_in->DCt+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_RSQrq_ptr(nlp_mem->qp_in->RSQrq+ii, nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxb_ptr(nlp_mem->qp_in->idxb[ii], nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxs_rev_ptr(nlp_mem->qp_in->idxs_rev[ii], nlp_mem->constraints[ii]);
        config->constraints[ii]->memory_set_idxe_ptr(nlp_mem->qp_in->idxe[ii], nlp_mem->constraints[ii]);
    }

    // alias to regularize memory
    config->regularize->memory_set_RSQrq_ptr(dims->regularize, nlp_mem->qp_in->RSQrq, nlp_mem->regularize_mem);
    config->regularize->memory_set_rq_ptr(dims->regularize, nlp_mem->qp_in->rqz, nlp_mem->regularize_mem);
    config->regularize->memory_set_BAbt_ptr(dims->regularize, nlp_mem->qp_in->BAbt, nlp_mem->regularize_mem);
    config->regularize->memory_set_b_ptr(dims->regularize, nlp_mem->qp_in->b, nlp_mem->regularize_mem);
    config->regularize->memory_set_idxb_ptr(dims->regularize, nlp_mem->qp_in->idxb, nlp_mem->regularize_mem);
    config->regularize->memory_set_DCt_ptr(dims->regularize, nlp_mem->qp_in->DCt, nlp_mem->regularize_mem);
    config->regularize->memory_set_ux_ptr(dims->regularize, nlp_mem->qp_out->ux, nlp_mem->regularize_mem);
    config->regularize->memory_set_pi_ptr(dims->regularize, nlp_mem->qp_out->pi, nlp_mem->regularize_mem);
    config->regularize->memory_set_lam_ptr(dims->regularize, nlp_mem->qp_out->lam, nlp_mem->regularize_mem);

    // copy sampling times into dynamics model
#if defined(ACADOS_WITH_OPENMP)
    #pragma omp for
#endif

    // NOTE(oj): this will lead in an error for irk_gnsf, T must be set in precompute;
    //    -> remove here and make sure precompute is called everywhere.
    for (ii = 0; ii < N; ii++)
    {
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii],
                                         nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
    }

#if defined(ACADOS_WITH_OPENMP)
    } // end of parallel region
#endif

    //
    // if (opts->initialize_t_slacks > 0)
    //     ocp_nlp_initialize_t_slacks(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    // initialize QP
    ocp_nlp_initialize_qp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);

    // main sqp loop
    int sqp_outer_iter = 0;
    int sqp_inner_iter = 0;

    nlp_mem->sqp_iter = &sqp_outer_iter;
    
    // store current solution
    int *nv = dims->nv;
    int *nx = dims->nx;
    int *nu = dims->nu;
    int *ni = dims->ni;
    int *nz = dims->nz;

    bool sens_forw = true;
    for (ii=0; ii<N; ii++)
        config->dynamics[ii]->opts_set(config->dynamics[ii], nlp_opts->dynamics[ii],"sens_forw", &sens_forw);

    for (sqp_outer_iter = 0; sqp_outer_iter < opts->max_outer_iter; sqp_outer_iter++)
    {

        // globalize outer iterates
        // double outer_alpha = 1.0; 
        // for (int i = 0; i <= N; i++)
        // {
        //     // step in primal variables
        //     blasfeo_daxpy(nv[i], 1.0-outer_alpha, nlp_out->ux + i, 0, nlp_work->tmp_nlp_out->ux + i, 0, nlp_out->ux + i, 0);
        
        //     // update dual variables
        //     if (i < N)
        //     {
        //         blasfeo_dvecsc(nx[i+1], 1.0-outer_alpha, nlp_out->pi+i, 0);
        //         blasfeo_daxpy(nx[i+1], outer_alpha, nlp_out->pi+i, 0, nlp_work->tmp_nlp_out->pi+i, 0, nlp_out->pi+i, 0);
        //     }

        //     blasfeo_dvecsc(2*ni[i], 1.0-outer_alpha, nlp_out->lam+i, 0);
        //     blasfeo_daxpy(2*ni[i], outer_alpha, nlp_out->lam+i, 0, nlp_work->tmp_nlp_out->lam+i, 0, nlp_out->lam+i, 0);

        //     // update slack values
        //     blasfeo_dvecsc(2*ni[i], 1.0-outer_alpha, nlp_out->t+i, 0);
        //     blasfeo_daxpy(2*ni[i], outer_alpha, nlp_out->t+i, 0, nlp_work->tmp_nlp_out->t+i, 0, nlp_out->t+i, 0);

        // }

        if (sqp_outer_iter > 0) 
        {
            // update outer iterations
            for (int i = 0; i <= N; i++)
                blasfeo_dveccp(nv[i], inner_nlp_out->ux+i, 0, nlp_out->ux+i, 0);
                // TODO(andrea): slacks????

            for (int i = 0; i < N; i++)
                blasfeo_dveccp(nx[i+1], inner_nlp_out->pi+i, 0, nlp_out->pi+i, 0);

            for (int i = 0; i <= N; i++)
            {
                blasfeo_dveccp(2*ni[i], inner_nlp_out->lam+i, 0, nlp_out->lam+i, 0);
                blasfeo_dveccp(2*ni[i], inner_nlp_out->t+i, 0, nlp_out->t+i, 0);
            }
        }

        bool sens_forw = true;
        for (ii=0; ii<N; ii++)
            config->dynamics[ii]->opts_set(config->dynamics[ii], nlp_opts->dynamics[ii],"sens_forw", &sens_forw);
        
        // switch off hotstart
        int warm_start = 2;
        // condense Hessian
        int hess = 1;
        // expand dual solution
        int dual_sol = 1;

#if !HPIPM_ZO
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "warm_start", &warm_start);
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "cond_hess", &hess);
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "cond_dual_sol", &dual_sol);

        ((dense_qp_qpoases_memory *)(nlp_mem->qp_solver_mem->solver_memory))->first_it = 1;
#endif
        
        printf("res before linearization:\n");
#if 0
        // compute nlp residuals
        ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_mem->nlp_res, nlp_mem);

        nlp_out->inf_norm_res = nlp_mem->nlp_res->inf_norm_res_stat;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_eq > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_eq :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_ineq > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_ineq :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_comp > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_comp :
                                    nlp_out->inf_norm_res;

        if (opts->print_level > 0)
        {
            printf("-%i\t%e\t%e\t%e\t%e.\n", sqp_outer_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                    nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                    nlp_mem->nlp_res->inf_norm_res_comp );
        }

        // printf("outer qp:\n");
        // print_ocp_qp_in(nlp_mem->qp_in);
        if (opts->print_level > sqp_outer_iter + 1)
            print_ocp_qp_in(nlp_mem->qp_in);

        // save statistics
        if (sqp_outer_iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*sqp_outer_iter+0] = nlp_mem->nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*sqp_outer_iter+1] = nlp_mem->nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*sqp_outer_iter+2] = nlp_mem->nlp_res->inf_norm_res_ineq;
            mem->stat[mem->stat_n*sqp_outer_iter+3] = nlp_mem->nlp_res->inf_norm_res_comp;
        }

        // exit conditions on residuals
        if ((nlp_mem->nlp_res->inf_norm_res_stat < opts->tol_stat) &
            (nlp_mem->nlp_res->inf_norm_res_eq < opts->tol_eq) &
            (nlp_mem->nlp_res->inf_norm_res_ineq < opts->tol_ineq) &
            (nlp_mem->nlp_res->inf_norm_res_comp < opts->tol_comp) & sqp_outer_iter != 0)
        {
            // save sqp iterations number
            mem->sqp_outer_iter = sqp_outer_iter;
            nlp_out->sqp_iter = sqp_outer_iter;

            // stop timer
            total_time += acados_toc(&timer0);

            // save time
            nlp_out->total_time = total_time;
            mem->time_tot = total_time;

#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            mem->status = ACADOS_SUCCESS;

            if (opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_outer_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                    nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                    nlp_mem->nlp_res->inf_norm_res_comp );
                printf("\n\n");
            }

            printf("  -> SOLVED OUTER PROBLEM.\n\n");
            return mem->status;
        }
#endif
        // linearize NLP and update QP matrices
        acados_tic(&timer1);
        ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
        mem->time_lin += acados_toc(&timer1);

        // update QP rhs for SQP (step prim var, abs dual var)
        ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);



        // regularize Hessian
        acados_tic(&timer1);
        config->regularize->regularize_hessian(config->regularize, dims->regularize,
                                               opts->nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // // (typically) no warm start at first iteration
        // if (sqp_outer_iter == 0 && !opts->warm_start_first_qp)
        // {
        //     int tmp_int = 0;
        //     config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts,
        //                                  "warm_start", &tmp_int);
        // }

        // solve qp
        acados_tic(&timer1);
        qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, nlp_mem->qp_in, nlp_mem->qp_out,
                                        opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
        mem->time_qp_sol += acados_toc(&timer1);

        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
        mem->time_qp_solver_call += tmp_time;
        qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
        mem->time_qp_xcond += tmp_time;

        // compute correct dual solution in case of Hessian regularization
        acados_tic(&timer1);
        config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                             opts->nlp_opts->regularize, nlp_mem->regularize_mem);
        mem->time_reg += acados_toc(&timer1);

        // // restore default warm start
        // if (sqp_outer_iter==0)
        // {
        //     config->qp_solver->opts_set(config->qp_solver, opts->nlp_opts->qp_solver_opts,
        //                                 "warm_start", &opts->qp_warm_start);
        // }

        // TODO move into QP solver memory ???
        qp_info *qp_info_;
        ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
        nlp_out->qp_iter = qp_info_->num_iter;
        // printf("\nqp_iter = %d, sqp_outer_iter = %d, max_sqp_outer_iter = %d\n", nlp_out->qp_iter, sqp_outer_iter, opts->max_outer_iter);
        qp_iter = qp_info_->num_iter;

        // save statistics of last qp solver call
        if (sqp_outer_iter+1 < mem->stat_m)
        {
            mem->stat[mem->stat_n*(sqp_outer_iter+1)+4] = qp_status;
            mem->stat[mem->stat_n*(sqp_outer_iter+1)+5] = qp_iter;
        }

        // compute external QP residuals (for debugging)
        if (opts->ext_qp_res)
        {
            ocp_qp_res_compute(nlp_mem->qp_in, nlp_mem->qp_out, work->qp_res, work->qp_res_ws);
            if (sqp_outer_iter+1 < mem->stat_m)
                ocp_qp_res_compute_nrm_inf(work->qp_res, mem->stat+(mem->stat_n*(sqp_outer_iter+1)+6));
        }


        if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
        {
            // print_ocp_qp_in(nlp_mem->qp_in);
            if (opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_outer_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                    nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                    nlp_mem->nlp_res->inf_norm_res_comp );
                printf("\n\n");
            }

            // save sqp iterations number
            mem->sqp_outer_iter = sqp_outer_iter;
            nlp_out->sqp_iter = sqp_outer_iter;

            // stop timer
            total_time += acados_toc(&timer0);

            // save time
            mem->time_tot = total_time;
            nlp_out->total_time = total_time;
#ifndef ACADOS_SILENT
            printf("QP solver returned error status %d in iteration %d\n", qp_status, sqp_outer_iter);
#endif
#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif

            if (opts->print_level > 1)
            {
                printf("\n Failed to solve the following QP:\n");
                if (opts->print_level > sqp_outer_iter + 1)
                    print_ocp_qp_in(nlp_mem->qp_in);
            }

            mem->status = ACADOS_QP_FAILURE;
            return mem->status;
        }

        // globalization
        acados_tic(&timer1);
        double alpha = ocp_nlp_line_search(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work);
        mem->time_glob += acados_toc(&timer1);

        // update variables
        printf("alpha = %f\n", alpha);
        ocp_nlp_update_variables_sqp(config, dims, nlp_in, nlp_out, nlp_opts, nlp_mem, nlp_work, alpha);
        printf("res after linearization and QP step:\n");
#if 0
        // compute nlp residuals
        ocp_nlp_res_compute(dims, nlp_in, nlp_out, nlp_mem->nlp_res, nlp_mem);

        nlp_out->inf_norm_res = nlp_mem->nlp_res->inf_norm_res_stat;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_eq > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_eq :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_ineq > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_ineq :
                                    nlp_out->inf_norm_res;
        nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_comp > nlp_out->inf_norm_res) ?
                                    nlp_mem->nlp_res->inf_norm_res_comp :
                                    nlp_out->inf_norm_res;

        if (opts->print_level > 0)
        {
            printf("-%i\t%e\t%e\t%e\t%e.\n", sqp_outer_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                    nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                    nlp_mem->nlp_res->inf_norm_res_comp );
        }

        // printf("outer qp:\n");
        // print_ocp_qp_in(nlp_mem->qp_in);
        if (opts->print_level > sqp_outer_iter + 1)
            print_ocp_qp_in(nlp_mem->qp_in);

        // save statistics
        if (sqp_outer_iter < mem->stat_m)
        {
            mem->stat[mem->stat_n*sqp_outer_iter+0] = nlp_mem->nlp_res->inf_norm_res_stat;
            mem->stat[mem->stat_n*sqp_outer_iter+1] = nlp_mem->nlp_res->inf_norm_res_eq;
            mem->stat[mem->stat_n*sqp_outer_iter+2] = nlp_mem->nlp_res->inf_norm_res_ineq;
            mem->stat[mem->stat_n*sqp_outer_iter+3] = nlp_mem->nlp_res->inf_norm_res_comp;
        }

        // exit conditions on residuals
        if ((nlp_mem->nlp_res->inf_norm_res_stat < opts->tol_stat) &
            (nlp_mem->nlp_res->inf_norm_res_eq < opts->tol_eq) &
            (nlp_mem->nlp_res->inf_norm_res_ineq < opts->tol_ineq) &
            (nlp_mem->nlp_res->inf_norm_res_comp < opts->tol_comp) & sqp_outer_iter != 0)
        {
            // save sqp iterations number
            mem->sqp_outer_iter = sqp_outer_iter;
            nlp_out->sqp_iter = sqp_outer_iter;

            // stop timer
            total_time += acados_toc(&timer0);

            // save time
            nlp_out->total_time = total_time;
            mem->time_tot = total_time;

#if defined(ACADOS_WITH_OPENMP)
            // restore number of threads
            omp_set_num_threads(num_threads_bkp);
#endif
            mem->status = ACADOS_SUCCESS;

            if (opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_outer_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                    nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                    nlp_mem->nlp_res->inf_norm_res_comp );
                printf("\n\n");
            }

            printf("  -> SOLVED OUTER PROBLEM.\n\n");
            return mem->status;
        }
#endif
        // ocp_nlp_dims_print(nlp_out->dims);
        // ocp_nlp_out_print(dims, nlp_out);
        // exit(1);

        // init inner iterations
        for (int i = 0; i <= N; i++)
            blasfeo_dveccp(nv[i], nlp_out->ux+i, 0, inner_nlp_out->ux+i, 0);
            // TODO(andrea): slacks????

        for (int i = 0; i < N; i++)
            blasfeo_dveccp(nx[i+1], nlp_out->pi+i, 0, inner_nlp_out->pi+i, 0);

        for (int i = 0; i <= N; i++)
        {
            blasfeo_dveccp(2*ni[i], nlp_out->lam+i, 0, inner_nlp_out->lam+i, 0);
            blasfeo_dveccp(2*ni[i], nlp_out->t+i, 0, inner_nlp_out->t+i, 0);
        }

        sens_forw = true;
        for (ii=0; ii < N; ii++)
            config->dynamics[ii]->opts_set(config->dynamics[ii], nlp_opts->dynamics[ii],"sens_forw", &sens_forw);

        // switch off hotstart
        warm_start = 2;
        // do not condense Hessian
        hess = 0;
        // do not expand dual solution
        dual_sol = 0;

#if !HPIPM_ZO
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "warm_start", &warm_start);
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "cond_hess", &hess);
        qp_solver->opts_set(qp_solver, opts->nlp_opts->qp_solver_opts, "cond_dual_sol", &dual_sol);
#endif

        for (sqp_inner_iter = 0; sqp_inner_iter < opts->max_inner_iter; sqp_inner_iter++)
        {

            // print_ocp_qp_out(nlp_mem->qp_out);
            // print_ocp_qp_in(nlp_mem->qp_in);
            // linearize NLP and update QP matrices
            acados_tic(&timer1);
            ocp_nlp_approximate_qp_matrices(config, dims, nlp_in, inner_nlp_out, nlp_opts, nlp_mem, nlp_work);
            mem->time_lin += acados_toc(&timer1);

            // update QP rhs for SQP (step prim var, abs dual var)
            ocp_nlp_approximate_qp_vectors_sqp(config, dims, nlp_in, inner_nlp_out, nlp_opts, nlp_mem, nlp_work);


            // regularize Hessian
            acados_tic(&timer1);
            config->regularize->regularize_hessian(config->regularize, dims->regularize,
                                                   opts->nlp_opts->regularize, nlp_mem->regularize_mem);
            mem->time_reg += acados_toc(&timer1);

            // solve qp
            acados_tic(&timer1);
            qp_status = qp_solver->evaluate(qp_solver, dims->qp_solver, nlp_mem->qp_in, nlp_mem->qp_out,
                                            opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);
            mem->time_qp_sol += acados_toc(&timer1);

            qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_solver_call", &tmp_time);
            mem->time_qp_solver_call += tmp_time;
            qp_solver->memory_get(qp_solver, nlp_mem->qp_solver_mem, "time_qp_xcond", &tmp_time);
            mem->time_qp_xcond += tmp_time;

            // compute correct dual solution in case of Hessian regularization
            acados_tic(&timer1);
            config->regularize->correct_dual_sol(config->regularize, dims->regularize,
                                                 opts->nlp_opts->regularize, nlp_mem->regularize_mem);
            mem->time_reg += acados_toc(&timer1);

            // TODO move into QP solver memory ???
            qp_info *qp_info_;
            ocp_qp_out_get(nlp_mem->qp_out, "qp_info", &qp_info_);
            inner_nlp_out->qp_iter = qp_info_->num_iter;
            // printf("\nqp_iter = %d, sqp_iter = %d, max_sqp_iter = %d\n", nlp_out->qp_iter, sqp_iter, opts->max_outer_iter);
            qp_iter = qp_info_->num_iter;


            // printf("inner qp:\n");
            // print_ocp_qp_in(nlp_mem->qp_in);

            if ((qp_status!=ACADOS_SUCCESS) & (qp_status!=ACADOS_MAXITER))
            {
                // print_ocp_qp_in(nlp_mem->qp_in);
                if (opts->print_level > 0)
                {
                    printf("%i\t%e\t%e\t%e\t%e.\n", sqp_inner_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                        nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                        nlp_mem->nlp_res->inf_norm_res_comp );
                    printf("\n\n");
                }

                // stop timer
                total_time += acados_toc(&timer0);

                // save time
                mem->time_tot = total_time;
                inner_nlp_out->total_time = total_time;
#ifndef ACADOS_SILENT
                printf("QP solver returned error status %d in iteration %d\n", qp_status, sqp_inner_iter);
#endif
#if defined(ACADOS_WITH_OPENMP)
                // restore number of threads
                omp_set_num_threads(num_threads_bkp);
#endif

                mem->status = ACADOS_QP_FAILURE;
                return mem->status;
            }

            // globalization
            acados_tic(&timer1);
            // double alpha = ocp_nlp_line_search(config, dims, nlp_in, inner_nlp_out, nlp_opts, nlp_mem, nlp_work);
            alpha = 1.0;
            mem->time_glob += acados_toc(&timer1);

            // update variables
            ocp_nlp_update_variables_sqp(config, dims, nlp_in, inner_nlp_out, nlp_opts, nlp_mem, nlp_work, alpha);
            // ocp_nlp_dims_print(inner_nlp_out->dims);
            // ocp_nlp_out_print(dims, inner_nlp_out);
            // exit(1);

#if 0
            // compute nlp residuals
            ocp_nlp_res_compute(dims, nlp_in, inner_nlp_out, nlp_mem->nlp_res, nlp_mem);

            inner_nlp_out->inf_norm_res = nlp_mem->nlp_res->inf_norm_res_stat;
            inner_nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_eq > inner_nlp_out->inf_norm_res) ?
                                        nlp_mem->nlp_res->inf_norm_res_eq :
                                        inner_nlp_out->inf_norm_res;
            inner_nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_ineq > inner_nlp_out->inf_norm_res) ?
                                        nlp_mem->nlp_res->inf_norm_res_ineq :
                                        inner_nlp_out->inf_norm_res;
            inner_nlp_out->inf_norm_res = (nlp_mem->nlp_res->inf_norm_res_comp > inner_nlp_out->inf_norm_res) ?
                                        nlp_mem->nlp_res->inf_norm_res_comp :
                                        inner_nlp_out->inf_norm_res;

            if (opts->print_level > 0)
            {
                printf("%i\t%e\t%e\t%e\t%e.\n", sqp_inner_iter, nlp_mem->nlp_res->inf_norm_res_stat,
                        nlp_mem->nlp_res->inf_norm_res_eq, nlp_mem->nlp_res->inf_norm_res_ineq,
                        nlp_mem->nlp_res->inf_norm_res_comp );
            }
            
            // exit conditions on residuals
            if ((nlp_mem->nlp_res->inf_norm_res_stat < opts->tol_stat) &
                (nlp_mem->nlp_res->inf_norm_res_eq < opts->tol_eq) &
                (nlp_mem->nlp_res->inf_norm_res_ineq < opts->tol_ineq) &
                (nlp_mem->nlp_res->inf_norm_res_comp < opts->tol_comp))
            {
                // save sqp iterations number

                printf("SOLVED INNER PROBLEM.\n");
                break; 
            }
#else

            // terminate based on primal-step for now

            double temp = 0.0;
            double step_norm = 0.0;
            for (int i=0; i<=N;i++){
                // step in primal variables
                blasfeo_dvecnrm_inf(nv[i], nlp_mem->qp_out->ux + i, 0, &temp);
                step_norm += temp*temp;
            }

            if (opts->print_level > 0)
                printf("%i\t%e.\n", sqp_inner_iter, step_norm);

            if (step_norm < opts->tol_stat && sqp_inner_iter == 0)
            {
                if (opts->print_level > 0)
                    printf("  ->  SOLVED OUTER PROBLEM.\n\n");

                mem->status = ACADOS_SUCCESS;
                return mem->status;
            }

            if (step_norm < 1e-10) {
                {
                    if (opts->print_level > 0)
                        printf("SOLVED INNER PROBLEM.\n");
                    break; 
                }
            }
#endif
        }
    }


    // stop timer
    total_time += acados_toc(&timer0);

    if (opts->print_level > 0)
        printf("\n\n");

    // ocp_nlp_out_print(inner_nlp_out);

    // save sqp iterations number
    mem->sqp_outer_iter = sqp_outer_iter;
    inner_nlp_out->sqp_iter = sqp_outer_iter;

    // save time
    mem->time_tot = total_time;
    inner_nlp_out->total_time = total_time;

    // maximum number of iterations reached
#if defined(ACADOS_WITH_OPENMP)
    // restore number of threads
    omp_set_num_threads(num_threads_bkp);
#endif
    // do not return ACADOS_MAXITER since we are interested in using 
    // this solver with a limited number of iterations
    mem->status = ACADOS_MAXITER;
    // mem->status = ACADOS_SUCCESS;
#ifndef ACADOS_SILENT
    printf("\n ocp_nlp_zo_sqp: maximum iterations reached\n");
#endif

    return mem->status;
}



int ocp_nlp_zo_sqp_precompute(void *config_, void *dims_, void *nlp_in_, void *nlp_out_,
                void *opts_, void *mem_, void *work_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_zo_sqp_memory *mem = mem_;
    ocp_nlp_in *nlp_in = nlp_in_;
    // ocp_nlp_out *nlp_out = nlp_out_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;

    ocp_nlp_zo_sqp_workspace *work = work_;
    ocp_nlp_zo_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    int N = dims->N;
    int status = ACADOS_SUCCESS;

    int ii;

    // TODO(all) add flag to enable/disable checks
    for (ii = 0; ii <= N; ii++)
    {
        int module_val;
        config->constraints[ii]->dims_get(config->constraints[ii], dims->constraints[ii], "ns", &module_val);
        if (dims->ns[ii] != module_val)
        {
            printf("ocp_nlp_zo_sqp_precompute: inconsistent dimension ns for stage %d with constraint module, got %d, module: %d.",
                   ii, dims->ns[ii], module_val);
            exit(1);
        }
    }

    // precompute
    for (ii = 0; ii < N; ii++)
    {
        // set T
        config->dynamics[ii]->model_set(config->dynamics[ii], dims->dynamics[ii],
                                        nlp_in->dynamics[ii], "T", nlp_in->Ts+ii);
        // dynamics precompute
        status = config->dynamics[ii]->precompute(config->dynamics[ii], dims->dynamics[ii],
                                                nlp_in->dynamics[ii], opts->nlp_opts->dynamics[ii],
                                                nlp_mem->dynamics[ii], nlp_work->dynamics[ii]);
        if (status != ACADOS_SUCCESS)
            return status;
    }
    return status;
}



void ocp_nlp_zo_sqp_eval_param_sens(void *config_, void *dims_, void *opts_, void *mem_, void *work_,
                                 char *field, int stage, int index, void *sens_nlp_out_)
{
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;
    ocp_nlp_zo_sqp_memory *mem = mem_;
    ocp_nlp_memory *nlp_mem = mem->nlp_mem;
    ocp_nlp_out *sens_nlp_out = sens_nlp_out_;

    ocp_nlp_zo_sqp_workspace *work = work_;
    ocp_nlp_zo_sqp_cast_workspace(config, dims, opts, mem, work);
    ocp_nlp_workspace *nlp_work = work->nlp_work;

    d_ocp_qp_copy_all(nlp_mem->qp_in, work->tmp_qp_in);
    d_ocp_qp_set_rhs_zero(work->tmp_qp_in);

    double one = 1.0;

    if ((!strcmp("ex", field)) & (stage==0))
    {
        d_ocp_qp_set_el("lbx", stage, index, &one, work->tmp_qp_in);
        d_ocp_qp_set_el("ubx", stage, index, &one, work->tmp_qp_in);

//        d_ocp_qp_print(work->tmp_qp_in->dim, work->tmp_qp_in);

        config->qp_solver->eval_sens(config->qp_solver, dims->qp_solver, work->tmp_qp_in, work->tmp_qp_out,
                               opts->nlp_opts->qp_solver_opts, nlp_mem->qp_solver_mem, nlp_work->qp_work);

//        d_ocp_qp_sol_print(work->tmp_qp_out->dim, work->tmp_qp_out);
//        exit(1);
        
        /* copy tmp_qp_out into sens_nlp_out */

        int i;

        int N = dims->N;
        int *nv = dims->nv;
        int *nx = dims->nx;
        // int *nu = dims->nu;
        int *ni = dims->ni;
        // int *nz = dims->nz;

        for (i = 0; i <= N; i++)
        {
            blasfeo_dveccp(nv[i], work->tmp_qp_out->ux + i, 0, sens_nlp_out->ux + i, 0);

            if (i < N)
                blasfeo_dveccp(nx[i + 1], work->tmp_qp_out->pi + i, 0, sens_nlp_out->pi + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->lam + i, 0, sens_nlp_out->lam + i, 0);

            blasfeo_dveccp(2 * ni[i], work->tmp_qp_out->t + i, 0, sens_nlp_out->t + i, 0);

        }

    }
    else
    {
        printf("\nerror: field %s at stage %d not available in ocp_nlp_zo_sqp_eval_param_sens\n", field, stage);
        exit(1);
    }

    return;
}



// TODO rename memory_get ???
void ocp_nlp_zo_sqp_get(void *config_, void *dims_, void *mem_, const char *field, void *return_value_)
{
    ocp_nlp_config *config = config_;
    ocp_nlp_dims *dims = dims_;
    ocp_nlp_zo_sqp_memory *mem = mem_;

    if (!strcmp("sqp_outer_iter", field))
    {
        int *value = return_value_;
        *value = mem->sqp_outer_iter;
    }
    if (!strcmp("sqp_inner_iter", field))
    {
        int *value = return_value_;
        *value = mem->sqp_inner_iter;
    }
    else if (!strcmp("status", field))
    {
        int *value = return_value_;
        *value = mem->status;
    }
    else if (!strcmp("time_tot", field) || !strcmp("tot_time", field))
    {
        double *value = return_value_;
        *value = mem->time_tot;
    }
    else if (!strcmp("time_qp_sol", field) || !strcmp("time_qp", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_sol;
    }
    else if (!strcmp("time_qp_solver", field) || !strcmp("time_qp_solver_call", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_solver_call;
    }
    else if (!strcmp("time_qp_xcond", field))
    {
        double *value = return_value_;
        *value = mem->time_qp_xcond;
    }
    else if (!strcmp("time_lin", field))
    {
        double *value = return_value_;
        *value = mem->time_lin;
    }
    else if (!strcmp("time_reg", field))
    {
        double *value = return_value_;
        *value = mem->time_reg;
    }
    else if (!strcmp("time_glob", field))
    {
        double *value = return_value_;
        *value = mem->time_glob;
    }
    else if (!strcmp("time_sim", field) || !strcmp("time_sim_ad", field) || !strcmp("time_sim_la", field))
    {
        double tmp = 0.0;
        double *ptr = return_value_;
        int N = dims->N;
        int ii;
        for (ii=0; ii<N; ii++)
        {
            config->dynamics[ii]->memory_get(config->dynamics[ii], dims->dynamics[ii], mem->nlp_mem->dynamics[ii], field, &tmp);
            *ptr += tmp;
        }
    }
    else if (!strcmp("stat", field))
    {
        double **value = return_value_;
        *value = mem->stat;
    }
    else if (!strcmp("statistics", field))
    {
        int n_row = mem->stat_m < mem->sqp_outer_iter+1 ? mem->stat_m : mem->sqp_outer_iter+1;
        double *value = return_value_;
        for (int ii=0; ii<n_row; ii++)
        {
            value[ii+0] = ii;
            for (int jj=0; jj<mem->stat_n; jj++)
                value[ii+(jj+1)*n_row] = mem->stat[jj+ii*mem->stat_n];
        }
    }
    else if (!strcmp("stat_m", field))
    {
        int *value = return_value_;
        *value = mem->stat_m;
    }
    else if (!strcmp("stat_n", field))
    {
        int *value = return_value_;
        *value = mem->stat_n;
    }
    else if (!strcmp("nlp_mem", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem;
    }
    else if (!strcmp("qp_xcond_dims", field))
    {
        void **value = return_value_;
        *value = dims->qp_solver->xcond_dims;
    }
    else if (!strcmp("nlp_res", field))
    {
        ocp_nlp_res **value = return_value_;
        *value = mem->nlp_mem->nlp_res;
    }
    else if (!strcmp("qp_xcond_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_in;
    }
    else if (!strcmp("qp_xcond_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_solver_mem->xcond_qp_out;
    }
    else if (!strcmp("qp_in", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_in;
    }
    else if (!strcmp("qp_out", field))
    {
        void **value = return_value_;
        *value = mem->nlp_mem->qp_out;
    }
    else if (!strcmp("qp_iter", field))
    {
        config->qp_solver->memory_get(config->qp_solver,
            mem->nlp_mem->qp_solver_mem, "iter", return_value_);
    }
    else if (!strcmp("res_stat", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_stat;
    }
    else if (!strcmp("res_eq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_eq;
    }
    else if (!strcmp("res_ineq", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_ineq;
    }
    else if (!strcmp("res_comp", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->nlp_res->inf_norm_res_comp;
    }
    else if (!strcmp("cost_value", field))
    {
        double *value = return_value_;
        *value = mem->nlp_mem->cost_value;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_zo_sqp_get\n", field);
        exit(1);
    }
}



void ocp_nlp_zo_sqp_opts_get(void *config_, void *dims_, void *opts_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_opts *opts = opts_;

    if (!strcmp("nlp_opts", field))
    {
        void **value = return_value_;
        *value = opts->nlp_opts;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_zo_sqp_opts_get\n", field);
        exit(1);
    }
}


void ocp_nlp_zo_sqp_work_get(void *config_, void *dims_, void *work_,
                          const char *field, void *return_value_)
{
    // ocp_nlp_config *config = config_;
    ocp_nlp_zo_sqp_workspace *work = work_;

    if (!strcmp("nlp_work", field))
    {
        void **value = return_value_;
        *value = work->nlp_work;
    }
    else
    {
        printf("\nerror: field %s not available in ocp_nlp_zo_sqp_work_get\n", field);
        exit(1);
    }
}


void ocp_nlp_zo_sqp_config_initialize_default(void *config_)
{
    ocp_nlp_config *config = (ocp_nlp_config *) config_;

    config->opts_calculate_size = &ocp_nlp_zo_sqp_opts_calculate_size;
    config->opts_assign = &ocp_nlp_zo_sqp_opts_assign;
    config->opts_initialize_default = &ocp_nlp_zo_sqp_opts_initialize_default;
    config->opts_update = &ocp_nlp_zo_sqp_opts_update;
    config->opts_set = &ocp_nlp_zo_sqp_opts_set;
    config->opts_set_at_stage = &ocp_nlp_zo_sqp_opts_set_at_stage;
    config->memory_calculate_size = &ocp_nlp_zo_sqp_memory_calculate_size;
    config->memory_assign = &ocp_nlp_zo_sqp_memory_assign;
    config->workspace_calculate_size = &ocp_nlp_zo_sqp_workspace_calculate_size;
    config->evaluate = &ocp_nlp_zo_sqp;
    config->eval_param_sens = &ocp_nlp_zo_sqp_eval_param_sens;
    config->config_initialize_default = &ocp_nlp_zo_sqp_config_initialize_default;
    config->precompute = &ocp_nlp_zo_sqp_precompute;
    config->get = &ocp_nlp_zo_sqp_get;
    config->opts_get = &ocp_nlp_zo_sqp_opts_get;
    config->work_get = &ocp_nlp_zo_sqp_work_get;

    return;
}
