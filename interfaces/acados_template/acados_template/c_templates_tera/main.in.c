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

// standard
#include <stdio.h>
#include <stdlib.h>
// acados
#include "acados/utils/print.h"
#include "acados_c/ocp_nlp_interface.h"
#include "acados_c/external_function_interface.h"
#include "acados_solver_{{model.name}}.h"

// ** global data **
ocp_nlp_in * nlp_in;
ocp_nlp_out * nlp_out;
ocp_nlp_solver * nlp_solver;
void * nlp_opts;
ocp_nlp_plan * nlp_solver_plan;
ocp_nlp_config * nlp_config;
ocp_nlp_dims * nlp_dims;
{% if solver_config.integrator_type == "ERK" %}
{% if dims.np < 1 %}
external_function_casadi * forw_vde_casadi;
{% else %}
external_function_param_casadi * forw_vde_casadi;
{% endif %}
{% if solver_config.hessian_approx == "EXACT" %} 
{% if dims.np < 1 %}
external_function_casadi * hess_vde_casadi;
{% else %}
external_function_param_casadi * hess_vde_casadi;
{% endif %}
{% endif %}
{% else %}
{% if solver_config.integrator_type == "IRK" %}
{% if dims.np < 1 %}
external_function_casadi * impl_dae_fun;
external_function_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_casadi * impl_dae_jac_x_xdot_u_z;
{% else %}
external_function_param_casadi * impl_dae_fun;
external_function_param_casadi * impl_dae_fun_jac_x_xdot_z;
external_function_param_casadi * impl_dae_jac_x_xdot_u_z;
{% endif %}
{% endif %}
{% endif %}
{% if dims.npd > 0 %}
external_function_casadi * p_constraint;
{% endif %}
{% if dims.npd_e > 0 %}
external_function_casadi * p_e_constraint;
{% endif %}
{% if dims.nh > 0 %}
external_function_casadi * h_constraint;
{% endif %}
{% if dims.nh_e > 0 %}
external_function_casadi h_e_constraint;
{% endif %}

int main() {

    int status = 0;
    status = acados_create();

    if (status) { 
        printf("acados_create() returned status %d. Exiting.\n", status); 
        exit(1); }

    // set initial condition
    double x0[{{dims.nx}}];
    {% for item in constraints.x0 %}
    x0[{{ loop.index0 }}] = {{ item }};
    {% endfor %}
    
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "lbx", x0);
    ocp_nlp_constraints_model_set(nlp_config, nlp_dims, nlp_in, 0, "ubx", x0);

    {% if dims.np > 0%}
    double p[{{dims.np}}];
    {% for item in constraints.p %}
    p[{{ loop.index0 }}] = {{ item }};
    {% endfor %}
    {% endif %}
    
    {% if dims.np > 0%}
    {% if solver_config.integrator_type == "IRK" %}
    for (int ii = 0; ii < {{dims.N}}; ii++) {
    impl_dae_fun[ii].set_param(impl_dae_fun+ii, p);
    impl_dae_fun_jac_x_xdot_z[ii].set_param(impl_dae_fun_jac_x_xdot_z+ii, p);
    impl_dae_jac_x_xdot_u_z[ii].set_param(impl_dae_jac_x_xdot_u_z+ii, p);
    }
    {% else %}
    for (int ii = 0; ii < {{dims.N}}; ii++) {
    forw_vde_casadi[ii].set_param(forw_vde_casadi+ii, p);
    }
    {% endif %}
    {% endif %}

    double kkt_norm_inf = 1e12, elapsed_time;

#if 1
    int NTIMINGS = 100;
    double min_time = 1e12;
    for (int ii = 0; ii < NTIMINGS; ii ++) {
        for (int i = 0; i <= nlp_dims->N; ++i) {
            blasfeo_dvecse(nlp_dims->nu[i]+nlp_dims->nx[i], 0.0, nlp_out->ux+i, 0);
        }
        status = acados_solve();
        elapsed_time = nlp_out->total_time;
        if (elapsed_time < min_time) min_time = elapsed_time;
    }
    elapsed_time = min_time;
#else
    status = acados_solve();
#endif
    kkt_norm_inf = nlp_out->inf_norm_res;
    elapsed_time = nlp_out->total_time;
    printf(" iterations %2d | time  %f |  KKT %e\n", nlp_out->sqp_iter, elapsed_time, kkt_norm_inf);

    printf("\n--- solution ---\n");
    ocp_nlp_out_print(nlp_solver->dims, nlp_out);
    if (status) { 
        printf("acados_solve() returned status %d.\n", status); 
    }

    status = acados_free();

    if (status) { 
        printf("acados_free() returned status %d. \n", status); 
    }

    return status;
}
