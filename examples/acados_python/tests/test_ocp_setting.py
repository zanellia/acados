#
# Copyright 2019 Gianluca Frison, Dimitris Kouzoupis, Robin Verschueren,
# Andrea Zanelli, Niels van Duijkeren, Jonathan Frey, Tommaso Sartor,
# Branimir Novoselnik, Rien Quirynen, Rezart Qelibari, Dang Doan,
# Jonas Koenemann, Yutao Chen, Tobias Schöls, Jonas Schlagenhauf, Moritz Diehl
#
# This file is part of acados.
#
# The 2-Clause BSD License
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# 1. Redistributions of source code must retain the above copyright notice,
# this list of conditions and the following disclaimer.
#
# 2. Redistributions in binary form must reproduce the above copyright notice,
# this list of conditions and the following disclaimer in the documentation
# and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.;
#

import sys
sys.path.insert(0, '../getting_started/common')

from acados_template import *
from export_pendulum_ode_model import export_pendulum_ode_model
import numpy as np
import scipy.linalg
from ctypes import *
import argparse

# set to 'True' to generate test data
GENERATE_DATA = False

LOCAL_TEST = False
TEST_TOL = 1e-8

if LOCAL_TEST is True:
    COST_MODULE = 'LS'
    COST_MODULE_N = 'LS'
    SOLVER_TYPE = 'SQP_RTI'
    QP_SOLVER = 'FULL_CONDENSING_QPOASES'
    INTEGRATOR_TYPE = 'IRK'
else:
    parser = argparse.ArgumentParser(description='test Python interface on pendulum example.')
    parser.add_argument('--COST_MODULE', dest='COST_MODULE',
                        default='LS',
                        help='COST_MODULE: linear least-squares (LS) or nonlinear \
                                least-squares (NLS) (default: LS), external (EXTERNAL)')

    parser.add_argument('--COST_MODULE_N', dest='COST_MODULE_N',
                        default='LS',
                        help='COST_MODULE_N: linear least-squares (LS) or nonlinear \
                                least-squares (NLS) (default: LS), external (EXTERNAL)')

    parser.add_argument('--QP_SOLVER', dest='QP_SOLVER',
                        default='PARTIAL_CONDENSING_HPIPM',
                        help='QP_SOLVER: PARTIAL_CONDENSING_HPIPM, FULL_CONDENSING_HPIPM, ' \
                                'FULL_CONDENSING_HPIPM (default: PARTIAL_CONDENSING_HPIPM)')

    parser.add_argument('--INTEGRATOR_TYPE', dest='INTEGRATOR_TYPE',
                        default='ERK',
                        help='INTEGRATOR_TYPE: explicit (ERK) or implicit (IRK) or GNSF-IRK (GNSF) ' \
                                ' Runge-Kutta (default: ERK)')

    parser.add_argument('--SOLVER_TYPE', dest='SOLVER_TYPE',
                        default='SQP_RTI',
                        help='SOLVER_TYPE: (full step) sequential quadratic programming (SQP) or ' \
                                ' real-time iteration (SQP-RTI) (default: SQP-RTI)')


    args = parser.parse_args()

    COST_MODULE = args.COST_MODULE
    COST_MODULE_values = ['LS', 'NLS', 'EXTERNAL']
    if COST_MODULE not in COST_MODULE_values:
        raise Exception('Invalid unit test value {} for parameter COST_MODULE. Possible values are' \
                ' {}. Exiting.'.format(COST_MODULE, COST_MODULE_values))

    COST_MODULE_N = args.COST_MODULE_N
    COST_MODULE_N_values = ['LS', 'NLS', 'EXTERNAL']
    if COST_MODULE_N not in COST_MODULE_N_values:
        raise Exception('Invalid unit test value {} for parameter COST_MODULE_N. Possible values are' \
                ' {}. Exiting.'.format(COST_MODULE_N, COST_MODULE_N_values))

    QP_SOLVER = args.QP_SOLVER
    QP_SOLVER_values = ['PARTIAL_CONDENSING_HPIPM', 'FULL_CONDENSING_HPIPM', 'FULL_CONDENSING_QPOASES']
    if QP_SOLVER not in QP_SOLVER_values:
        raise Exception('Invalid unit test value {} for parameter QP_SOLVER. Possible values are' \
                ' {}. Exiting.'.format(QP_SOLVER, QP_SOLVER_values))

    INTEGRATOR_TYPE = args.INTEGRATOR_TYPE
    INTEGRATOR_TYPE_values = ['ERK', 'IRK', 'GNSF']
    if INTEGRATOR_TYPE not in INTEGRATOR_TYPE:
        raise Exception('Invalid unit test value {} for parameter INTEGRATOR_TYPE. Possible values are' \
                ' {}. Exiting.'.format(INTEGRATOR_TYPE, INTEGRATOR_TYPE_values))

    SOLVER_TYPE = args.SOLVER_TYPE
    SOLVER_TYPE_values = ['SQP', 'SQP-RTI']
    if SOLVER_TYPE not in SOLVER_TYPE:
        raise Exception('Invalid unit test value {} for parameter SOLVER_TYPE. Possible values are' \
                ' {}. Exiting.'.format(SOLVER_TYPE, SOLVER_TYPE_values))


# print test setting
print("Running test with:\n\tcost module:", COST_MODULE, "\n\tqp solver: ", QP_SOLVER,\
      "\n\tintergrator: ", INTEGRATOR_TYPE, "\n\tsolver: ", SOLVER_TYPE)

# create ocp object to formulate the OCP
ocp = AcadosOcp()

# set model
model = export_pendulum_ode_model()
ocp.model = model

Tf = 1.0
nx = model.x.size()[0]
nu = model.u.size()[0]
ny = nx + nu
ny_e = nx
N = 20

# set dimensions
ocp.dims.nx  = nx
ocp.dims.ny  = ny
ocp.dims.ny_e = ny_e
ocp.dims.nbx = 0
ocp.dims.nbu = nu 
ocp.dims.nu  = nu
ocp.dims.N   = N

# set cost
Q = 2*np.diag([1e3, 1e3, 1e-2, 1e-2])
R = 2*np.diag([1e-2])

ocp.cost.W_e = Q
ocp.cost.W = scipy.linalg.block_diag(Q, R)

x = ocp.model.x
u = ocp.model.u

if COST_MODULE == 'LS':
    ocp.cost.cost_type = 'LINEAR_LS'

    ocp.cost.Vx = np.zeros((ny, nx))
    ocp.cost.Vx[:nx,:nx] = np.eye(nx)

    Vu = np.zeros((ny, nu))
    Vu[4,0] = 1.0
    ocp.cost.Vu = Vu

elif COST_MODULE == 'NLS':
    ocp.cost.cost_type = 'NONLINEAR_LS'
    ocp.model.cost_y_expr = vertcat(x, u)

elif COST_MODULE == 'EXTERNAL':
    ocp.cost.cost_type = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost = vertcat(x, u).T @ ocp.cost.W @ vertcat(x, u)

else:
    raise Exception('Unknown COST_MODULE. Possible values are \'LS\', \'NLS\', \'EXTERNAL\'.')


if COST_MODULE_N == 'LS':
    ocp.cost.cost_type_e = 'LINEAR_LS'
    ocp.cost.Vx_e = np.eye(nx)

elif COST_MODULE_N == 'NLS':
    ocp.cost.cost_type_e = 'NONLINEAR_LS'
    ocp.model.cost_y_expr_e = x

elif COST_MODULE_N == 'EXTERNAL':
    ocp.cost.cost_type_e = 'EXTERNAL'
    ocp.model.cost_expr_ext_cost_e = x.T @ Q @ x

else:
    raise Exception('Unknown COST_MODULE_N. Possible values are \'LS\', \'NLS\', \'EXTERNAL\'.')


ocp.cost.yref  = np.zeros((ny, ))
ocp.cost.yref_e = np.zeros((ny_e, ))

# set constraints
Fmax = 80
ocp.constraints.constr_type = 'BGH'
ocp.constraints.lbu = np.array([-Fmax])
ocp.constraints.ubu = np.array([+Fmax])
ocp.constraints.x0 = np.array([0.0, np.pi, 0.0, 0.0])
ocp.constraints.idxbu = np.array([0])

# set options
ocp.solver_options.qp_solver = QP_SOLVER
ocp.solver_options.hessian_approx = 'GAUSS_NEWTON'
ocp.solver_options.integrator_type = INTEGRATOR_TYPE
ocp.solver_options.sim_method_num_stages = 2
ocp.solver_options.sim_method_num_steps = 5
ocp.solver_options.sim_method_newton_iter = 3

ocp.solver_options.nlp_solver_tol_stat = TEST_TOL
ocp.solver_options.nlp_solver_tol_eq = TEST_TOL
ocp.solver_options.nlp_solver_tol_ineq = TEST_TOL
ocp.solver_options.nlp_solver_tol_comp = TEST_TOL

ocp.solver_options.qp_solver_cond_N = 10
ocp.solver_options.nlp_solver_max_iter = 80
ocp.solver_options.qp_solver_iter_max = 50

# set prediction horizon
ocp.solver_options.tf = Tf
ocp.solver_options.nlp_solver_type = SOLVER_TYPE

if ocp.solver_options.integrator_type == 'GNSF':
    with open('../getting_started/common/' + model.name + '_gnsf_functions.json', 'r') as f:
        gnsf_dict = json.load(f)
    ocp.gnsf_model = gnsf_dict

ocp_solver = AcadosOcpSolver(ocp, json_file = 'acados_ocp.json')

simX = np.ndarray((N+1, nx))
simU = np.ndarray((N, nu))

status = ocp_solver.solve()

if status != 0:
    raise Exception('acados returned status {}. Exiting.'.format(status))

# get solution
for i in range(N):
    simX[i,:] = ocp_solver.get(i, "x")
    simU[i,:] = ocp_solver.get(i, "u")
simX[N,:] = ocp_solver.get(N, "x")

if COST_MODULE in {'LINEAR_LS', 'NONLINEAR_LS'}:
    # update reference
    for j in range(N):
        ocp_solver.cost_set(j, "yref", np.array([0, 0, 0, 0, 0]))
    ocp_solver.cost_set(N, "yref", np.array([0, 0, 0, 0]))

# dump result to JSON file for unit testing
test_file_name = 'test_data/pendulum_ocp_formulations/test_ocp_' + COST_MODULE + '_' + COST_MODULE_N + '_' + QP_SOLVER + '_' + \
            INTEGRATOR_TYPE + '_' + SOLVER_TYPE + '.json'

if GENERATE_DATA:
    with open(test_file_name, 'w') as f:
        json.dump({"simX": simX.tolist(), "simU": simU.tolist()}, f, indent=4, sort_keys=True)
else:
    with open(test_file_name, 'r') as f:
        test_data = json.load(f)
    simX_error = np.linalg.norm(test_data['simX'] - simX)
    simU_error = np.linalg.norm(test_data['simU'] - simU)

    if simX_error > TEST_TOL or simU_error > TEST_TOL:
        raise Exception("Python acados test failure with accuracies" +
                        " {:.2E} and {:.2E} ({:.2E} required)".format(simX_error, simU_error, TEST_TOL) +
                        " on pendulum example! Exiting.\n")
    else:
        print('Python test passed with accuracy {:.2E}'.format(max(simU_error, simX_error)))
