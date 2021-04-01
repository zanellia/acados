from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from reduced_model import reduced_model
import numpy as np
import os


def reduced_settings(N, Ts, solver_type="SQP_RTI", qp_solver="PARTIAL_CONDENSING_HPIPM", soft_constr=False, print_level=0):
    # create render arguments
    ocp = AcadosOcp()

    # export model
    model = reduced_model()

    # define acados ODE
    model_ac = AcadosModel()
    model_ac.f_impl_expr = model.f_impl_expr
    model_ac.f_expl_expr = model.f_expl_expr
    model_ac.x = model.x
    model_ac.xdot = model.xdot
    model_ac.u = model.u
    model_ac.z = model.z
    model_ac.name = model.name
    ocp.model = model_ac

    # set dimensions
    nx = model.x.size()[0]
    nu = model.u.size()[0]

    ocp.dims.nx = nx
    ocp.dims.np = 0
    ocp.dims.nbx = nx
    ocp.dims.nbu = nu
    ocp.dims.nu = nu
    ocp.dims.N = N
    ocp.dims.nh = 0
    ocp.dims.nsh = 0
    ocp.dims.ns = nx + nu

    if False:
        ny = 0
        ny_e = 0
        ocp.dims.ny = ny
        ocp.dims.ny_e = ny_e

        ocp.cost.cost_type = "EXTERNAL"
        ocp.cost.cost_type_e = "EXTERNAL"
        ocp.model.cost_expr_ext_cost = model.cost_expr_ext_cost
        ocp.model.cost_expr_ext_cost_e = 0
    else:
        ny = 1
        ny_e = 1
        ocp.dims.ny = ny
        ocp.dims.ny_e = ny_e

        ocp.cost.cost_type = 'LINEAR_LS'
        ocp.cost.cost_type_e = 'LINEAR_LS'
        ocp.cost.W = np.diag([1])
        ocp.cost.W_e = np.diag([1])

        Vx = np.zeros((ny, nx))
        Vx[0, 0] = 1
        ocp.cost.Vx = Vx
        ocp.cost.Vx_e = Vx

        Vu = np.zeros((ny, nu))
        ocp.cost.Vu = Vu

        Y = np.ones((ny,))
        ocp.cost.yref = Y
        ocp.cost.yref_e = Y
    ####

    if soft_constr:
        ocp.dims.nsbx = nx
        ocp.dims.nsbu = nu

        ocp.constraints.lsbx = np.zeros(ocp.dims.nsbx)
        ocp.constraints.usbx = np.zeros(ocp.dims.nsbx)
        ocp.constraints.idxsbx = np.array([0])

        ocp.constraints.lsbu = np.zeros(ocp.dims.nsbu)
        ocp.constraints.usbu = np.zeros(ocp.dims.nsbu)
        ocp.constraints.idxsbu = np.array([0])

        ocp.cost.zl = 100*np.ones(2)
        ocp.cost.zu = 100*np.ones(2)
        ocp.cost.Zl = np.zeros(2)
        ocp.cost.Zu = np.zeros(2)
    else:
        ocp.dims.nsbx = nx
        ocp.dims.nsbu = nu

    # setting constraints

    ocp.constraints.lbx = np.array(     # lower bounds on x
        [
            model.x_min
        ]
    )

    ocp.constraints.ubx = np.array(     # upper bounds on x
        [
            model.x_max
        ]
    )

    ocp.constraints.idxbx = np.array([0])   # indices of bounds on x

    ocp.constraints.lbu = np.array(     # lower bounds on u
        [
            model.v_min
        ]
    )
    ocp.constraints.ubu = np.array(     # upper bounds on u
        [
            model.v_max
        ]
    )
    ocp.constraints.idxbu = np.array([0])  # indices of bounds on u

    # set intial condition
    ocp.constraints.x0 = np.zeros(nx)
    ocp.constraints.idxbx_0 = np.array([0])
    ocp.parameter_values = np.zeros(ocp.dims.np)
    # set QP solver and integration
    ocp.solver_options.Tsim = Ts
    ocp.solver_options.tf = Ts * N
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.nlp_solver_type = solver_type
    ocp.solver_options.hessian_approx = "GAUSS_NEWTON"
    ocp.solver_options.integrator_type = "ERK"
    ocp.solver_options.nlp_solver_step_length = 1.0
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.sim_method_num_steps = 3
    ocp.solver_options.print_level = print_level

    ocp.solver_options.tol = 1e-4

    acados_source_path = os.environ['ACADOS_SOURCE_DIR']
    ocp.acados_include_path = acados_source_path + '/include'
    ocp.acados_lib_path = acados_source_path + '/lib'

    # create solver with agent specific code files
    filename = "acados_mpcc_pacejka.json"
    acados_solver = AcadosOcpSolver(ocp, json_file=filename)
    acados_integrator = AcadosSimSolver(ocp, json_file=filename)

    return model, acados_solver, acados_integrator, nx, nu, ocp.dims.np
