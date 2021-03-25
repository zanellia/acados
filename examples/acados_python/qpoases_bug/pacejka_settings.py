from acados_template import AcadosModel, AcadosOcp, AcadosOcpSolver, AcadosSimSolver
from pacejka_model import pacejka_model
import numpy as np
import os


def pacejka_settings(N, Ts, solver_type="SQP_RTI", qp_solver="PARTIAL_CONDENSING_HPIPM", print_level=0):
    # create render arguments
    ocp = AcadosOcp()

    # export model
    model, constraint = pacejka_model()

    # define acados ODE
    model_ac = AcadosModel()
    model_ac.f_impl_expr = model.f_impl_expr
    model_ac.f_expl_expr = model.f_expl_expr
    model_ac.x = model.x
    model_ac.xdot = model.xdot
    model_ac.u = model.u
    model_ac.z = model.z
    model_ac.p = model.p
    model_ac.name = model.name
    ocp.model = model_ac

    # define constraint
    model_ac.con_h_expr = constraint.expr

    # set dimensions
    nx = model.x.size()[0]
    nu = model.u.size()[0]
    ny = 0
    ny_e = 0

    ocp.dims.nx = nx
    ocp.dims.np = model.p.size()[0]
    ocp.dims.ny = ny
    ocp.dims.ny_e = ny_e
    ocp.dims.nbx = nx
    ocp.dims.nbu = nu
    ocp.dims.nu = nu
    ocp.dims.N = N
    ocp.dims.nh = 1
    ocp.dims.nsh = 1
    ocp.dims.ns = nx + nu + 1

    ocp.dims.nsbx = nx
    ocp.dims.nsbu = nu

    ocp.cost.cost_type = "EXTERNAL"
    ocp.cost.cost_type_e = "EXTERNAL"
    ocp.model.cost_expr_ext_cost = model.cost_expr_ext_cost
    ocp.model.cost_expr_ext_cost_e = 0

    ocp.constraints.lh = np.array([0.0])
    ocp.constraints.uh = np.array([0.23 * 0.23])

    ocp.constraints.lsh = np.zeros(ocp.dims.nsh)
    ocp.constraints.ush = np.zeros(ocp.dims.nsh)
    ocp.constraints.idxsh = np.array([0])

    ocp.constraints.lsbx = np.zeros(ocp.dims.nsbx)
    ocp.constraints.usbx = np.zeros(ocp.dims.nsbx)
    ocp.constraints.idxsbx = np.array([0,1,2,3,4,5,6,7,8])

    ocp.constraints.lsbu = np.zeros(ocp.dims.nsbu)
    ocp.constraints.usbu = np.zeros(ocp.dims.nsbu)
    ocp.constraints.idxsbu = np.array([0,1,2])

    ocp.cost.zl = 100*np.ones(13)
    ocp.cost.zu = 100*np.ones(13)
    ocp.cost.Zl = np.zeros(13)
    ocp.cost.Zu = np.zeros(13)

    ocp.solver_options.levenberg_marquardt = 0.0

    # setting constraints

    ocp.constraints.lbx = np.array(     # lower bounds on x
        [
            model.x_min,
            model.y_min,
            model.yaw_min,
            model.vx_min,
            model.vy_min,
            model.dyaw_min,
            model.T_min,
            model.delta_min,
            model.theta_min
        ]
    )

    ocp.constraints.ubx = np.array(     # upper bounds on x
        [
            model.x_max,
            model.y_max,
            model.yaw_max,
            model.vx_max,
            model.vy_max,
            model.dyaw_max,
            model.T_max,
            model.delta_max,
            model.theta_max
        ]
    )

    ocp.constraints.idxbx = np.array([0, 1, 2, 3, 4, 5 ,6, 7, 8])   # indices of bounds on x

    ocp.constraints.lbu = np.array(     # lower bounds on u
        [
            model.dT_min,
            model.ddelta_min,
            model.dtheta_min
        ]
    )
    ocp.constraints.ubu = np.array(     # upper bounds on u
        [
            model.dT_max,
            model.ddelta_max,
            model.dtheta_max
        ]
    )
    ocp.constraints.idxbu = np.array([0, 1, 2])  # indices of bounds on u

    # set intial condition
    ocp.constraints.x0 = np.zeros(nx)
    ocp.constraints.idxbx_0 = np.array([0, 1, 2, 3, 4, 5, 6, 7, 8])
    ocp.parameter_values = np.zeros(ocp.dims.np)
    # set QP solver and integration
    ocp.solver_options.Tsim = Ts
    ocp.solver_options.tf = Ts * N
    ocp.solver_options.qp_solver = qp_solver
    ocp.solver_options.nlp_solver_type = solver_type
    ocp.solver_options.hessian_approx = "EXACT"
    ocp.solver_options.integrator_type = "ERK"
    ocp.solver_options.nlp_solver_step_length = 1.0
    ocp.solver_options.qp_solver_iter_max = 50
    ocp.solver_options.sim_method_num_stages = 4
    ocp.solver_options.sim_method_num_steps = 3
    ocp.solver_options.print_level = print_level
    ocp.solver_options.qp_solver_iter_max = 1000

    ocp.solver_options.tol = 1e-4
    # ocp.solver_options.levenberg_marquardt = 1e-1

    acados_source_path = os.environ['ACADOS_SOURCE_DIR']
    ocp.acados_include_path = acados_source_path + '/include'
    ocp.acados_lib_path = acados_source_path + '/lib'

    # create solver with agent specific code files
    filename = "acados_mpcc_pacejka.json"
    acados_solver = AcadosOcpSolver(ocp, json_file=filename)
    acados_integrator = AcadosSimSolver(ocp, json_file=filename)

    return model, acados_solver, acados_integrator, nx, nu, ocp.dims.np
