import math
import numpy as np
from casadi import SX, types, vertcat


def pacejka_model():

    # define structs
    model = types.SimpleNamespace()
    constraint = types.SimpleNamespace()
    model_name = "pacejka_model"

    """ Constraints """
    model.x_min = -2.0
    model.x_max = 2.0

    model.y_min = -2.0
    model.y_max = 2.0

    model.yaw_min = -1000.0
    model.yaw_max = 1000.0

    model.vx_min = 0.0
    model.vx_max = 3.5

    model.vy_min = -2.0
    model.vy_max = 2.0

    model.dyaw_min = -4.0
    model.dyaw_max = 4.0

    model.delta_min = -0.4
    model.delta_max = 0.4

    model.T_min = 0.0
    model.T_max = 1.0

    model.theta_min = 0.0
    model.theta_max = 1000.0

    model.ddelta_min = -2.0
    model.ddelta_max = 2.0

    model.dT_min = -15.0
    model.dT_max = 15.0

    model.dtheta_min = 0.0
    model.dtheta_max = 3.5

    # CasADi - states
    xp = SX.sym("xp")
    yp = SX.sym("yp")
    yaw = SX.sym("yaw")
    vx = SX.sym("vx")
    vy = SX.sym("vy")
    omega = SX.sym("omega")
    T = SX.sym("T")
    delta = SX.sym("delta")
    theta = SX.sym("theta")
    x = vertcat(xp, yp, yaw, vx, vy, omega, T, delta, theta)

    # CasADi - input
    dT = SX.sym("dT")
    ddelta = SX.sym("ddelta")
    dtheta = SX.sym("dtheta")
    u = vertcat(dT, ddelta, dtheta)

    # xdot
    xpdot = SX.sym("xpdot")
    ypdot = SX.sym("ypdot")
    yawdot = SX.sym("yawdot")
    vxdot = SX.sym("vxdot")
    vydot = SX.sym("vydot")
    omegadot = SX.sym("omegadot")
    xdot = vertcat(xpdot, ypdot, yawdot, vxdot, vydot,
                   omegadot, dT, ddelta, dtheta)

    # algebraic variables
    z = vertcat([])

    # define params
    xd = SX.sym("xd")
    yd = SX.sym("yd")
    grad_xd = SX.sym("grad_xd")
    grad_yd = SX.sym("grad_yd")
    theta_hat = SX.sym("theta_hat")
    phi_d = SX.sym("phi_d")
    Q1 = SX.sym("Q1")
    Q2 = SX.sym("Q2")
    R1 = SX.sym("R1")
    R2 = SX.sym("R2")
    R3 = SX.sym("R3")
    q = SX.sym("q")
    lr = SX.sym("lr")
    lf = SX.sym("lf")
    m = SX.sym("m")
    I = SX.sym("I")
    Df = SX.sym("Df")
    Cf = SX.sym("Cf")
    Bf = SX.sym("Bf")
    Dr = SX.sym("Dr")
    Cr = SX.sym("Cr")
    Br = SX.sym("Br")
    Cm1 = SX.sym("Cm1")
    Cm2 = SX.sym("Cm2")
    Cd = SX.sym("Cd")
    Croll = SX.sym("Croll")

    # parameters
    p = vertcat(xd, yd, grad_xd, grad_yd, theta_hat, phi_d, Q1, Q2, R1,
                R2, R3, q, lr, lf, m, I, Df, Cf, Bf, Dr, Cr, Br, Cm1, Cm2, Cd, Croll)

    # dynamics
    Fx = (Cm1 - Cm2 * vx) * T - Cd * vx * vx - Croll
    beta = np.arctan(vy / vx)
    ar = -beta + lr * omega / vx
    af = delta - beta - lf * omega / vx
    Fr = Dr * np.sin(Cr * np.arctan(Br * ar))
    Ff = Df * np.sin(Cf * np.arctan(Bf * af))
    f_expl = vertcat(
        vx * np.cos(yaw) - vy * np.sin(yaw),
        vx * np.sin(yaw) + vy * np.cos(yaw),
        omega,
        1 / m * (Fx - Ff * np.sin(delta) + m * vy * omega),
        1 / m * (Fr + Ff * np.cos(delta) - m * vx * omega),
        1 / I * (Ff * lf * np.cos(delta) - Fr * lr),
        dT,
        ddelta,
        dtheta,
    )

    # cost
    eC = np.sin(phi_d)*(xp-xd-grad_xd*(theta-theta_hat)) - \
        np.cos(phi_d)*(yp-yd-grad_yd*(theta-theta_hat))
    eL = -np.cos(phi_d)*(xp-xd-grad_xd*(theta-theta_hat)) - \
        np.sin(phi_d)*(yp-yd-grad_yd*(theta-theta_hat))

    c_eC = eC*eC*Q1
    c_eL = eL*eL*Q2
    c_theta = -q*theta
    c_dT = dT*dT*R1
    c_ddelta = ddelta*ddelta*R2
    c_dtheta = dtheta*dtheta*R3

    model.cost_expr_ext_cost = c_eC + c_eL + c_theta + c_dT + c_ddelta + c_dtheta

    # nonlinear track constraints
    radius_sq = (xp - xd)*(xp - xd) + (yp - yd)*(yp - yd)
    constraint.expr = vertcat(radius_sq)

    params = types.SimpleNamespace()
    params.xd = xd
    params.yd = yd
    params.grad_xd = grad_xd
    params.grad_yd = grad_yd
    params.phi_d = phi_d
    params.theta_hat = theta_hat
    params.Q1 = Q1
    params.Q2 = Q2
    params.R1 = R1
    params.R2 = R2
    params.R3 = R3
    params.q = q
    params.lr = lr
    params.lf = lf
    params.m = m
    params.I = I
    params.Df = Df
    params.Cf = Cf
    params.Bf = Bf
    params.Dr = Dr
    params.Cr = Cr
    params.Br = Br
    params.Cm1 = Cm1
    params.Cm2 = Cm2
    params.Cd = Cd
    params.Croll = Croll

    model.f_impl_expr = xdot - f_expl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.p = p
    model.name = model_name
    model.params = params

    # return model, constraint
    return model, constraint
