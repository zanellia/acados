import math
import numpy as np
from casadi import SX, types, vertcat


def reduced_model():

    # define structs
    model = types.SimpleNamespace()
    model_name = "reduced_model"

    """ Constraints """
    model.x_min = -2.0
    model.x_max = 2.0

    model.v_min = -2.0
    model.v_max = 2.0

    # CasADi - states
    xp = SX.sym("xp")
    x = vertcat(xp)

    # CasADi - input
    v = SX.sym("v")
    u = vertcat(v)

    # xdot
    xpdot = SX.sym("xpdot")
    xdot = vertcat(xpdot)

    # algebraic variables
    z = vertcat([])

    # dynamics
    f_expl = vertcat(
        v
    )

    # cost
    model.cost_expr_ext_cost = (xp - 1) * (xp - 1)

    model.f_impl_expr = xdot - f_expl
    model.f_expl_expr = f_expl
    model.x = x
    model.xdot = xdot
    model.u = u
    model.z = z
    model.name = model_name

    # return model, constraint
    return model
