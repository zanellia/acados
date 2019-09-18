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

from casadi import *
class acados_dae():
    def __init__(self):
        self.f_impl_expr = None #: CasADi expression for the implicit dynamics :math:`F(\dot{x}, x, u, z) = 0`
        self.f_expl_expr = None #: CasADi expression for the explicit dynamics :math:`\dot{x} = f(x, u)`
        self.x = None           #: CasADi variable describing the state of the system
        self.xdot = None        #: CasADi variable describing the derivative of the state wrt time
        self.u = None           #: CasADi variable describing the input of the system
        self.z = None           #: CasADi variable describing the algebraic variables of the DAE
        self.p = None           #: CasADi variable describing parameters of the DAE
        self.name = None        #: name associated with the function

class acados_constraint():
    def __init__(self):
        self.expr = None #: CasADi expression for the constraint
        self.x = None    #: CasADi variable describing the state of the system
        self.u = None    #: CasADi variable describing the input of the system
        self.z = None    #: CasADi variable describing the algebraic variables of the DAE
        self.nc = None   #: number of constraints
        self.name = None #: name associated with the function

def acados_dae_strip_non_num(acados_constraint):
    out = acados_constraint
    del out['f_impl_expr']
    del out['f_expl_expr']
    del out['x']
    del out['xdot']
    del out['u']
    del out['z']
    del out['p']
    return out

def acados_constraint_strip_non_num(acados_constraint):
    out = acados_constraint
    del out['x']
    del out['u']
    del out['z']
    del out['expr']
    del out['nc']
    return out

