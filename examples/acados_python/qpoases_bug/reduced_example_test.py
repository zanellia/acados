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
import os
sys.path.insert(0, 'common')

from reduced_settings import reduced_settings
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from time import perf_counter as clock


# Important options
soft_constraints = True
N = 12  # NOTE: N=12 breaks QPOASES with soft constraints, N=11 works
qp_solver = "FULL_CONDENSING_QPOASES"  # "PARTIAL_CONDENSING_HPIPM", "FULL_CONDENSING_QPOASES"
solver_type = "SQP_RTI"  # "SQP_RTI", "SQP", "ZO_SQP"
print_level = 0
plot_update_time = 0.5
niter = 50

load_previous_trajectory = False
initial_velocity = 1.0

track_dt = 0.01
track_N = 157
track = np.zeros(track_N)
for ii in range(track_N):
    track[ii] = ii * track_dt

xtrack = np.array(track)

# plot the track and its boundaries
plt.ion()
plt.show()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(0.0 * np.ones(N), '--k')
ax.plot(1.0 * np.ones(N), '--k')

# optimization problem starts here...
Ts = 1.0/(N-1)
model, acados_solver, acados_integrator, nx, nu, npar = reduced_settings(N, Ts, solver_type, qp_solver, soft_constraints, print_level)

# initialize trajectory
if load_previous_trajectory:
    X = np.load('Xconv.npy')
    U = np.load('Uconv.npy')
    acados_solver.set(0, "lbx", X[0, :])
    acados_solver.set(0, "ubx", X[0, :])
    for i in range(N):
        xt = X[i, :]
        ut = U[i, :]
        acados_solver.set(i, "x", xt)
        acados_solver.set(i, "u", ut)
else:
    x0 = np.array([0])
    u0 = np.zeros(nu)
    acados_solver.set(0, "lbx", x0)
    acados_solver.set(0, "ubx", x0)
    for i in range(N):
        t = i * Ts
        xt = np.array([t * initial_velocity])
        acados_solver.set(i, "x", xt)
        acados_solver.set(i, "u", u0)


# find initial set of linearization points
plotx0, = ax.plot([ii for ii in range(N)], np.zeros(N), '-bo')
finalX = np.zeros((N, nx))
finalU = np.zeros((N, nu))
for i in range(niter+1):
    print(i)
    if i > 0:
        t_before = clock()
        status = acados_solver.solve()
        t_run = clock() - t_before
        print('solve time,', t_run)
        if status != 0:
            print("Failed on: ", i)
            #raise RuntimeError(
            #    "acados returned status {0} in closed loop iteration {1}. Exiting.".format(status, i))

    x_horizon = []
    for j in range(N):

        x = acados_solver.get(j, "x")
        u = acados_solver.get(j, "u")

        x_horizon.append(x[0])

        if i == niter:
            finalX[j, :] = x
            finalU[j, :] = u

        if j == N-1:
            plotx0.set_ydata(x_horizon)

    plt.draw()
    plt.pause(plot_update_time)

np.save('Xconv.npy', finalX)
np.save('Uconv.npy', finalU)

input("Press [enter] to quit.")


