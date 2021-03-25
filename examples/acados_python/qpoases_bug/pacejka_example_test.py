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

from pacejka_settings import pacejka_settings
import numpy as np
import matplotlib.pyplot as plt
from math import sin, cos, pi
from time import perf_counter as clock


# Important options
load_previous_trajectory = True
# load_previous_trajectory = False 
track_name = "arc"  # "arc", "straight"
# qp_solver = "PARTIAL_CONDENSING_HPIPM"  # "PARTIAL_CONDENSING_HPIPM", "FULL_CONDENSING_QPOASES"
qp_solver = "FULL_CONDENSING_QPOASES"
solver_type = "ZO_SQP"  # "SQP_RTI", "SQP", "ZO_SQP"
# solver_type = "SQP"
print_level = 1
plot_update_time = 0.2
niter = 100
if not load_previous_trajectory:
    initial_guess_type = "straight_line"  # "stationary", "follow_track", "straight_line"
    initial_velocity = [1.0, 0.0]


def straight_track_p(t):
    return np.array([t, 0, 1, 0, t, 0])


def arc_track_p(t): #t = arclength
    return np.array([sin(t), 1-cos(t), cos(t), sin(t), t, t])


# define track
def track_p(t):
    if track_name == "arc":
        return arc_track_p(t)
    elif track_name == "straight":
        return straight_track_p(t)


track_dt = 0.01
track_N = 157
track = np.zeros((track_N, 6))
for ii in range(track_N):
    track[ii, :] = track_p(ii * track_dt)

xtrack = np.array(track[:, 0])
ytrack = np.array(track[:, 1])
xrate = np.array(track[:, 2])
yrate = np.array(track[:, 3])
arcLength = np.array(track[:, 4])
tangentAngle = np.array(track[:, 5])
half_track_width = 0.23

# plot the track and its boundaries
plt.ion()
plt.show()
fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xtrack, ytrack, 'b--')
ax.plot(xtrack + half_track_width * yrate, ytrack - half_track_width * xrate, 'k')
ax.plot(xtrack - half_track_width * yrate, ytrack + half_track_width * xrate, 'k')

# optimization problem starts here...
N = 10
Ts = 1.0/N
model, acados_solver, acados_integrator, nx, nu, npar = pacejka_settings(N, Ts, solver_type, qp_solver, print_level)

# cost params
p_Q1 = 1.0  # contouring cost
p_Q2 = 1000.0  # lag cost
p_R1 = 0.3  # dtorque cost
p_R2 = 0.3  # dsteer cost
p_R3 = 0.3  # darclength cost
p_q = 3.0  # arclength cost
# size params
p_lr = 0.038
p_lf = 0.052
p_m = 0.181
p_I = 0.000505
# lateral force params
p_Df = 0.65
p_Cf = 1.5
p_Bf = 5.2
p_Dr = 1.0
p_Cr = 1.45
p_Br = 8.5
# longitudinal force params
p_Cm1 = 0.98028992
p_Cm2 = 0.01814131
p_Cd = 0.02750696
p_Croll = 0.08518052
p = np.array([0, 0, 1, 0, 0, 0,
              p_Q1, p_Q2, p_R1, p_R2, p_R3, p_q,
              p_lr, p_lf, p_m, p_I,
              p_Df, p_Cf, p_Bf, p_Dr, p_Cr, p_Br,
              p_Cm1, p_Cm2, p_Cd, p_Croll])

# initialize trajectory
if load_previous_trajectory:
    X = np.load('X.npy')
    U = np.load('U.npy')
    acados_solver.set(0, "lbx", X[0, :])
    acados_solver.set(0, "ubx", X[0, :])
    for i in range(N):
        xt = X[i, :]
        ut = U[i, :]
        p[0:6] = track_p(xt[8])
        acados_solver.set(i, "p", p)
        acados_solver.set(i, "x", xt)
        acados_solver.set(i, "u", ut)
else:
    x0 = np.array([0, 0, 0, initial_velocity[0], initial_velocity[1], 0, 0, 0, 0])
    u0 = np.zeros(nu)
    acados_solver.set(0, "lbx", x0)
    acados_solver.set(0, "ubx", x0)
    for i in range(N):
        t = i * Ts
        if initial_guess_type == "follow_track":
            p[0:6] = track_p(t)
            xt = np.array([p[0], p[1], p[5], p[2], p[3], 0, 0, 0, t])
        elif initial_guess_type == "stationary":
            p[0:6] = track_p(0)
            xt = np.array([0, 0, 0, 0, 0, 0, 0, 0, 0])
        elif initial_guess_type == "straight_line":
            p[0:6] = track_p(0)
            xt = np.array([t * x0[3], t * x0[4], 0, x0[3], x0[4], 0, 0, 0, t * x0[3]])
        acados_solver.set(i, "p", p)
        acados_solver.set(i, "x", xt)
        acados_solver.set(i, "u", u0)


# find initial set of linearization points
plot, = ax.plot([], [], 'r')
finalX = np.zeros((N, nx))
finalU = np.zeros((N, nu))
for i in range(niter+1):
    print("iteration", i)
    # if i > 0:
    t_before = clock()
    status = acados_solver.solve()
    t_run = clock() - t_before
    import pdb; pdb.set_trace()
    print('solve time,', t_run)
    if status != 0:
        print("Failed on: ", i)
            #raise RuntimeError(
            #    "acados returned status {0} in closed loop iteration {1}. Exiting.".format(status, i))

    x_horizon = []
    y_horizon = []
    for j in range(N):

        x = acados_solver.get(j, "x")
        u = acados_solver.get(j, "u")

        p[0:6] = track_p(x[8])

        # set params for specific stage
        acados_solver.set(j, "p", p)

        x_horizon.append(x[0])
        y_horizon.append(x[1])

        plot.set_xdata(x_horizon)
        plot.set_ydata(y_horizon)

        if i == niter:
            finalX[j, :] = x
            finalU[j, :] = u

    plt.draw()
    plt.pause(plot_update_time)

np.save('X.npy', finalX)
np.save('U.npy', finalU)

input("Press [enter] to quit.")


