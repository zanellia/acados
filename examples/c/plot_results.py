from numpy import *
import matplotlib.pyplot as plt

SCALE = 'LIN'
NN_v = [10, 20, 30, 40, 50, 60, 80, 100, 120, 140, 160, 180, 200, 220, 240, 260, 280]
NN_v = [10, 20, 30, 40, 50, 60, 80, 100]
NN_v = [10, 20, 30]
data_points = len(NN_v)
NM = 7
NS = 4
INIT_DISCARD = 0 # number of initial points to be discarded

zero_rti = zeros((data_points, 1))
std_rti = zeros((data_points, 1))
for i in range(data_points):
    for j in range(2):
        print(i)
        NN = NN_v[i]
        RTI = j 
        print('')
        print('---------------------------------------------------------------')
        print('Parsing log file for NN = {0}, NM = {1}, NS = {2} and {3}-RTI scheme'.format(NN, NS, NM, RTI))
        print('---------------------------------------------------------------')

        # import logged data
        if RTI == 0:
            OFFLINE_COND = 1
        else:
            OFFLINE_COND = 0
        NM_PAR = NM - 1
        file_name = 'LOG_FILE_OFFLINE_COND_{}_NN_{}_NUM_FREE_MASSES_{}_IRK_STAGES_NUM_{}.txt'.format(OFFLINE_COND, NN, NM_PAR, NS)

        data = genfromtxt(file_name, delimiter=',')

        NSIM = int(data[4])

        CPU_time  = data[5+INIT_DISCARD:5+NSIM]
        iters  = data[5+NSIM+INIT_DISCARD:]

        max_CPU = max(CPU_time)
        avg_CPU = mean(CPU_time)
        min_CPU = min(CPU_time)
        max_iters = int(max(iters))

        print('')
        print('max CPU time = {0}, avg CPU time = {1}, max iters = {2}'.format(max_CPU, avg_CPU, max_iters))
        print('')

        # import pdb; pdb.set_trace()
        if j == 0:
            # zero_rti[i] = min_CPU 
            # zero_rti[i] = avg_CPU 
            zero_rti[i] = max_CPU 
        else:
            # std_rti[i] = min_CPU
            # std_rti[i] = avg_CPU
            std_rti[i] = max_CPU

# fit cubic
alpha = polyfit(NN_v, std_rti, 3)

N_grid = linspace(0, NN_v[data_points-1], 50) 
cubic_fit = alpha[3]*N_grid + alpha[2]*N_grid + alpha[1]*N_grid**2 + alpha[0]*N_grid**3

plt.figure()
if SCALE == 'LOGLOG':
    plt.loglog(NN_v, zero_rti)
    plt.loglog(NN_v, std_rti)
    # plt.loglog(NN_v, cubic_fit)

if SCALE == 'SEMILOG':
    plt.semilogy(NN_v, zero_rti)
    plt.semilogy(NN_v, std_rti)
    # plt.semilogy(N_grid, cubic_fit)

if SCALE == 'LIN':
    plt.plot(NN_v, zero_rti)
    plt.plot(NN_v, std_rti)
    # plt.plot(N_grid, cubic_fit, '--')

plt.ylabel('CPU time [ms]')
plt.xlabel('N')
plt.legend(['0-RTI','std RTI','cubic fit'])
plt.grid()
plt.show()
