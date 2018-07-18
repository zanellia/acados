import argparse
from numpy import *
imoport matlplotlib as plt

# command line arguments parser
parser = argparse.ArgumentParser(description='Process some integers.')
parser.add_argument('-NN',  type=int, help='horizon length')
parser.add_argument('-NM',  type=int, help='number of free masses')
parser.add_argument('-NS',  type=int, help='number of RK stages')
parser.add_argument('-RTI', type=int, help='RTI scheme: 0 - zero-RTI, 1 - stadard-RTI')

args = parser.parse_args()
NN = args.NN
NM = args.NM
NS = args.NS
RTI = args.RTI
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

CPU_time  = data[5:5+NSIM]
iters  = data[5+NSIM:]

max_CPU = max(CPU_time)
avg_CPU = mean(CPU_time)
max_iters = int(max(iters))

print('')
print('max CPU time = {0}, avg CPU time = {1}, max iters = {2}'.format(max_CPU, avg_CPU, max_iters))
print('')

