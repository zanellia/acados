import argparse
from numpy import *

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
print('Parsing log file for NN = {0}, NM = {1}, NS = {2} and {3}-RTI scheme'.format(NN, NS, NM, RTI))

# import logged data
if RTI == 0:
    OFFLINE_COND = 1
else:
    OFFLINE_COND = 0
NM_PAR = NM - 1
file_name = 'LOG_FILE_OFFLINE_COND_{}_NN_{}_NUM_FREE_MASSES_{}_IRK_STAGES_NUM_{}.txt'.format(OFFLINE_COND, NN, NM, NS)
data = genfromtxt(file_name, delimiter=',')
import pdb; pdb.set_trace()
NSIM = data[4]
print(data)
