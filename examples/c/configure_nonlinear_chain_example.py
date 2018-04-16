#!/usr/bin/env python3
import os

OFFLINE_COND = ['1','0']
NN = ['10', '20', '30']
NUM_FREE_MASSES = ['2', '3', '4', '5']
IRK_STAGES_NUM = ['2','3', '4', '6', '8']

for i in range(len(OFFLINE_COND)):
    for j in range(len(NN)):
        for k in range(len(NUM_FREE_MASSES)):
            for l in range(len(IRK_STAGES_NUM)):
                print('OFFLINE_COND = {0}, NN = {1}  NUM_FREE_MASSES = {2}, IRK_STAGES_NUM = {3}'.format(OFFLINE_COND[i], NN[j], NUM_FREE_MASSES[k], IRK_STAGES_NUM[l]))
                source_code_file  = open('nonlinear_chain_ocp_nlp_no_interface_nmpc_code_generation.tmp', 'r')
                source_code = source_code_file.read()
                
                source_code = source_code.replace('@OFFLINE_COND@', OFFLINE_COND[i])
                source_code = source_code.replace('@NN@', NN[j])
                source_code = source_code.replace('@NUM_FREE_MASSES@', NUM_FREE_MASSES[k])
                source_code = source_code.replace('@IRK_STAGES_NUM@',  IRK_STAGES_NUM[l])

                configured_source_code_file = open('nonlinear_chain_ocp_nlp_no_interface_nmpc_code_generation.c','w') 
                configured_source_code_file.write(source_code)
                configured_source_code_file.close()

                os.chdir('../../build')
                curr_path = os.getcwd()
                print('In %s' %curr_path)
                os.system('make -j4 >/dev/null 2>&1')
                os.chdir('examples/c')
                curr_path = os.getcwd()
                print('In %s' %curr_path)
                os.system('./nonlinear_chain_ocp_nlp_no_interface_nmpc_code_generation_example')
                os.chdir('../../../examples/c')
                curr_path = os.getcwd()
                print('In %s' %curr_path)
                
