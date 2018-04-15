#!/usr/bin/env python3

OFFLINE_COND = ['0','1']
NN = ['10', '20', '30', '40', '60', '80', '100']
NUM_FREE_MASSES = ['2', '3', '4', '5']
IRK_STAGES_NUM = ['2','3', '4', '6', '8']

for i in range(len(NN_v)):
    for j in range(len(NN_v)):
        for k in range(len(NN_v)):
            for l in range(len(NN_v)):
                source_code_file  = open('nonlinear_chain_ocp_nlp_no_interface_nmpc_code_generation.tmp', 'r')
                source_code = source_code_file.read()
                import pdb; pdb.set_trace()
                
                source_code = source_code.replace('@OFFLINE_COND@', OFFLINE_COND[i])
                source_code = source_code.replace('@NN@', NN[i])
                source_code = source_code.replace('@NUM_FREE_MASSES@', NUM_FREE_MASSES[i])
                source_code = source_code.replace('@IRK_STAGES_NUM@',  IRK_STAGES_NUM[i])

                configured_source_code_file = open('nonlinear_chain_ocp_nlp_no_interface_nmpc_code_generation.c','w') 
                configured_source_code_file.write(source_code)
                configured_source_code_file.close()
