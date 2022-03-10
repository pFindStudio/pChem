from blind_search import blind_search 
from close_identify import new_close_search 
from utils import parameter_file_read, pfind_path_find, delete_file

import os 
import shutil 
from argparse import ArgumentParser 
import time 

def run(): 
    current_path = os.getcwd() 
    print('Welcome to use pChem! (version 1.1)')
    
    parser = ArgumentParser() 
    parser.add_argument("--fasta_path", type=str, default='None', help='path to the fasta file')
    parser.add_argument("--msms_path", type=str, default='None', help='path to the msms file') 
    parser.add_argument("--output_path", type=str, default='None', help='path to the result file') 

    args = parser.parse_args() 
    
    cfg_path = os.path.join(current_path, 'pChem.cfg') 
    parameter_dict = parameter_file_read(cfg_path) 

    '''
    target_path = [os.getcwd()] 
    # print(parameter_dict)
    pfind_path_find(parameter_dict['pfind_install_path'], target_path)
    if len(target_path) == 1: 
        print('please set the correct pfind install path!') 
        return 
    if len(target_path) > 2: 
        print('More than one pFind are installed in current path, please set in detail!')
        return  
    parameter_dict['pfind_install_path'] = target_path[1]
    '''

    
    # 判断参数的正确性
    '''
    if os.path.exists(parameter_dict['pfind_install_path']) == False: 
        print('pfind install path is not exist!') 
        return
    '''
    if os.path.exists(parameter_dict['output_path']) == False: 
        print('output path is not exist!')
        return 
    if os.path.exists(parameter_dict['fasta_path']) == False: 
        print('fasta path is not exist!') 
        return 
    for ms_path in parameter_dict['msms_path']: 
        ms_path = ms_path.split('=')[1].strip()
        if os.path.exists(ms_path) == False: 
            print('msms path is not exist!') 
            return 

    
    pparse_output_path = os.path.join(os.path.join(parameter_dict['output_path'], 'source'), 'pParse')
    if os.path.exists(pparse_output_path): 
        shutil.rmtree(pparse_output_path)
 
    start_time = time.time() 
    blind_search(current_path) 
    blind_time = time.time()  
    if parameter_dict['use_close_search'] == 'True':
        new_close_search(current_path) 
        close_time = time.time() 
        print('blind search cost time (s): ', blind_time - start_time)
        print('restricted search cost time (s): ', close_time - blind_time)
    else:
        print('blind search cost time (s): ', blind_time - start_time)
    
    delete_file(current_path, 'psite.txt')
    
    

if __name__ == "__main__": 
    run()
    os.system("pause")