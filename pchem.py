from blind_search import blind_search 
from close_identify import new_close_search 
from utils import parameter_file_read 

import os 
import shutil 
from argparse import ArgumentParser 

if __name__ == "__main__":
    current_path = os.getcwd() 
    print('Welcome to use pChem!')
    
    parser = ArgumentParser() 
    parser.add_argument("--fasta_path", type=str, default='None', help='path to the fasta file')
    parser.add_argument("--msms_path", type=str, default='None', help='path to the msms file') 
    parser.add_argument("--output_path", type=str, default='None', help='path to the result file') 

    args = parser.parse_args() 
   
    cfg_path = os.path.join(current_path, 'pChem.cfg') 
    parameter_dict = parameter_file_read(cfg_path)
    
    pparse_output_path = os.path.join(parameter_dict['output_path'], 'pParse') 
    if os.path.exists(pparse_output_path): 
        shutil.rmtree(pparse_output_path)
 
    
    blind_search(current_path) 
    if parameter_dict['use_close_search'] == 'True':
        new_close_search(current_path) 
    
    os.system("pause")
    