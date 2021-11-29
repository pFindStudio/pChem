from utils import parameter_pick, parameter_modify, parameter_file_read, pfind_path_find
import os 


def parameter_select(cfg_path): 
    parameter_dict = {} 
    parameter_dict['msmspath'] = []
    with open(cfg_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    for line in lines: 
        if 'msmstype' in line: 
            parameter_dict['msmstype'] = parameter_pick(line) 
        if 'msmspath' in line: 
            parameter_dict['msmspath'].append(parameter_pick(line))
        if 'output_path' in line: 
            parameter_dict['output_path'] = parameter_pick(line) 
    return parameter_dict 


def data_preprocess(cfg_path, current_path): 
    pparse_file_path = os.path.join(current_path, 'bin') 
    pfind_file_path = os.path.join(pparse_file_path, 'pFind')
    pparse_file_path = os.path.join(pparse_file_path, 'pParse')
    pparse_exe_path = os.path.join(pparse_file_path, 'pParse.exe') 

    parameter_dict = parameter_file_read(cfg_path) 
    
    
    pparse_output_path = os.path.join(os.path.join(parameter_dict['output_path'], 'source'), 'pParse')
    if not os.path.exists(pparse_output_path): 
        os.makedirs(pparse_output_path) 
    
    # 对raw文件进行预处理
    for msms_path in parameter_dict['msms_path']: 
        _, t_msms_path = msms_path.split('=')
        t_msms_path = t_msms_path[:-1]
        
        cmd = pparse_exe_path + ' -D ' + t_msms_path + ' -O ' + pparse_output_path 
        # 切换环境 
        os.chdir(pparse_file_path)
        receive = os.system(cmd) 
        print(receive)
    
    pf2_list = os.listdir(pparse_output_path) 
    pf2_name_list = []
    for pf2 in pf2_list: 
        if pf2[-4:] == '.pf2':
            pf2_name_list.append(os.path.join(pparse_output_path, pf2)) 
    
    print(parameter_dict['msms_path']) 
    print(pf2_name_list)
    parameter_dict['msmstype'] = 'PF2'
    # assert len(parameter_dict['msms_path']) == len(pf2_name_list) 
    if len(parameter_dict['msms_path']) != len(pf2_name_list): 
        print('number of pf2 not equal to msms_path, there may be exist bug.')
    pf2_path_list = []
    i = 0 
    for msms_path in parameter_dict['msms_path']:
        prefix, _ = msms_path.split('=') 
        pf2_path_list.append(prefix + '=' + pf2_name_list[i] + '\n') 
        i += 1
    while i < len(pf2_name_list): 
        pf2_path_list.append('msmspath' + str(i+1)+'='+pf2_name_list[i] + '\n')
        i += 1
    parameter_dict['msms_path'] = pf2_path_list 
    # print(parameter_dict) 
    parameter_dict['close_mass_diff_number'] = 100 
    parameter_dict['msms_num'] = len(pf2_name_list)
    parameter_dict['mass_calibration'] = 'True' 

    # target_path = [os.getcwd()]
    # print(parameter_dict)
    #pfind_path_find(parameter_dict['pfind_install_path'], target_path)
    # parameter_dict['pfind_install_path'] = target_path[1] 
    parameter_dict['pfind_install_path'] = pfind_file_path
    print(parameter_dict)
    return parameter_dict 


if __name__ == "__main__": 
    current_path = os.getcwd() 
    cfg_path = os.path.join(current_path, 'pChem.cfg')
    data_preprocess(cfg_path, current_path) 
