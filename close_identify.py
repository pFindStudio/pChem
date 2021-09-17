# 限定式搜索
# 运行限定式搜索输出最终鉴定结果

import os 
from utils import parameter_file_read, mass_diff_read, expand_modification_ini, modification_ini_path, close_cfg_write, search_exe_path,\
    delete_file, remove_file
from mass_diff_correction import mass_correct, small_delta_filter, mass_diff_diff_filter, mass_static, summary_write, \
    mass_select, explain_dict_generate, new_summary_write, update_identification_efficiency    
from pparse import data_preprocess


def new_close_search(current_path): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    close_cfg_path = os.path.join(os.path.join(os.path.join(current_path, 'bin'), 'template'), 'close.cfg') 
    # parameter_dict = parameter_file_read(pchem_cfg_path) 
    parameter_dict = data_preprocess(pchem_cfg_path, current_path)

    close_output_path = os.path.join(parameter_dict['output_path'], 'close')
    if not os.path.exists(close_output_path): 
        os.makedirs(close_output_path) 

    
    # 进行限定式搜索 
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + close_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd)
    print(receive) 
    
    
    os.chdir(current_path)
    # 对限定式结果进行分析 
    mass_diff_list = [] 
    with open('mod_list.txt', 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    
    for m in lines:
        if len(m) < 2:
            break 
        mass_diff_list.append(m.strip().split('.')[0])
    
    print('mass', mass_diff_list)
    
    close_path = os.path.join(parameter_dict['output_path'], 'close') 
    close_res_path = os.path.join(close_path, 'pFind-Filtered.spectra') 
    mod_static_dict, mod_number_dict = mass_static(close_res_path, current_path, mass_diff_list, parameter_dict['side_position'])
    mass_diff_dict = mass_correct(current_path, close_res_path, mass_diff_list, system_correct='mean', mod_correct='mean') 
    explain_dict = explain_dict_generate(current_path)
    close_pfind_path = os.path.join(close_path, 'pFind.summary') 
    mod2pep = mod2pep_generate(close_pfind_path, mass_diff_list)
    print(mass_diff_dict)
    
    mass_diff_pair_rank = new_summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, explain_dict, pattern='close') 
    # 对雷达图指标进行更新 
    update_identification_efficiency(current_path, parameter_dict) 
    
    os.chdir(current_path)
    
    reporting_result_path = os.path.join(parameter_dict['output_path'], 'reporting_summary')
    close_result_path = os.path.join(parameter_dict['output_path'], 'close')

    delete_file(current_path, 'modification-null.ini')
    delete_file(current_path, 'modification-new.ini')
    delete_file(current_path, 'mass_diff_list.txt')
    delete_file(current_path, 'mod_list.txt')
    delete_file(current_path, 'pChem.metric')
    delete_file(current_path, 'pChem.summary')
    delete_file(current_path, 'modification-new.ini') 
    delete_file(current_path, 'heat_map.pdf') 
    remove_file(current_path, 'pChem-close.summary', close_result_path)
    remove_file(current_path, 'radar.pdf', reporting_result_path)





def mod2pep_generate(close_pfind_path, mass_diff_list): 
    mod2pep = {}

    with open(close_pfind_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 

    for m in mass_diff_list:
        for line in lines: 
            if m in line: 
                mod2pep[m] = int(line.split('\t')[1].split()[0])
    return mod2pep

def close_search(current_path):
    
    
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    close_cfg_path = os.path.join(os.path.join(os.path.join(current_path, 'bin'), 'template'), 'close.cfg')

    # 读取位置修饰质量数的列表
    mass_diff_list, _ = mass_diff_read(current_path)

    
    # 盲搜的结果文件
    parameter_dict = parameter_file_read(pchem_cfg_path)
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    # print(mass_diff_list)
    
    # 对得到的未知质量数进行过滤和统计 
    name2mass, mass_diff_list = small_delta_filter(mass_diff_list, parameter_dict['min_mass_modification'])
    mod_static_dict, mod_number_dict = mass_static(blind_path, current_path, mass_diff_list) 

    # 将统计结果写入结果文件 
    # summary_write(current_path, mod_static_dict, mod_number_dict) 

    
    if 'mass_diff_diff' in parameter_dict.keys() and parameter_dict['mass_diff_diff'] != -1.0:
        mass_diff_list = mass_diff_diff_filter(name2mass, mass_diff_list, parameter_dict['mass_diff_diff'])
    
    mass_diff_list = mass_select(mass_diff_list, parameter_dict['close_mass_diff_number'], name2mass)
    
    
    # 对未知质量数质量做校正
    mass_diff_dict = mass_correct(current_path, blind_path, mass_diff_list)
    print(mass_diff_dict)

    
    # 统计未知质量发生的位点 
    # mod_static_dict = mass_static(blind_path, current_path, mass_diff_list) 
    # print(mod_static_dict)
    # print(mod_static_dict['PFIND_DELTA_252.12'].most_common()[0][0])


    # 修改的modification-new.ini文件，加入质量数修饰
    #ini_path = modification_ini_path(parameter_dict)
    #new_ini_path = expand_modification_ini(mass_diff_dict, mod_static_dict, current_path, ini_path)
    

    
    # 生成限定式参数文件 
    res_path = close_cfg_write(close_cfg_path, current_path, parameter_dict, mass_diff_dict)


    # 调用pfind进行限定式搜索 
    
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + close_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd)
    print(receive)
    

if __name__ == "__main__":
    current_path = os.getcwd() 
    new_close_search(current_path)
