# 限定式搜索
# 运行限定式搜索输出最终鉴定结果

import os 
from utils import parameter_file_read, mass_diff_read, expand_modification_ini, modification_ini_path, close_cfg_write, search_exe_path,\
    delete_file, remove_file
from mass_diff_correction import close_mass_correct, small_delta_filter, mass_diff_diff_filter, mass_static, summary_write, \
    mass_select, explain_dict_generate, new_summary_write, update_identification_efficiency    
from pparse import data_preprocess
from confidence_set import accurate_mass_for_result_file

def new_close_search(current_path): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    close_cfg_path = os.path.join(os.path.join(os.path.join(current_path, 'bin'), 'template'), 'close.cfg') 
    # parameter_dict = parameter_file_read(pchem_cfg_path) 
    parameter_dict = data_preprocess(pchem_cfg_path, current_path) 
    original_output_path = parameter_dict['output_path']
    parameter_dict['output_path'] = os.path.join(original_output_path, 'source')

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
    mass_diff_blind_dict = {} 
    total_name_dict = {}

    with open('mod_list.txt', 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    
    for m in lines:
        if len(m) < 2:
            break 
        name, mass = m.strip().split()
        mass_diff_blind_dict[name] = float(mass) 
        simple_name = name.split('.')[0]
        mass_diff_list.append(simple_name) 
        total_name_dict[simple_name] = name 
    # print(mass_diff_blind_dict, total_name_dict)
    
    # print('mass', mass_diff_list)
    
    close_path = os.path.join(parameter_dict['output_path'], 'close') 
    pfind_summary_path = os.path.join(close_path, 'pFind.summary')
    close_res_path = os.path.join(close_path, 'pFind-Filtered.spectra') 
    mod_static_dict, mod_number_dict = mass_static(close_res_path, current_path, mass_diff_list, parameter_dict['side_position'])
    mass_diff_dict, sim_mod_dict = close_mass_correct(current_path, close_res_path, mass_diff_list, parameter_dict, \
                                    mass_diff_blind_dict, total_name_dict) 
    print('most', mass_diff_dict) 
    accurate_mass_for_result_file(pfind_summary_path, mass_diff_dict)

    explain_dict = explain_dict_generate(current_path)
    close_pfind_path = os.path.join(close_path, 'pFind.summary') 
    mod2pep = mod2pep_generate(close_pfind_path, mass_diff_list)
    # print(mass_diff_dict)
    
    mass_diff_pair_rank = new_summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, explain_dict, sim_mod_dict, pattern='close') 
    # 对雷达图指标进行更新 
    update_identification_efficiency(current_path, parameter_dict) 
    
    os.chdir(current_path)
    
    reporting_result_path = os.path.join(original_output_path, 'reporting_summary')
    close_result_path = os.path.join(parameter_dict['output_path'], 'close')

    
    delete_file(current_path, 'modification-null.ini')
    delete_file(current_path, 'modification-new.ini')
    delete_file(current_path, 'mass_diff_list.txt')
    delete_file(current_path, 'mod_list.txt')
    delete_file(current_path, 'pChem.metric')
    delete_file(current_path, 'pChem.summary')
    delete_file(current_path, 'modification-new.ini') 
    delete_file(current_path, 'heat_map.pdf') 
    remove_file(current_path, 'pChem-close.summary', reporting_result_path)
    remove_file(current_path, 'radar.pdf', reporting_result_path)
    # 出现[]可能是因为有个1000ppm的限制 
    
    # 结果文件合并 
    close_summary_path = os.path.join(reporting_result_path, 'pChem-close.summary') 
    summary_path = os.path.join(reporting_result_path, 'pChem.summary') 
    if os.path.exists(close_summary_path) and os.path.exists(summary_path): 
        summary_file_combine(close_summary_path, summary_path) 
        delete_file(reporting_result_path, 'pChem-close.summary')



def summary_file_combine(close_summary_path, original_path): 
    with open(close_summary_path, 'r', encoding='utf-8') as f: 
        close_lines = f.readlines() 
    with open(original_path, 'r', encoding='utf-8') as f: 
        blind_lines = f.readlines() 
    new_lines = [close_lines[0]]

    info_dict = {}
    for i in range(1, len(close_lines)):
        close_elem = close_lines[i].split('\t') 
        info_dict[close_elem[1]] = [close_elem[2], close_elem[5], close_elem[6]]
    # print(info_dict)
    for i in range(1, len(blind_lines)): 
        blind_elem = blind_lines[i].split('\t') 
        if blind_elem[1] in info_dict.keys():
            info_list = info_dict[blind_elem[1]]
            # print(info_list)
            blind_elem[2] = info_list[0]
            blind_elem[5] = info_list[1]
            blind_elem[6] = info_list[2]
        line = ''
        for token in blind_elem:
            line += token + '\t' 
        new_lines.append(line[:-1]) 

    sort_lines = sorted(new_lines[1:], key=lambda k: int(k.strip().split('\t')[5]), reverse=True) 
    idx = 0 
    new_sort_lines = [new_lines[0]]
    for line in sort_lines:
        new_line = str(idx) + '\t' + line.split('\t', 1)[1] 
        idx += 1
        new_sort_lines.append(new_line)

    with open(original_path, 'w', encoding='utf-8') as f:
        for line in new_sort_lines: 
            f.write(line)




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
    name2mass, mass_diff_list = small_delta_filter(mass_diff_list, parameter_dict)
    mod_static_dict, mod_number_dict = mass_static(blind_path, current_path, mass_diff_list) 

    # 将统计结果写入结果文件 
    # summary_write(current_path, mod_static_dict, mod_number_dict) 

    
    if 'mass_diff_diff' in parameter_dict.keys() and parameter_dict['mass_diff_diff'] != -1.0:
        mass_diff_list = mass_diff_diff_filter(name2mass, mass_diff_list, parameter_dict['mass_diff_diff'])
    
    mass_diff_list = mass_select(mass_diff_list, parameter_dict['close_mass_diff_number'], name2mass)
    
    
    # 对未知质量数质量做校正
    mass_diff_dict = close_mass_correct(current_path, blind_path, mass_diff_list)
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
