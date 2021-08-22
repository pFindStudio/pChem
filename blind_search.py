# 在人工筛选常见修饰后，对数据进行盲搜鉴定

import os 
from utils import parameter_file_read, modification_ini_path, modification_ini_dict, \
    modification_ini_generation, blind_cfg_write, search_exe_path, mass_diff_list_generate, \
    modification_ini_generation_from_param, mass_diff_read
from mass_diff_correction import mass_correct, small_delta_filter, mass_diff_diff_filter, \
    mass_static, summary_write, mass_select, new_summary_write, unimod_dict_generate   


def blind_search(current_path):
    
    # 路径参数 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    blind_cfg_path = os.path.join(os.path.join(current_path, 'template'), 'blind.cfg')

    parameter_dict = parameter_file_read(pchem_cfg_path) 

    
    # 读取所有的modification作为查找字典 
    modification_path = modification_ini_path(parameter_dict)
    modification_dict = modification_ini_dict(modification_path)
    
    unimod_dict = unimod_dict_generate(modification_dict) 

    # 重新生成modification.ini文件
    if parameter_dict['open_flag'] == 'True':
        # 读取common_modification_list.txt 文件来生成
        common_modification_list = modification_ini_generation(current_path, modification_dict)
    else: 
        # 使用参数文件中的列表来生成 
        common_modification_list = modification_ini_generation_from_param(current_path, modification_dict, parameter_dict)
    #common_modification_list = common_modification_list[:2]
    #print(common_modification_list)
    
    # 重新生成blind.cfg文件
    res_path = blind_cfg_write(blind_cfg_path, current_path, parameter_dict, common_modification_list)
    
    
    # 运行盲搜search.exe进行搜索
    bin_path, exe_path = search_exe_path(parameter_dict)
    cmd = exe_path + ' ' + blind_cfg_path 
    os.chdir(bin_path)
    receive = os.system(cmd) 
    print(receive)
    

    
    # 读取鉴定结果，生成位置修饰的候选列表
    mass_diff_list_generate(res_path, current_path)
    mass_diff_list, mod2pep = mass_diff_read(current_path) 

    
    # 盲搜的结果文件
    
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    #print(mass_diff_list)
    
    # 对得到的未知质量数进行过滤和统计 

    # mass_diff_list 修饰名称的列表  mod_static_dict 修饰名字-> Counter
    # mod_number_dict 修饰名字 -> PSM 
    # mass_diff_dict 修饰名字 -> 精准质量 
    # mod2pep 修饰名字 -> peptide 

    # 过滤小于200和复数修饰 
    name2mass, mass_diff_list = small_delta_filter(mass_diff_list, parameter_dict['min_mass_modification']) 
    mod_static_dict, mod_number_dict = mass_static(blind_path, current_path, mass_diff_list) 

    # 计算精确质量 
    #system_correct={mean, median}, mod_correct={mean, median, weight}
    mass_diff_dict = mass_correct(current_path, blind_path, mass_diff_list, system_correct='mean', mod_correct='mean') 
    #mass_diff_dict = mass_correct(current_path, blind_path, mass_diff_list, system_correct='median', mod_correct='median') 

    # 将统计结果写入结果文件 
    new_summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict) 
    


if __name__ == "__main__": 
    current_path = os.getcwd() 
    blind_search(current_path)
