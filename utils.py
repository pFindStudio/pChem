import os
import shutil  


# 读取参数行内容并返回
def parameter_pick(line):
    eq_idx = line.find('=')
    parameter_content = line[eq_idx+1:].strip()
    return parameter_content


# 改写参数行的内容并返回
def parameter_modify(line, content):
    eq_idx = line.find('=')
    line = line[:eq_idx+1]
    line += content
    line += '\n'
    return line


# 读取参数文件
def parameter_file_read(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    parameter_dict = {} 
    parameter_dict['psite_run'] = 'True'

    for i in range(len(lines)):
        if 'pfind_install' in lines[i]:
            parameter_dict['pfind_install_path'] = parameter_pick(lines[i])
        if 'fasta_path' in lines[i]:
            parameter_dict['fasta_path'] = parameter_pick(lines[i])
        if 'msmsnum' in lines[i]:
            msmsnum = int(parameter_pick(lines[i]))
            parameter_dict['msms_num'] = msmsnum 
            parameter_dict['msms_path'] = []
            i += 1
            for _ in range(msmsnum):
                parameter_dict['msms_path'].append(lines[i])
                i += 1
            # parameter_dict['msms_path'].reverse()
        if 'output_path' in lines[i]:
            parameter_dict['output_path'] = parameter_pick(lines[i])
            if not os.path.exists(parameter_dict['output_path']):
                os.mkdir(parameter_dict['output_path'])
        if 'open_flag' in lines[i]:
            parameter_dict['open_flag'] = parameter_pick(lines[i])
        if 'common_modification_list' in lines[i]:
            parameter_dict['common_modification_list'] = parameter_pick(lines[i]) 
        if 'mass_calibration' in lines[i]:
            parameter_dict['mass_calibration'] = parameter_pick(lines[i]) 
        if 'mass_of_diff_diff' in lines[i]:
            parameter_dict['mass_of_diff_diff'] = float(parameter_pick(lines[i])) 
        if 'common_modification_number' in lines[i]:
            parameter_dict['common_modification_number'] = int(parameter_pick(lines[i]))
        if 'close_mass_diff_number' in lines[i]:
            parameter_dict['close_mass_diff_number'] = int(parameter_pick(lines[i])) 
        if 'min_mass_modification' in lines[i]:
            parameter_dict['min_mass_modification'] = float(parameter_pick(lines[i]))
        if 'max_mass_modification' in lines[i]:
            parameter_dict['max_mass_modification'] = float(parameter_pick(lines[i]))
        if 'mass_diff_diff_range' in lines[i]:
            parameter_dict['mass_diff_diff_range'] = float(parameter_pick(lines[i])) 
        if 'filter_frequency' in lines[i]: 
            parameter_dict['filter_frequency'] = float(parameter_pick(lines[i])) 
        if 'side_position' in lines[i]: 
            parameter_dict['side_position'] = parameter_pick(lines[i]) 
        if 'activation_type' in lines[i]:
            parameter_dict['activation_type'] = parameter_pick(lines[i]) 
        #if 'use_close_search' in lines[i]: 
        #    parameter_dict['use_close_search'] = parameter_pick(lines[i]) 
        if 'msmstype' in lines[i]: 
            parameter_dict['msmstype'] = parameter_pick(lines[i]) 
        if 'report_statistics' in lines[i]: 
            parameter_dict['report_statistical'] = parameter_pick(lines[i]) 
        if 'isotope_labeling' in lines[i]: 
            parameter_dict['isotope_labeling'] = parameter_pick(lines[i]) 
        if 'p_value_threshold' in lines[i]:
            parameter_dict['p_value_threshold'] = float(parameter_pick(lines[i]))
        parameter_dict['close_mass_diff_number'] = 10 
    if parameter_dict['isotope_labeling'] == 'True':
        parameter_dict['use_close_search'] = 'True'
    else:
        parameter_dict['use_close_search'] = 'False'
    return parameter_dict


# 写open参数文件
def open_cfg_write(cfg_path, parameter_dict):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 指定参数内容修改
    new_lines = []
    for i in range(len(lines)): 
        if 'msmstype' in lines[i]: 
            lines[i] = parameter_modify(lines[i], parameter_dict['msmstype'])
        if 'activation_type' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['activation_type'])
        if 'modpath' in lines[i]:
            mod_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
            mod_path = os.path.join(mod_path, 'modification.ini')
            lines[i] = parameter_modify(lines[i], mod_path)
        if 'fastapath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['fasta_path'])
        if 'outputpath' in lines[i]:
            res_path = os.path.join(parameter_dict['output_path'], 'open')
            if not os.path.exists(res_path):
                os.mkdir(res_path)
            lines[i] = parameter_modify(lines[i], res_path)
        if 'outputname' in lines[i]:
            lines[i] = parameter_modify(lines[i], 'open') 
        if 'maxdelta' in lines[i]: 
            lines[i] = parameter_modify(lines[i], str(parameter_dict['max_mass_modification']))
        if 'msmsnum' in lines[i]:
            lines[i] = parameter_modify(lines[i], str(parameter_dict['msms_num']))
            new_lines.append(lines[i])
            for path in parameter_dict['msms_path']:
                new_lines.append(path)
            continue
        # 之前的地址需要清除 
        if 'msmspath' in lines[i]:
            continue
        new_lines.append(lines[i])
    
    # 写入参数文件
    with open(cfg_path, 'w', encoding='utf-8') as f:
        for line in new_lines:
            f.write(line)
    return res_path


# 写blind参数文件
def blind_cfg_write(cfg_path, current_path, parameter_dict, common_modification_list):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    mod_line = ""
    for mod in common_modification_list:
        mod_line += (mod + ';')
    mod_line = mod_line[:-1]

    # 指定参数内容修改
    new_lines = []
    for i in range(len(lines)): 
        if 'msmstype' in lines[i]: 
            lines[i] = parameter_modify(lines[i], parameter_dict['msmstype'])
        if 'activation_type' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['activation_type'])
        if 'selectmod' in lines[i]:
            lines[i] = parameter_modify(lines[i], mod_line)
        if 'modpath' in lines[i]:
            mod_path = os.path.join(current_path, 'modification-null.ini')
            lines[i] = parameter_modify(lines[i], mod_path)
        if 'fastapath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['fasta_path'])
        if 'outputpath' in lines[i]:
            res_path = os.path.join(parameter_dict['output_path'], 'blind')
            if not os.path.exists(res_path):
                os.mkdir(res_path)
            lines[i] = parameter_modify(lines[i], res_path)
        if 'outputname' in lines[i]:
            lines[i] = parameter_modify(lines[i], 'blind') 
        if 'maxdelta' in lines[i]: 
            lines[i] = parameter_modify(lines[i], str(parameter_dict['max_mass_modification']))
        if 'msmsnum' in lines[i]:
            lines[i] = parameter_modify(lines[i], str(parameter_dict['msms_num']))
            new_lines.append(lines[i])
            for path in parameter_dict['msms_path']:
                new_lines.append(path)
            continue
        if 'msmspath' in lines[i]:
            continue

        new_lines.append(lines[i])
    
    # 写入参数文件
    with open(cfg_path, 'w', encoding='utf-8') as f:
        for line in new_lines:
            f.write(line)
    return res_path


# 写close参数文件
def close_cfg_write(cfg_path, current_path, parameter_dict, mass_diff_pair_rank):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    mod_line, flag = combine_common_list(current_path, parameter_dict, mass_diff_pair_rank)
    # print(mod_line)
    # 指定参数内容修改 
    new_lines = [] 
    for i in range(len(lines)): 
        if 'msmstype' in lines[i]: 
            lines[i] = parameter_modify(lines[i], parameter_dict['msmstype'])
        if 'activation_type' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['activation_type'])
        if 'selectmod' in lines[i]:
            lines[i] = parameter_modify(lines[i], mod_line)
        #if 'fixmod' in lines[i] and flag == True:
        #    lines[i] = parameter_modify(lines[i],'Carbamidomethyl[C]')
        if 'modpath' in lines[i]:
            mod_path = os.path.join(current_path, 'modification-new.ini')
            lines[i] = parameter_modify(lines[i], mod_path)
        if 'fastapath' in lines[i]:
            lines[i] = parameter_modify(lines[i], parameter_dict['fasta_path'])
        if 'outputpath' in lines[i]:
            res_path = os.path.join(parameter_dict['output_path'], 'close')
            if not os.path.exists(res_path):
                os.mkdir(res_path)
            lines[i] = parameter_modify(lines[i], res_path)
        if 'outputname' in lines[i]:
            lines[i] = parameter_modify(lines[i], 'close')
        if 'maxdelta' in lines[i]: 
            lines[i] = parameter_modify(lines[i], str(parameter_dict['max_mass_modification']))
        if 'msmsnum' in lines[i]:
            lines[i] = parameter_modify(lines[i], str(parameter_dict['msms_num']))
            new_lines.append(lines[i])
            for path in parameter_dict['msms_path']:
                new_lines.append(path)
            continue
        if 'msmspath' in lines[i]:
            continue

        new_lines.append(lines[i])
    
    # 写入参数文件
    with open(cfg_path, 'w', encoding='utf-8') as f:
        for line in new_lines:
            f.write(line)




# 盲搜之后生成合并的修饰列表
def combine_common_list(current_path, parameter_dict, mass_diff_pair_rank):
    mod_line = ""
    '''
    if parameter_dict['open_flag'] == 'True': 
        common_path = os.path.join(current_path, 'common_modification_list.txt')
        with open(common_path, 'r', encoding='utf-8') as f:
            lines = f.readlines()
        lines = lines[1:]
        for line in lines:
            if len(line) < 4:
                break
            mod = line.split('\t')[0]
            mod_line += mod + ';' 
    else:
        mod_line += parameter_dict['common_modification_list'] 
        if mod_line[-1] != ';':
            mod_line += ';'
    '''

    for key in mass_diff_pair_rank: 
        if parameter_dict['isotope_labeling'] == 'False':
            short_key = 'PFIND_DELTA_' + key.split('.')[0]
        else:
            short_key = key.split('.')[0]
        mod_line += short_key + ';'
    
    flag = False 
    '''
    if 'Carbamidomethyl[C]' in mod_line:
        mod_list = mod_line.split(';')
        mod_line = ''
        for mod in mod_list:
            if mod == '':
                continue
            if mod == 'Carbamidomethyl[C]':
                flag = True 
                continue
            mod_line += mod + ';'
    '''

    return mod_line, flag


# 返回search.exe地址
def search_exe_path(parameter_dict):
    bin_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
    exe_path = os.path.join(bin_path, 'Searcher.exe')
    return bin_path, exe_path 


# 返回modification.ini地址
def modification_ini_path(parameter_dict):
    bin_path = os.path.join(parameter_dict['pfind_install_path'], 'bin')
    modification_ini_path = os.path.join(bin_path, 'modification.ini')
    return modification_ini_path


# 读取谱图结果文件, 将常见的修饰列表写出来提供删选
def spectra_result_read(spectra_res_path, target_path, common_modification_number):
    # 保存常见的修饰
    model_res = []
    with open(spectra_res_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    for i in range(len(lines)):
        if 'Modifications:' in lines[i]:
            i += 1
            model_res.append(lines[i])
            i += 1
            '''
            while True:
                freq = lines[i].split('\t')[1]
                idx_left = freq.find('(')
                idx_right = freq.find('%')
                freq = float(freq[idx_left+1:idx_right])
                if freq < 1.5:
                    break
                model_res.append(lines[i])
                i += 1
            '''
            num = 0 
            while num < common_modification_number:
                model_res.append(lines[i])
                i += 1
                num += 1
    # print(model_res)
    target_path = os.path.join(target_path, 'common_modification_list.txt')
    with open(target_path, 'w', encoding='utf-8') as f:
        for line in model_res:
            f.write(line)


# 读取所有候选修饰，并返回dict
def modification_ini_dict(path):
    with open(path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    modification_dict = {}
    i = 1
    while i < len(lines):
        # 防止文件后面的空行
        if len(lines[i]) < 4:
            break
        modification_name = lines[i].split()[0]
        eq_idx = modification_name.find('=')
        modification_name = modification_name[eq_idx+1:]
        modification_dict[modification_name] = lines[i] + lines[i+1]
        i += 2
    return modification_dict


# 生成modification-null.ini用于检索
def modification_ini_generation(path, modification_dict):
    # 读取之前确定的常见修饰列表
    common_list_path = os.path.join(path, 'common_modification_list.txt')
    with open(common_list_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    lines = lines[1:]
    common_modification_list = []
    modification_ini_lines = []
    for line in lines:
        if len(line) < 4:
            break
        mod_name = line.split('\t')[0]
        common_modification_list.append(mod_name)
        modification_ini_lines.append(modification_dict[mod_name])
    # print(common_modification_list)
    # print(modification_ini_lines)
    
    # 写入新的modification-null.ini
    new_ini_path = os.path.join(path, 'modification-null.ini')
    i = 1
    with open(new_ini_path, 'w', encoding='utf-8') as f:
        f.write('@NUMBER_MODIFICATION=' + str(len(common_modification_list)) + '\n')
        for line in modification_ini_lines:
            if 'name' in line:
                eq_idx = line.find('=')
                new_line = 'name' + str(i) + line[eq_idx:]
                i += 1
                f.write(new_line)
            else:
                f.write(line)
    return common_modification_list


# 从参数文件中生成modification-null.ini用于检索
def modification_ini_generation_from_param(path, modification_dict, parameter_dict):
    modification_list = parameter_dict['common_modification_list'].split(';')
    common_modification_list = []
    modification_ini_lines = []
    for modification in modification_list: 
        if modification in modification_dict.keys():
            common_modification_list.append(modification)
            modification_ini_lines.append(modification_dict[modification])
    
    # 写入新的modification-null.ini
    new_ini_path = os.path.join(path, 'modification-null.ini')
    i = 1
    with open(new_ini_path, 'w', encoding='utf-8') as f:
        f.write('@NUMBER_MODIFICATION=' + str(len(common_modification_list)) + '\n')
        for line in modification_ini_lines:
            if 'name' in line:
                eq_idx = line.find('=')
                new_line = 'name' + str(i) + line[eq_idx:]
                i += 1
                f.write(new_line)
            else:
                f.write(line)
    return common_modification_list


# 读取盲搜结果并生成未知质量数修饰
def mass_diff_list_generate(res_path, current_path):
    summary_path = os.path.join(res_path, 'pFind.summary')
    with open(summary_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    mass_diff_lines = []
    i = 0
    while i < len(lines):
        if 'Modifications:' in lines[i]:
            i += 1
            mass_diff_lines.append(lines[i])
            i += 1
            while i < len(lines):
                if '------' in lines[i]:
                    break 
                if 'PFIND' in lines[i]:
                    mass_diff_lines.append(lines[i])
                i += 1
        i += 1
    # 写入txt文件方便筛选
    mass_diff_path = os.path.join(current_path, 'mass_diff_list.txt')
    with open(mass_diff_path, 'w', encoding='utf-8') as f:
        for line in mass_diff_lines:
            f.write(line)


# 读取mass_diff并返回列表
def mass_diff_read(path):
    mass_diff_txt = os.path.join(path, 'mass_diff_list.txt')
    with open(mass_diff_txt, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    lines = lines[1:]
    mass_diff_list = []
    mod2pep = {} 
    for line in lines: 
        # print(line)
        if len(line) < 2:
            break 
        mod, pep = line.split('\t')[0], line.split('\t')[1].split()[0] 
        mass_diff_list.append(mod)
        mod2pep[mod] = pep 
    return mass_diff_list, mod2pep 



# 在盲搜确定修饰质量后，将其加入modification-new.ini文件
def expand_modification_ini(mass_diff_pair_rank, mass_diff_dict, mod_static_dict, current_path, ini_path, parameter_dict, refine_ion_list=None, exist_ion_flag_list=None):
    # new_mass_list = new_mass_list_generate(mass_diff_pair_rank, mass_diff_list)
    with open(ini_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    while True:
        if len(lines[-1]) < 4:
            lines = lines[:-1]
        else:
            break
    # 修改第一行总数
    lines[0], num = lines[0].split('=')
    total_num = int(num) + len(mass_diff_pair_rank)
    lines[0] = lines[0] + '=' + str(total_num) + '\n' 
    lines[-1] = lines[-1].strip() + '\n'
    # 末尾添加新增的修饰
    num = int(num) + 1 
    ion_num = 0 

    # print(mass_diff_pair_rank)
    #print(exist_ion_flag_list)
    #print(refine_ion_list)

    for key in mass_diff_pair_rank: 
        if parameter_dict['isotope_labeling'] == 'False':
            short_key = 'PFIND_DELTA_' + key.split('.')[0]
        else:
            short_key = key.split('.')[0]
        
        line1 = 'name' + str(num) + '=' + short_key + ' 0 \n'
        lines.append(line1)
        # 目前只有2种情况，N端的话就是全部；其他则是选择最高频率的氨基酸种类
        if mod_static_dict[key].most_common()[0][0] == 'N-SIDE': 
            line2 = short_key + '=ABCDEFGHIJKLMNOPQRSTUVWXYZ PEP_N ' + str(mass_diff_dict[key])\
                + ' ' + str(mass_diff_dict[key])
        else:
            line2 = short_key + '=' + mod_static_dict[key].most_common()[0][0] + ' NORMAL ' +  str(mass_diff_dict[key])\
                + ' ' + str(mass_diff_dict[key]) 
             
            #line2 = key + '=ABCDEFGHIJKLMNOPQRSTUVWXYZ NORMAL ' +  str(mass_diff_dict[key])\
            #    + ' ' + str(mass_diff_dict[key]) 
        # 是否需要加入中性丢失 
        
        if exist_ion_flag_list is not None and exist_ion_flag_list[ion_num//2] == True and ion_num < len(refine_ion_list) * 2 and len(refine_ion_list[ion_num//2][0]) >= 1: 
            if ion_num % 2 == 0: 
                line2 += ' 1 ' + str(refine_ion_list[ion_num//2][0][0]) + ' ' + str(refine_ion_list[ion_num//2][0][0]) + ' pFindDELTA \n' 
            else:
                line2 += ' 1 ' + str(refine_ion_list[ion_num//2][1][0]) + ' ' + str(refine_ion_list[ion_num//2][1][0]) + ' pFindDELTA \n'
            ion_num += 1 
        else:
            line2 += ' 0 pFindDELTA \n' 
        lines.append(line2)
        num += 1
    new_ini_path = os.path.join(current_path, 'modification-new.ini')
    with open(new_ini_path, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line)
    return new_ini_path 

def new_mass_list_generate(mass_diff_pair_rank, mass_diff_list): 
    new_mass_list = []
    for m in mass_diff_pair_rank: 
        for n in mass_diff_list:
            if m in n: 
                new_mass_list.append(n)
                break 
    return new_mass_list

def add_mass_list(mass_diff_pair_rank, parameter_dict):
    new_mass_pair_rank = [] 
    mass_diff = int(parameter_dict['mass_of_diff_diff'])
    for m in mass_diff_pair_rank:
        new_mass_pair_rank.append(m) 
        mm = str(int(m) + mass_diff)
        new_mass_pair_rank.append(mm) 
    return new_mass_pair_rank 


def delete_file(current_path, file_name):
    file_path = os.path.join(current_path, file_name)
    if os.path.exists(file_path): 
        os.remove(file_path)

def remove_file(current_path, file_name, reporting_res_path): 
    file_path = os.path.join(current_path, file_name) 
    target_file_path = os.path.join(reporting_res_path, file_name)
    if os.path.exists(target_file_path):
        os.remove(target_file_path)
    if os.path.exists(file_path): 
        shutil.move(file_path, reporting_res_path)


# 搜索pFind所在路劲
def pfind_path_find(path, target_path, name='pFind.exe.config'):
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path):
            pfind_path_find(item_path, target_path, name) 
        elif os.path.isfile(item_path):
            if name in item: 
                target_path.append(item_path[:-21])
                


# 对结果输出文件去掉第2列
def summary_remove_second_col(path): 
    with open(path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    new_lines = []
    for line in lines: 
        line_list = line.split('\t')
        new_line = '' 
        for i in range(len(line_list)): 
            if i == 1:
                continue  
            new_line +=  line_list[i].strip() + '\t' 
        new_line += '\n' 
        new_lines.append(new_line) 
    with open(path, 'w', encoding='utf-8') as f: 
        for l in new_lines: 
            f.write(l) 

