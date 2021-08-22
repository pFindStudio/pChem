import os 

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
        if 'mass_diff_diff' in lines[i]:
            parameter_dict['mass_diff_diff'] = float(parameter_pick(lines[i])) 
        if 'common_modification_number' in lines[i]:
            parameter_dict['common_modification_number'] = int(parameter_pick(lines[i]))
        if 'close_mass_diff_number' in lines[i]:
            parameter_dict['close_mass_diff_number'] = int(parameter_pick(lines[i])) 
        if 'min_mass_modification' in lines[i]:
            parameter_dict['min_mass_modification'] = float(parameter_pick(lines[i]))
        if 'mass_diff_diff_range' in lines[i]:
            parameter_dict['mass_diff_diff_range'] = int(parameter_pick(lines[i]))
    return parameter_dict


# 写open参数文件
def open_cfg_write(cfg_path, parameter_dict):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    # 指定参数内容修改
    new_lines = []
    for i in range(len(lines)):
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
def close_cfg_write(cfg_path, current_path, parameter_dict, mass_diff_dict):
    with open(cfg_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    mod_line, flag = combine_common_list(current_path, parameter_dict, mass_diff_dict)

    # 指定参数内容修改 
    new_lines = [] 
    for i in range(len(lines)):
        if 'selectmod' in lines[i]:
            lines[i] = parameter_modify(lines[i], mod_line)
        if 'fixmod' in lines[i] and flag == True:
            lines[i] = parameter_modify(lines[i],'Carbamidomethyl[C]')
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



# 盲搜之后生成合并的修饰列表
def combine_common_list(current_path, parameter_dict, mass_diff_dict):
    mod_line = ""
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
    for key in mass_diff_dict.keys():
        mod_line += key + ';'
    
    flag = False 
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
def expand_modification_ini(mass_diff_dict, mod_static_dict, current_path, ini_path):
    with open(ini_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    while True:
        if len(lines[-1]) < 4:
            lines = lines[:-1]
        else:
            break
    # 修改第一行总数
    lines[0], num = lines[0].split('=')
    total_num = int(num) + len(mass_diff_dict)
    lines[0] = lines[0] + '=' + str(total_num) + '\n' 
    lines[-1] = lines[-1].strip() + '\n'
    # 末尾添加新增的修饰
    num = int(num) + 1
    for key in mass_diff_dict:
        line1 = 'name' + str(num) + '=' + key + ' 0 \n'
        lines.append(line1)
        # 目前只有2种情况，N端的话就是全部；其他则是选择最高频率的氨基酸种类
        if mod_static_dict[key].most_common()[0][0] == 'N-SIDE': 
            line2 = key + '=ACDEFGHIKLMNPQRSTVWXY PEP_N ' + str(mass_diff_dict[key])\
                + ' ' + str(mass_diff_dict[key]) + ' 0 pFindDELTA \n' 
        else:
            line2 = key + '=' + mod_static_dict[key].most_common()[0][0] + ' NORMAL ' +  str(mass_diff_dict[key])\
                + ' ' + str(mass_diff_dict[key]) + ' 0 pFindDELTA \n' 
        lines.append(line2)
        num += 1
    new_ini_path = os.path.join(current_path, 'modification-new.ini')
    with open(new_ini_path, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line)
    return new_ini_path



