# 对盲搜鉴定到的修饰，进行质量校正
import os 
import numpy as np 
from collections import Counter 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns 
import math 


element_dict={
    "C": 12.0000000,
    "H": 1.00782503207,
    "Pm": 1.00727647012,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "S": 31.972071
}

amino_acid_dict={
    "A" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "B" : element_dict["C"]*0,
    "C" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "D" : element_dict["C"]*4 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "E" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "F" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "G" : element_dict["C"]*2 + element_dict["H"]*3 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "H" : element_dict["C"]*6 + element_dict["H"]*7 + element_dict["N"]*3 + element_dict["O"]*1 + element_dict["S"]*0,
    "I" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "J" : element_dict["C"]*0,
    "K" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "L" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "M" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "N" : element_dict["C"]*4 + element_dict["H"]*6 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "O" : element_dict["C"]*0,
    "P" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Q" : element_dict["C"]*5 + element_dict["H"]*8 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "R" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*4 + element_dict["O"]*1 + element_dict["S"]*0,
    "S" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "T" : element_dict["C"]*4 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "U" : element_dict["C"]*0,
    "V" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "W" : element_dict["C"]*11 + element_dict["H"]*10 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "X" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Y" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "Z" : element_dict["C"]*0, 
}

h2o_mass = element_dict["H"]*2 + element_dict["O"]*1
proton_mass = element_dict['Pm']


# 读取选择的常见修饰，返回dict
def common_dict_create(current_path):
    modification_path = os.path.join(current_path, 'modification-null.ini')
    with open(modification_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    i = 1
    common_dict = {} 
    while i < len(lines):
        if len(lines[i]) < 2:
            break
        mod_name = lines[i].split()[0]
        eq_idx = mod_name.find('=')
        mod_name = mod_name[eq_idx+1:]
        mod_mass = lines[i+1].split()[2]
        common_dict[mod_name] = float(mod_mass)
        i += 2
    return common_dict


# 计算系统误差(绝对值)
def system_shift_compute(lines, system_correct='mean'):
    mass_shift = []
    for line in lines:
        line = line.split('\t')
        if len(line[10]):
            continue
        else:
            mass_shift.append(float(line[7]))
    if system_correct == 'mean':
        system_shift = np.mean(mass_shift) 
    else:
        system_shift = np.median(mass_shift) 
    return system_shift


# 计算系统误差(ppm) 
def ppm_system_shift_compute(lines, system_correct='mean'):
    mass_shift = []
    for line in lines:
        line = line.split('\t')
        if len(line[10]):
            continue
        else:
            ppm = float(line[7]) * 1000000 / float(line[6]) 
            mass_shift.append(ppm) 
    if system_correct == 'mean': 
        # 去掉大于3倍标准差的值 
        # print(len(mass_shift))
        m = np.mean(mass_shift)
        sd = (sum([(v - m)**2 for v in mass_shift])/(len(mass_shift)-1))**0.5
        new_mass_shift = [v for v in mass_shift if abs(v - m)/sd <= 3]
        system_shift = np.mean(new_mass_shift) 
        # print(len(new_mass_shift))
    else:
        system_shift = np.median(mass_shift) 
    return system_shift


# 计算修饰的精确质量 
def accurate_mass_compute(lines, mass, common_dict, factor_shift, mod_correct='mean'):
    mass_list = []
    for line in lines:
        if mass not in line:
            continue 
        line = line.split('\t') 
        mod_list = line[10].split(';')[:-1] 
        #if len(mod_list) > 1:
        #    continue
        parent_mass = float(line[2]) * factor_shift
        sequence = line[5]
        amino_mass = 0.0
        for a in sequence: 
            if a in amino_acid_dict.keys():
                amino_mass += amino_acid_dict[a]
        mod_mass = parent_mass - amino_mass - proton_mass - h2o_mass
        if len(mod_list) > 1:
            for mod in mod_list:
                if mass in mod:
                    continue
                mod = mod.split(',')[1]
                mod_mass -= common_dict[mod]
        mass_list.append(mod_mass)
    
    if len(mass_list) == 0: 
        return 0.0 
    else:
        if mod_correct == 'mean': 
            return np.mean(mass_list) 
        else: 
            return np.median(mass_list)



# 计算分数加权后修饰的精确质量 
def weight_accurate_mass_compute(lines, mass, common_dict):
    mass_sum = 0.0 
    score_sum = 0.0 
    for line in lines:
        if mass not in line:
            continue
        line = line.split('\t')
        mod_list = line[10].split(';')[:-1]
        score = math.log10(float(line[9]))
        parent_mass = float(line[2])
        sequence = line[5]
        amino_mass = 0.0
        for a in sequence: 
            if a in amino_acid_dict.keys():
                amino_mass += amino_acid_dict[a]
        mod_mass = parent_mass - amino_mass - proton_mass - h2o_mass
        if len(mod_list) > 1:
            for mod in mod_list:
                if mass in mod:
                    continue
                mod = mod.split(',')[1]
                mod_mass -= common_dict[mod]
        mass_sum += mod_mass * score
        score_sum += score 
    # return np.mean(mass_list) 
    if score_sum == 0.0: 
        return 0.0
    else: 
        return mass_sum / score_sum 


# 只选择FDR千分之一的谱图
def q_value_filter(lines): 
    for i in range(len(lines)): 
        q_value = float(lines[i].split('\t')[4]) 
        if q_value > 0.001: 
            return lines[:i]
    return lines 


# system_correct={mean, median}, mod_correct={mean, median, weight}
def mass_correct(current_path, blind_path, mass_diff_list, system_correct='mean', mod_correct='mean'):
    # 读取常见修饰列表
    common_dict = common_dict_create(current_path)

    # 读取盲搜鉴定结果文件
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    lines = lines[1:] 
    
    # 只选择q-value < 0.001的谱图 
    # lines = q_value_filter(lines) 

    # 计算系统误差
    system_shift = ppm_system_shift_compute(lines, system_correct) 
    factor_shift = 1.0 / (1.0 + system_shift/ 1000000.0) 

    # 计算未知修饰的精确质量
    mass_dict = {}
    for mass_diff in mass_diff_list: 
        if mod_correct == 'weight': 
            accurate_mass_diff = weight_accurate_mass_compute(lines, mass_diff, common_dict) 
        else: 
            accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict, factor_shift, mod_correct) 
        # 使用绝对值做校准
        # mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff - system_shift)) 

        # 使用ppm做校准 
        # print(accurate_mass_diff)
        mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff))
    return mass_dict



# 删除质量小于200Da的偏差和负数的修饰 
def small_delta_filter(mass_difference_list, min_mass_modification):
    name2mass = {}
    new_mass_diff_list = []
    for mass_diff in mass_difference_list:
        mass = mass_diff.split('_')[-1]
        if mass[0] == '-':
            continue
        mass = float(mass)
        if mass < min_mass_modification:
            continue
        name2mass[mass_diff] = mass
        new_mass_diff_list.append(mass_diff)
    return name2mass, new_mass_diff_list  



# 筛选修饰质量差满足设定的插值的修饰 
# mass_diff_diff = 6.0201 
# 应该是20ppm，但是实际太小完全没有删选度
def mass_diff_diff_filter(name2mass, mass_diff_list, mass_diff_diff):
    refined_list = []
    for i in range(len(mass_diff_list)):
        for j in range(i+1, len(mass_diff_list)):
            mass_left = name2mass[mass_diff_list[i]]
            mass_right = name2mass[mass_diff_list[j]]
            delta_mass = abs(mass_right -  mass_left)
            if abs(delta_mass - mass_diff_diff) < 0.1:
                refined_list.append(mass_diff_list[i])
                refined_list.append(mass_diff_list[j])
    new_refined_list = []
    for mod_name in refined_list:
        flag_rep = False 
        for ref_mod in new_refined_list:
            if abs(name2mass[mod_name] - name2mass[ref_mod]) < 0.05:
                flag_rep = True
        if flag_rep == False:
            new_refined_list.append(mod_name)
    return new_refined_list 



# 统计修饰发生的位点分布 
def mass_static(blind_path, current_path, mass_diff_list): 
    mod_position_dict = {} 
    mod_number_dict = {} 
    for mass in mass_diff_list:
        mod_position_dict[mass] = [] 
        mod_number_dict[mass] = 0 
    
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    spectra_num = len(lines)
    for i in range(1, spectra_num):
        if len(lines[i]) < 4:
            break
        sequence = lines[i].split('\t')[5]
        mod_list = lines[i].split('\t')[10].split(';')[:-1]
        for mod in mod_list:
            pos, mod_name = mod.split(',')
            if mod_name in mass_diff_list:
                pos = int(pos) 
                mod_number_dict[mod_name] += 1
                if pos == 0 or pos == 1:
                    mod_position_dict[mod_name].append(sequence[0])
                    mod_position_dict[mod_name].append('N-SIDE')
                elif pos >= len(sequence):
                    mod_position_dict[mod_name].append(sequence[-1])
                    mod_position_dict[mod_name].append('C-SIDE')
                else:
                    mod_position_dict[mod_name].append(sequence[pos-1])
    
    mod_static_dict = {}
    for mod_name in mass_diff_list:
        counter = Counter(mod_position_dict[mod_name])
        mod_static_dict[mod_name] = counter 
    
    return mod_static_dict, mod_number_dict  



# 写pchem.summary结果文件
# mass_diff_list 修饰名称的列表  mod_static_dict 修饰名字-> Counter
# mod_number_dict 修饰名字 -> PSM 
# mass_diff_dict 修饰名字 -> 精准质量 
# mod2pep 修饰名字 -> peptide 
# simple_dict 修饰名字 -> 简化版 
def summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict): 
    mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, rank_dict, label_dict, mass_diff_rank, _  = \
        mass_refine(mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict)
    lines = []
    lines.append('Rank \tLabel \tMass modification \tPeptide \tPSM  \tSite (Location Score/Probability)-Top1 \tSite (Location Score/Probability)-Others \tAccurate mass \n') 

    # 按照出现频率进行排序 
    for mod in mass_diff_rank:
        local_list = mod_static_dict[mod].most_common() 
        Top1 = local_list[0][0] + '(' + str(round(local_list[0][1]/mod_number_dict[mod],3)) + '); ' 
        Others = ''
        for j in range(1, len(local_list)):
            Others += local_list[j][0] + '(' + str(round(local_list[j][1]/mod_number_dict[mod],3)) + '); ' 
        line = str(rank_dict[mod]) + '\t' + label_dict[mod] + '\t' + 'PFIND_DELTA_' + mod + '\t' + str(mod2pep[mod]) + '\t' + str(mod_number_dict[mod]) + '\t' \
            + Top1 + '\t' + Others + '\t' + str(mass_diff_dict[mod]) + '\n'
        lines.append(line)
    with open(os.path.join(current_path, 'pChem.summary'), 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line)
    # print(mod_static_dict)



# 只输出轻标修饰到结果文件 
# 根据mass_diff_pair_pair输出结果
def new_summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict): 
    mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, unimod_list, ppm_error_dict = \
        mass_refine(mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict) 
    lines = [] 
    lines.append('Rank \tMass Shift \tIsotopic Label \tPeptide Total \tPSM Total \tPeptide L|H \tPSM L|H \
                \tTop1 Site \tTop1 Probability \tAccurate Mass \tSite (Location /Probability)-Others \tAdditional Modification\n')
    
    idx = 1
    for mass_pair in mass_diff_pair_rank: 
        light_mod = mass_pair[0]
        heavy_mod = mass_pair[1] 
        local_list = mod_static_dict[light_mod].most_common() 
        Top1_site = local_list[0][0]
        Top1_pro =  str(round(local_list[0][1]/mod_number_dict[light_mod],3)) 
        Others = '' 
        for j in range(1, len(local_list)):
            Others += local_list[j][0] + '(' + str(round(local_list[j][1]/mod_number_dict[light_mod],3)) + '); ' 
        line = str(idx) + '\t' + 'PFIND_DELTA_' + light_mod + '\t' + 'Yes' + '\t' + str(mod2pep[light_mod] + mod2pep[heavy_mod]) + '\t' + \
            str(mod_number_dict[light_mod] + mod_number_dict[heavy_mod]) + '\t' + str(mod2pep[light_mod]) + '|' + str(mod2pep[heavy_mod]) + '\t' + \
            str(mod_number_dict[light_mod]) + '|' + str(mod_number_dict[heavy_mod]) + '\t' + Top1_site + '\t' + Top1_pro + '\t' + str(mass_diff_dict[light_mod]) + '|' + str(mass_diff_dict[heavy_mod]) + '(' + str(int(ppm_error_dict[light_mod])) + ' ppm)' + \
            '\t' + Others + '\t' + unimod_list[idx-1] + '\n' 
        lines.append(line) 
        idx += 1 
    with open(os.path.join(current_path, 'pChem.summary'), 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line) 
    # 同时保存热力图 
    heat_map_plot(current_path, mass_diff_pair_rank, mod_static_dict, mod_number_dict)



# 绘制结果热力图
def heat_map_plot(current_path, mass_diff_pair_rank, mod_static_dict, mod_number_dict): 
    y_stick = ["N-SIDE", "C-SIDE", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "X", "Y"] 
    y_dict = {}
    for i in range(len(y_stick)):
        y_dict[y_stick[i]] = i 
    
    x_stick = [mod[0] for mod in mass_diff_pair_rank] 
    mod_map = {}
    
    for i in range(len(x_stick)):
        freq_list = [0.0] * len(y_stick) 
        local_list = mod_static_dict[x_stick[i]].most_common()
        for j in range(len(local_list)): 
            cur_position = local_list[j][0]
            cur_freq = round(local_list[j][1]/mod_number_dict[x_stick[i]],3)
            freq_list[y_dict[cur_position]] = cur_freq 
        mod_map[x_stick[i]] = freq_list 
    
    pd_mod_map = pd.DataFrame(mod_map, index=y_stick, columns=x_stick) 
    # print(pd_mod_map)
    ax = sns.heatmap(pd_mod_map, vmin=0.0, vmax=1.0, cmap='YlGnBu', annot=True, annot_kws={"size":4})
    plt.ylabel('amino acid selectivity')
    plt.xlabel('probes')
    png_path = os.path.join(current_path, 'heat_map.png') 
    plt.savefig(png_path, dpi=200)
    


# 保留整数，其余信息进行合并 
def mass_refine(mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict): 
    new_mod_static_dict, new_mod_number_dict, new_mod2pep, new_mass_diff_dict = {}, {}, {}, {}
    new_mass_list = [] 
    int_mass_list = []
    for mass in mod_static_dict: 
        if mass_diff_dict[mass] < 200.0: 
            continue 
        simple_mass = str(int(float(mass.split('_')[2]))) 
        # 不是第一次出现
        if simple_mass in new_mass_list: 
            new_mod_static_dict[simple_mass] = new_mod_static_dict[simple_mass] + mod_static_dict[mass] 
            new_mod_number_dict[simple_mass] = new_mod_number_dict[simple_mass] + mod_number_dict[mass]
            new_mod2pep[simple_mass] = new_mod2pep[simple_mass] + int(mod2pep[mass]) 
        else:
            new_mass_list.append(simple_mass)
            int_mass_list.append(int(simple_mass))
            new_mod_static_dict[simple_mass] = mod_static_dict[mass] 
            new_mod_number_dict[simple_mass] = mod_number_dict[mass]
            new_mod2pep[simple_mass] = int(mod2pep[mass]) 
            new_mass_diff_dict[simple_mass] = mass_diff_dict[mass] 
    # rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank = label_determine(new_mod_number_dict, int_mass_list, int(parameter_dict['mass_diff_diff']))
    rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, ppm_error_dict = accurate_label_determine(new_mod_number_dict, new_mass_diff_dict, parameter_dict)
    unimod_list = unimod_match(unimod_dict, mass_diff_pair_rank, new_mass_diff_dict)
    return new_mod_static_dict, new_mod_number_dict, new_mod2pep, new_mass_diff_dict, rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, unimod_list, ppm_error_dict   


# 精确质量法确定轻重标记
# parameter_dict['mass_diff_diff'] = 6.020132
# parameter_dict['mass_diff_diff_range'] = 100 
def accurate_label_determine(mod_number_dict, mass_diff_dict, parameter_dict): 
    rank_dict = {} 
    label_dict = {} 
    ppm_error_dict = {} 
    # 按照出现频率进行排序 
    rank_tuple = sorted(mod_number_dict.items(), key=lambda x: -x[1]) 
    # 将数字化修饰按照频率从高到低  
    freq_list = [int(mod[0]) for mod in rank_tuple] 
    # 记录已经配对成功的修饰
    mass_diff_rank = [] 
    mass_diff_pair_rank = [] 
    id = 1 
    for i in range(len(freq_list)):
        mod = rank_tuple[i][0]
        if mod_number_dict[mod] < 4:
            break 
        if mod in mass_diff_rank: 
            continue 
        # 只报出找到匹配的修饰对
        for j in range(i+1, len(freq_list)): 
            if rank_tuple[j][0] in mass_diff_rank:
                continue 
            pair_mod = rank_tuple[j][0] 
            ppm_error = ppm_calculate(mass_diff_dict[mod], mass_diff_dict[pair_mod])
            if ppm_error < parameter_dict['mass_diff_diff_range']*1.5: 
                if mass_diff_dict[mod] < mass_diff_dict[pair_mod]: 
                    label_dict[mod] = 'L'
                    label_dict[pair_mod] = 'H'
                    mass_diff_pair_rank.append([mod, pair_mod])
                    ppm_error_dict[mod] = ppm_error 
                else:
                    label_dict[mod] = 'H'
                    label_dict[pair_mod] = 'L'
                    mass_diff_pair_rank.append([pair_mod, mod]) 
                    ppm_error_dict[pair_mod] = ppm_error
                
                mass_diff_rank.append(mod)
                mass_diff_rank.append(pair_mod) 
                rank_dict[mod] = id 
                rank_dict[pair_mod] = id 
                id += 1 
                break 
    return rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, ppm_error_dict  


def ppm_calculate(a, b):
    return abs(abs(b-a)-6.020131)/6.020132*1000000 


# 根据质量差判断未知修饰之间的可能
def unimod_match(unimod_dict, mass_diff_pair_rank, mass_diff_dict): 
    unimod_list = []
    unimod_list.append('') 
    if len(mass_diff_pair_rank) <= 1:
        return unimod_list 
    baseline = mass_diff_pair_rank[0][0] 
    baseline_mass = mass_diff_dict[baseline]
    for mass_diff_pair in mass_diff_pair_rank[1:]:
        t_unimod = ''
        cur_mod = mass_diff_pair[0]
        cur_mod_mass = mass_diff_dict[cur_mod]
        target_diff = baseline_mass - cur_mod_mass 
        for k, v in unimod_dict.items():
            if abs(v - target_diff) < 0.001: 
                t_unimod += k + ';' 
        target_diff = cur_mod_mass - baseline_mass 
        for k, v in unimod_dict.items():
            if abs(v - target_diff) < 0.001: 
                t_unimod += k + ';'
        unimod_list.append(t_unimod)
    return unimod_list 


# 找到轻重标记 
def label_determine(mod_number_dict, int_mass_list, mass_diff): 
    rank_dict = {} 
    label_dict = {} 
    # 按照出现频率进行排序 
    rank_tuple = sorted(mod_number_dict.items(), key=lambda x: -x[1]) 
    
    # 将数字化修饰按照频率从高到低  
    freq_list = [int(mod[0]) for mod in rank_tuple]
    
    mass_diff_rank = [] 
    mass_diff_pair_rank = [] 
    id = 1 
    for i in range(len(freq_list)): 
        # 小于3的数目不显示 
        mod = rank_tuple[i][0]
        if mod_number_dict[mod] < 4:
            break 
        if mod in mass_diff_rank:
            continue
        rank_dict[mod] = id 
        mass_diff_rank.append(mod) 
        mod_num = int(mod) 
        
        find_pair = False 
        for j in range(i+1, len(freq_list)): 
            if rank_tuple[j][0] in mass_diff_rank:
                continue 
            if mod_num - freq_list[j] == mass_diff: 
                label_dict[mod] = 'H' 
                light_mod = rank_tuple[j][0]
                rank_dict[light_mod] = id  
                label_dict[light_mod] = 'L' 
                mass_diff_rank.append(light_mod) 
                mass_diff_pair = [light_mod, mod]
                mass_diff_pair_rank.append(mass_diff_pair)
                find_pair = True 
                break 
            if freq_list[j] - mod_num == mass_diff:
                label_dict[mod] = 'L'
                heavy_mod = rank_tuple[j][0] 
                rank_dict[heavy_mod] = id  
                label_dict[heavy_mod] = 'H' 
                mass_diff_rank.append(heavy_mod) 
                mass_diff_pair = [mod, heavy_mod]
                mass_diff_pair_rank.append(mass_diff_pair)
                find_pair = True 
                break 

        if find_pair == False:
            label_dict[mod] = ' ' 
            # mass_diff_pair_rank.append([mod])
        id += 1

        '''
        # 使用字典直接搜的话，会出现高频优先匹配到小质量低频，出现错误匹配。
        if (mod_num - mass_diff) in int_mass_list and str(mod_num - mass_diff) not in mass_diff_rank: 
            label_dict[mod] = 'H'
            light_mod = str(mod_num - mass_diff) 
            rank_dict[light_mod] = i 
            label_dict[light_mod] = 'L'
            mass_diff_rank.append(light_mod) 
            i += 1 
            continue 
        if (mod_num + mass_diff) in int_mass_list and str(mod_num + mass_diff) not in mass_diff_rank:
            label_dict[mod] = 'L' 
            h_mod = str(mod_num + mass_diff) 
            label_dict[h_mod] = 'H' 
            rank_dict[h_mod] = i 
            mass_diff_rank.append(h_mod) 
            i += 1 
            continue 
        label_dict[mod] = ' ' 
        i += 1 
        '''
    # print(mass_diff_rank)
    # print(rank_dict) 
    # print(mass_diff_pair_rank)
    return rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank  



# 选择不重复的topk来做限定式搜索
def mass_select(mass_diff_list, k, name2mass):
    selected_list = []
    i = 0 
    for mod in mass_diff_list: 
        flag = True
        for reference in selected_list:
            if abs(name2mass[mod]-name2mass[reference]) < 0.03:
                flag = False
        if flag == True:
            selected_list.append(mod)
            i += 1 
        if i >= k:
            break
    return selected_list 
        

# 读取unimod的质量，返回字典 
def unimod_dict_generate(modification_dict): 
    unimod_dict = {}
    for key in modification_dict.keys():
        unimod_mass = float(modification_dict[key].split('\n')[1].split()[2]) 
        unimod_name = key.split('[')[0]
        unimod_dict[unimod_name] = unimod_mass  
    return unimod_dict 


if __name__ == "__main__": 
    P = 'NAHSATTWSGQYVGGAEAR' 
    mass = 0.0 
    for a in P:
        mass += amino_acid_dict[a]
    mass += proton_mass + h2o_mass
    print(mass)