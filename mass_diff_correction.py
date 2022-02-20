# 对盲搜鉴定到的修饰，进行质量校正
import os 
import numpy as np 
from collections import Counter 
import matplotlib.pyplot as plt 
import pandas as pd 
import seaborn as sns 
import math 
from ion_type_learning import ion_type_determine, ion_filter  
from parameter import element_dict, amino_acid_dict, common_dict_create, h2o_mass, proton_mass
from matplotlib.backends.backend_pdf import PdfPages

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



# 生成所有未知修饰的质量列表
def origin_mass_list_generate(lines, common_dict, factor_shift): 
    origin_mass_list = [] 
    for line in lines: 
        if 'PFIND_DELTA' not in line:
            continue 
        line = line.split('\t') 
        mod_list = line[10].split(';')[:-1] 
        # if len(mod_list) > 1:
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
                mod = mod.split(',')[1]
                if 'PFIND_DELTA' in mod:
                    continue
                mod_mass -= common_dict[mod]
        origin_mass_list.append(float('%.6f'%(mod_mass))) 
    return origin_mass_list 



# 通过多次迭代优化 
def iterative_mass_compute(mass_diff, original_unknown_list, light_heavy_dict=None, mass_diff_diff=None, mode='blind'): 
    # print(mass_diff_diff)
    if mode == 'blind':
        target_mass = float(mass_diff.split('_')[2]) 
    else:
        target_mass = mass_diff 
    #if mass_diff in light_heavy_dict.keys():
    #    mode = light_heavy_dict[mass_diff] 
    

    data = [] 
    for mass in original_unknown_list: 
        if mass >= target_mass - 0.01 and mass <= target_mass + 0.01: 
            data.append(mass) 
        #if mode == 'light':
        #    if mass - mass_diff_diff >= target_mass - 0.01 and mass - mass_diff_diff <= target_mass + 0.01: 
        #        data.append(mass- mass_diff_diff) 
        #if mode == 'heavy':
        #    if mass + mass_diff_diff >= target_mass - 0.01 and mass + mass_diff_diff <= target_mass + 0.01: 
        #        data.append(mass + mass_diff_diff) 
    print(data)
    spectrum_num = len(data)
    if len(data) == 0: 
        return 0.0, 0
    mu0 = np.mean(data) 
    times = 1 

    while True: 
        if len(data) <= 100:
            break
        a = [] 
        for mass in original_unknown_list:
            if mass >= mu0 - 0.01 and mass <= mu0 + 0.01: 
                a.append(mass) 
            #if mode == 'light':
            #    if mass - mass_diff_diff >= mu0 - 0.01 and mass - mass_diff_diff <= mu0 + 0.01: 
            #        data.append(mass- mass_diff_diff) 
            #if mode == 'heavy':
            #    if mass + mass_diff_diff >= mu0 - 0.01 and mass + mass_diff_diff <= mu0 + 0.01: 
            #        data.append(mass + mass_diff_diff) 
        spectrum_num = max(spectrum_num, len(a))
        mu, sigma = np.mean(a), np.std(a) 
        data = a 
        if mu == mu0 or times == 3: 
            break 
        mu0 = mu
        times += 1 
    return mu0, spectrum_num



# 计算修饰的精确质量 
def accurate_mass_compute(lines, mass, common_dict, factor_shift, mod_correct='mean'):
    mass_list = []
    for line in lines:
        if mass not in line:
            continue 
        line = line.split('\t') 
        mod_list = line[10].split(';')[:-1] 
        if len(mod_list) > 1:
            continue
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


# 改写pfind结果文件 
def pfind_result_rewrite(blind_path, origin_lines, common_dict, factor_shift): 
    new_lines = []
    new_lines.append(origin_lines[0][:-1] + '\t Accurate modification mass \n')
    for line in origin_lines[1:]: 
        line_list = line.split('\t') 
        mod_list = line_list[10]
        if 'PFIND_DELTA' in mod_list: 
            mod_list = mod_list.split(';')[:-1] 
            parent_mass = float(line_list[2]) * factor_shift
            sequence = line_list[5]
            amino_mass = 0.0
            for a in sequence: 
                if a in amino_acid_dict.keys():
                    amino_mass += amino_acid_dict[a]
            mod_mass = parent_mass - amino_mass - proton_mass - h2o_mass
            if len(mod_list) > 1:
                for mod in mod_list:
                    if 'PFIND_DELTA' in mod:
                        continue
                    mod = mod.split(',')[1]
                    mod_mass -= common_dict[mod] 
            new_lines.append(line[:-1] + '\t' + str(float('%.6f'%(mod_mass))) + '\n') 
        else:
            new_lines.append(line) 
    with open(blind_path, 'w', encoding='utf-8') as f: 
        for line in new_lines:
            f.write(line)




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
def mass_correct(current_path, blind_path, mass_diff_list, parameter_dict, system_correct='mean', mod_correct='mean', pattern='open'):
    # 读取常见修饰列表
    common_dict = common_dict_create(current_path)

    # 读取盲搜鉴定结果文件
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    origin_lines = lines 
    lines = lines[1:] 
    
    # 只选择q-value < 0.001的谱图 
    # lines = q_value_filter(lines) 

    # 计算系统误差 
    # 是否需要引入误差校准
    if parameter_dict['mass_calibration'] == 'True': 
        system_shift = ppm_system_shift_compute(lines, system_correct) 
        factor_shift = 1.0 / (1.0 + system_shift/ 1000000.0) 
    else:
        factor_shift = 1.0 
        print('no mass calibration is employed')
    
    # 计算未知修饰的精确质量 
    original_unknown_list = origin_mass_list_generate(lines, common_dict, factor_shift)
    mass_dict = {} 
    mod_number_dict = {}
    # 如果是盲搜，先做一轮粗删选，方便加速 -> 新版是需要的
    if pattern != 'close': 
        new_mass_diff_list, light_heavy_dict = coarse_filter(mass_diff_list, parameter_dict) 
        # print(new_mass_diff_list, light_heavy_dict) 
        mass_diff_list = new_mass_diff_list
    for mass_diff in mass_diff_list: 
        if mod_correct == 'weight': 
            accurate_mass_diff = weight_accurate_mass_compute(lines, mass_diff, common_dict) 
        else: 
            if pattern == 'close': 
                accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict, factor_shift, mod_correct)
            else:
                accurate_mass_diff, spectrum_num = iterative_mass_compute(mass_diff, original_unknown_list) 
                mod_name = str(int(float(mass_diff.split('_')[2]))) 
                if mod_name in mod_number_dict.keys(): 
                    mod_number_dict[mod_name] += spectrum_num 
                else:
                    mod_number_dict[mod_name] = spectrum_num
                # print(accurate_mass_diff)
                # accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict, factor_shift, mod_correct) 
        # 使用绝对值做校准
        # mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff - system_shift)) 

        # 使用ppm做校准 
        # print(accurate_mass_diff)
        mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff)) 
    # print(mass_dict)
    pfind_result_rewrite(blind_path, origin_lines, common_dict, factor_shift)
    return mass_dict, mod_number_dict


# system_correct={mean, median}, mod_correct={mean, median, weight}
# mass_diff_blind_dict['PFIND_DELTA_252.12'] = 252.12222
# total_name_dict['PFIND_DELTA_252'] = 'PFIND_DELTA_252.12'
def close_mass_correct(current_path, blind_path, mass_diff_list, parameter_dict, mass_diff_blind_dict, total_name_dict):
    # 读取常见修饰列表
    common_dict = common_dict_create(current_path)

    # 读取盲搜鉴定结果文件
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    origin_lines = lines 
    lines = lines[1:] 
    
    # 只选择q-value < 0.001的谱图 
    # lines = q_value_filter(lines) 

    # 计算系统误差 
    # 是否需要引入误差校准
    if parameter_dict['mass_calibration'] == 'True': 
        system_shift = ppm_system_shift_compute(lines) 
        factor_shift = 1.0 / (1.0 + system_shift/ 1000000.0) 
    else:
        factor_shift = 1.0 
        print('no mass calibration is employed')
    
    # 计算未知修饰的精确质量 
    original_unknown_list = origin_mass_list_generate(lines, common_dict, factor_shift)
    mass_dict = {} 
    mod_number_dict = {}
    # 如果是盲搜，先做一轮粗删选，方便加速 
    #if pattern != 'close': 
    #    new_mass_diff_list, light_heavy_dict = coarse_filter(mass_diff_list, parameter_dict) 
        # print(new_mass_diff_list, light_heavy_dict) 
    #    mass_diff_list = new_mass_diff_list
    for mass_diff in mass_diff_list: 
        print(mass_diff)
        #if mod_correct == 'weight': 
        #    accurate_mass_diff = weight_accurate_mass_compute(lines, mass_diff, common_dict) 
        #else: 
        #    if pattern == 'close': 
        #        accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict, factor_shift, mod_correct)
        #    else: 
        total_name = total_name_dict[mass_diff] 
        target_mass = mass_diff_blind_dict[total_name]
        accurate_mass_diff, spectrum_num = iterative_mass_compute(target_mass, original_unknown_list, mode='close') 
        mod_name = str(int(float(mass_diff.split('_')[2]))) 
        if mod_name in mod_number_dict.keys(): 
            mod_number_dict[mod_name] += spectrum_num 
        else:
            mod_number_dict[mod_name] = spectrum_num
                # print(accurate_mass_diff)
                # accurate_mass_diff = accurate_mass_compute(lines, mass_diff, common_dict, factor_shift, mod_correct) 
        # 使用绝对值做校准
        # mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff - system_shift)) 

        # 使用ppm做校准 
        
        mass_dict[mass_diff] = float('%.6f'%(accurate_mass_diff)) 
    # print(mass_dict)
    pfind_result_rewrite(blind_path, origin_lines, common_dict, factor_shift)
    return mass_dict, mod_number_dict



# 根据质量差先做一轮筛选 
def coarse_filter(mass_diff_list, parameter_dict):
    light_heavy_dict = {} 
    new_mass_diff_list = [] 
    delta_mass = int(parameter_dict['mass_of_diff_diff']) 
    int_mass_list = [] 
    light_list = []
    heavy_list = []
    compet_list = []
    for mass in mass_diff_list: 
        int_mass_list.append(int(float(mass.split('_')[2]))) 
    for mass in mass_diff_list: 
        int_mass = int(float(mass.split('_')[2]))  
        if int_mass in compet_list:
            continue 
        if int_mass in light_list:
            compet_list.append(int_mass)
            new_mass_diff_list.append(mass)
            light_heavy_dict[mass] = 'light'
            continue
        if int_mass in heavy_list:
            compet_list.append(int_mass)
            new_mass_diff_list.append(mass)
            light_heavy_dict[mass] = 'heavy'
            continue 
        if int_mass + delta_mass in int_mass_list: 
            heavy_list.append(int_mass + delta_mass) 
            compet_list.append(int_mass)
            new_mass_diff_list.append(mass)
            light_heavy_dict[mass] = 'light' 
            continue 
        if int_mass - delta_mass in int_mass_list:
            light_list.append(int_mass - delta_mass)
            compet_list.append(int_mass)
            light_heavy_dict[mass] = 'heavy'
            new_mass_diff_list.append(mass)
            continue
    return new_mass_diff_list, light_heavy_dict 



# 删除质量小于200Da的偏差和负数的修饰 
def small_delta_filter(mass_difference_list, parameter_dict):
    name2mass = {}
    new_mass_diff_list = []
    for mass_diff in mass_difference_list:
        mass = mass_diff.split('_')[-1]
        if mass[0] == '-':
            continue
        mass = float(mass)
        if mass < parameter_dict['min_mass_modification'] or mass > parameter_dict['max_mass_modification']:
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
def mass_static(blind_path, current_path, mass_diff_list, side_position='False'): 
    mod_position_dict = {} 
    mod_number_dict = {} 
    if side_position == 'False': 
        side_flag = False
    else:
        side_flag = True
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
                    if side_flag == True:
                        mod_position_dict[mod_name].append('N-SIDE')
                elif pos >= len(sequence):
                    mod_position_dict[mod_name].append(sequence[-1]) 
                    if side_flag == True:
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


def mod_number_update(origin_dict, new_dict): 
    combine_dict = {} 
    for k in origin_dict.keys(): 
        combine_dict[k] = max(origin_dict[k], new_dict[k]) 
    return combine_dict 


# 根据mass_diff_pair_pair输出结果
def new_summary_write(current_path, mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict, sim_mod_number=None, pattern='blind'): 
    
    
    mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, unimod_list, ppm_error_dict = \
        mass_refine(mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict, pattern) 
    # print(mod_number_dict)
    
    new_mod_number_dict = mod_number_update(mod_number_dict, sim_mod_number) 
    mod_number_dict = sim_mod_number
    lines = [] 
    lines.append('Rank \tPDM \tAccurate Mass \
                \tTop1 Site|Probability \t Others \t #PSM \t #PSM L|H \n')
    # print('matched:', mass_diff_pair_rank)
    idx = 1 
    
    prob_dict = {}
    for mass_pair in mass_diff_pair_rank: 
        light_mod = mass_pair[0]
        heavy_mod = mass_pair[1] 
        local_list = mod_static_dict[light_mod].most_common() 
        
        prob_dict[light_mod] = []
        Top1_site = local_list[0][0]
        
        Others = '' 
        residual_pro = 0.0 
        factor = (new_mod_number_dict[light_mod] - mod_number_dict[light_mod] + local_list[0][1]) * mod_number_dict[light_mod] / new_mod_number_dict[light_mod] /  local_list[0][1] 

        if Top1_site == 'N-SIDE' or Top1_site == 'C-SIDE': 
            parent_num = 0 
            for j in range(1, len(local_list)): 
                parent_num += local_list[j][1] 
        else: 
            parent_num = new_mod_number_dict[light_mod] 

        for j in range(1, len(local_list)):
            temp_pro = local_list[j][1] / parent_num 
            if local_list[j][0] == 'N-SIDE' or local_list[j][0] == 'C-SIDE':
                temp_pro = temp_pro * parent_num * factor / mod_number_dict[light_mod] 
            prob_dict[light_mod].append([local_list[j][0], temp_pro])
            Others += local_list[j][0] + '(' + str(round(temp_pro, 3)) + '); '
            if local_list[j][0] == 'N-SIDE' or local_list[j][0] == 'C-SIDE': 
                continue
            residual_pro += temp_pro  
        
        if Top1_site == 'N-SIDE' or Top1_site == 'C-SIDE': 
            Top1_pro = str(round(local_list[0][1] / parent_num, 3)) 
        else: 
            Top1_pro =  str(round(1 - residual_pro, 3)) 
        prob_dict[light_mod].append([Top1_site, Top1_pro])
        line = str(idx) + '\t' + 'PFIND_DELTA_' + light_mod + '\t'+ str(mass_diff_dict[light_mod]) +'\t' + Top1_site +'|'+ Top1_pro + '\t' + Others + '\t' + \
            str(mod_number_dict[light_mod] + mod_number_dict[heavy_mod]) + '\t' + str(mod_number_dict[light_mod]) + '|' +  str(mod_number_dict[heavy_mod])  + '\n' 
        lines.append(line) 
        idx += 1 
    
    # print(lines)

    if pattern == 'blind':
        summary_path = os.path.join(current_path, 'pChem.summary')
    else:
        summary_path = os.path.join(current_path, 'pChem-close.summary')

    
    with open(summary_path, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line) 
    # 删除PSM 
    filter_mod = [mod[0] for mod in mass_diff_pair_rank] 
    if pattern == 'blind': 
        filter_mod = summary_filter(current_path, parameter_dict, filter_mod, pattern, summary_path, prob_dict) 
    # 同时保存热力图 
    if len(lines) < 2: 
        print('The number of unknown modification is none, please expand the error range.')
    else: 
        if pattern == 'blind':
            heat_map_plot(current_path, filter_mod, mod_static_dict, mod_number_dict, new_mod_number_dict) 
    
    
    # 对所有的筛选后的未知修饰做中性丢失
    refine_ion_list = [] 
    exist_ion_flag_list = []
    for i in range(len(filter_mod)): 
        # 中性丢失确认，函数输入形式确认 
        modification_list, modification_dict = ion_function_input(filter_mod[i], mass_diff_dict, int(parameter_dict['mass_of_diff_diff']))
        # 计算中性丢失的质量 
        mod2ion, ion_list, exist_ion_flag = ion_type_determine(current_path, modification_list, modification_dict, parameter_dict)
        # print(mod2ion)
        t_refine_ion_list = ion_filter(ion_list, int(parameter_dict['mass_of_diff_diff'])) 
        exist_ion_flag_list.append(exist_ion_flag) 
        refine_ion_list.append(t_refine_ion_list) 
 
    # 将中性丢失结果写入结果文件  
    # print(refine_ion_list, exist_ion_flag_list)
    add_ion_summary(summary_path, refine_ion_list, exist_ion_flag_list)

    return filter_mod, refine_ion_list, exist_ion_flag_list  




def add_ion_summary(summary_path, refine_ion_list, exist_ion_flag_list): 
    # summary_path = os.path.join(current_path, 'pChem.summary')
    with open(summary_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    if 'DFLs' not in  lines[0]: 
        lines[0] = lines[0][:-1] + ' \tDFLs \n' 
    for i in range(len(exist_ion_flag_list)): 
        if exist_ion_flag_list[i] == True and len(refine_ion_list[i][0]) >= 1: 
            add_info = list2string(refine_ion_list[i][0]) + ' |' + list2string(refine_ion_list[i][1]) + '\n'
            lines[i+1] = lines[i+1][:-1] + '\t' + add_info 
    sort_lines = sorted(lines[1:], key=lambda k: int(k.strip().split('\t')[5]), reverse=True) 
    idx = 0 
    new_sort_lines = [lines[0]]
    for line in sort_lines:
        new_line = str(idx) + '\t' + line.split('\t', 1)[1] 
        idx += 1
        new_sort_lines.append(new_line)

    with open(summary_path, 'w', encoding='utf-8') as f: 
        for line in new_sort_lines: 
            f.write(line)


def list2string(num_list): 
    str_num = ''
    for num in num_list: 
        str_num += str(num)+', ' 
    return str_num[:-2]


def ion_function_input(light_mass, mass_diff_dict, mass_diff): 
    heavy_mass = str(int(light_mass) + mass_diff) 
    t_light_mass = 'PFIND_DELTA_' + light_mass 
    t_heavy_mass = 'PFIND_DELTA_' + heavy_mass 
    modification_list = [t_light_mass, t_heavy_mass] 
    modification_dict = {} 
    modification_dict[t_light_mass] = mass_diff_dict[light_mass] 
    modification_dict[t_heavy_mass] = mass_diff_dict[heavy_mass] 
    return modification_list, modification_dict 
    


# 绘制结果热力图
def heat_map_plot(current_path, filter_mod, mod_static_dict, mod_number_dict, new_mod_number_dict): 
    y_stick = ["N-SIDE", "C-SIDE", "A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"] 
    y_dict = {}
    for i in range(len(y_stick)):
        y_dict[y_stick[i]] = i 
    
    x_stick = filter_mod
    mod_map = {}
    
    
    for i in range(len(x_stick)):
        freq_list = [0.0] * len(y_stick) 
        local_list = mod_static_dict[x_stick[i]].most_common() 
        Top1_site = local_list[0][0]
        residual_pro = 0.0 
        factor = (new_mod_number_dict[x_stick[i]] - mod_number_dict[x_stick[i]] + local_list[0][1]) * mod_number_dict[x_stick[i]] / new_mod_number_dict[x_stick[i]] /  local_list[0][1] 


        if Top1_site == 'N-SIDE' or Top1_site == 'C-SIDE': 
            parent_num = 0 
            for j in range(1, len(local_list)): 
                parent_num += local_list[j][1] 
        else: 
            parent_num = new_mod_number_dict[x_stick[i]] 
        
        for j in range(1, len(local_list)): 
            cur_position = local_list[j][0] 
            if cur_position not in y_stick: 
                continue 
            cur_freq = round(local_list[j][1]/parent_num, 3) 
            if cur_position == 'N-SIDE' or cur_position == 'C-SIDE': 
                cur_freq = cur_freq * parent_num * factor / mod_number_dict[x_stick[i]]
            freq_list[y_dict[cur_position]] = cur_freq 
            if cur_position == 'N-SIDE' or cur_position == 'C-SIDE': 
                continue 
            residual_pro += cur_freq 
        
        if Top1_site == 'N-SIDE' or Top1_site == 'C-SIDE': 
            Top1_pro = round(local_list[0][1] / parent_num, 3)
        else: 
            Top1_pro = round(1 - residual_pro, 3) 
        freq_list[y_dict[Top1_site]] = Top1_pro 

        mod_map[x_stick[i]] = freq_list 
    
    pd_mod_map = pd.DataFrame(mod_map, index=y_stick, columns=x_stick) 
    # print(pd_mod_map)
    # ax = sns.heatmap(pd_mod_map, vmin=0.0, vmax=1.0, cmap='YlGnBu', annot=True, annot_kws={"size":4})
    ax = sns.heatmap(pd_mod_map, vmin=0.0, vmax=1.0, cmap='YlGnBu')
    plt.ylabel('amino acid selectivity')
    plt.xlabel('modifications')
    png_path = os.path.join(current_path, 'heat_map.pdf') 
    # plt.show()
    plt.savefig(png_path, dpi=200, bbox_inches='tight')
    plt.close()


# 保留整数，其余信息进行合并 
def mass_refine(mod_static_dict, mod_number_dict, mod2pep, mass_diff_dict, parameter_dict, unimod_dict, pattern): 
    
    new_mod_static_dict, new_mod_number_dict, new_mod2pep, new_mass_diff_dict = {}, {}, {}, {}
    new_mass_list = [] 
    int_mass_list = []
    for mass in mod_static_dict: 
        if mass not in mass_diff_dict.keys():
            continue
        if mass_diff_dict[mass] < parameter_dict['min_mass_modification']: 
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
            t_mass = int(mod2pep[mass])
            new_mod2pep[simple_mass] = t_mass 
            new_mass_diff_dict[simple_mass] = mass_diff_dict[mass] 
    # rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank = label_determine(new_mod_number_dict, int_mass_list, int(parameter_dict['mass_diff_diff'])) 
    # min_num = modification_filter_frequency(new_mod_number_dict, parameter_dict) 
    
    rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, ppm_error_dict = accurate_label_determine(new_mod_number_dict, new_mass_diff_dict, parameter_dict, pattern)
    # print('new', mass_diff_pair_rank)
    unimod_list = unimod_match(unimod_dict, mass_diff_pair_rank, new_mass_diff_dict)
    # print('total: ', mass_diff_rank)
    return new_mod_static_dict, new_mod_number_dict, new_mod2pep, new_mass_diff_dict, rank_dict, label_dict, mass_diff_rank, mass_diff_pair_rank, unimod_list, ppm_error_dict   


# 删选低于频率低于filter_frequency的修饰
def modification_filter_frequency(mod_number_dict, parameter_dict): 
    total_sum = 0 
    for mod in mod_number_dict.keys():
        total_sum += mod_number_dict[mod] 
    if parameter_dict['filter_frequency'] <0 or parameter_dict['filter_frequency'] > 99:
        parameter_dict['filter_frequency'] = 0
    min_num = int(total_sum * parameter_dict['filter_frequency'] * 0.01) 
    return min_num 


# 精确质量法确定轻重标记
# parameter_dict['mass_of_diff_diff'] = 6.020132
# parameter_dict['mass_diff_diff_range'] = 100 
def accurate_label_determine(mod_number_dict, mass_diff_dict, parameter_dict, pattern): 
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
        # 最小数目限制，认为噪声 
        if mod_number_dict[mod] < 4:
            break 
        if mod in mass_diff_rank: 
            continue 
        
        # 只报出找到匹配的修饰对
        for j in range(i+1, len(freq_list)): 
            if rank_tuple[j][0] in mass_diff_rank:
                continue 
            pair_mod = rank_tuple[j][0] 
            #if mod_number_dict[pair_mod] < min_num: 
            #    continue 
            # 这里改成了绝对值限制 
            ppm_error = abs_ppm_calculate(mass_diff_dict[mod], mass_diff_dict[pair_mod], parameter_dict['mass_of_diff_diff'])
            # print(pattern)
            if ppm_error <= parameter_dict['mass_diff_diff_range'] or (pattern != 'blind' and ppm_error < 0.01): 
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


# 防止分母为0报错 
def ppm_calculate(a, b, mass_diff_diff): 
    return abs(abs(b-a)-mass_diff_diff)/(mass_diff_diff+0.000001)*1000000 

# 绝对值限制 
def abs_ppm_calculate(a, b, mass_diff_diff): 
    return abs(abs(b-a) - mass_diff_diff) 


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
            if abs(v - target_diff) < 0.1: 
                t_unimod += k + ';' 
        target_diff = cur_mod_mass - baseline_mass 
        for k, v in unimod_dict.items():
            if abs(v - target_diff) < 0.1: 
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


# 读取解释性的修饰 
def explain_dict_generate(current_path): 
    bin_path = os.path.join(current_path, 'bin')
    template_path = os.path.join(bin_path, 'template')
    explain_file = os.path.join(template_path, 'explain_modification.ini') 
    with open(explain_file, 'r', encoding='utf') as f: 
        lines = f.readlines() 
    explain_dict = {}
    for i in range(len(lines)): 
        if len(lines[i]) < 4:
            break 
        if 'name' in lines[i]: 
            mod_name = lines[i].split()[0].split('=')[1] 
            mod_mass = float(lines[i+1].split()[2]) 
            explain_dict[mod_name] = mod_mass 
    return explain_dict 
            


# 删选PSM低于指定阈值百分比的输出 
def summary_filter(current_path, parameter_dict, filter_mod, pattern, summary_path, prob_dict): 
    if parameter_dict['filter_frequency'] < 0.0:
        print('filter_frequency out of range!')
        return filter_mod
    if parameter_dict['filter_frequency'] > 100.0:
        print('filter_frequency out of range!') 
        return filter_mod 
    
    # summary_path = os.path.join(current_path, 'pChem.summary')
    new_filter_mod = []
    with open(summary_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    total_psm = 0 
    for line in lines[1:]: 
        if len(line) < 5:
            break 
        total_psm += int(line.split('\t')[5]) 
    print(total_psm)
    if pattern == 'blind':
        min_psm = parameter_dict['filter_frequency'] * total_psm * 0.01 
    else:
        # 限定式不需要过滤PSM
        min_psm = 0
    new_lines = []
    new_lines.append(lines[0]) 
    i = 1 
    print('Current pChem filter out ', str(int(parameter_dict['filter_frequency'])), '% PSM')
    for line in lines[1:]: 
        if int(line.split('\t')[5]) >= min_psm: 
            new_filter_mod.append(line.split('\t')[1][12:]) 
            t_idx = line.find('\t')
            new_lines.append(str(i) + line[t_idx:])
            i += 1 
    # print(new_lines)  
    
    with open(summary_path, 'w', encoding='utf-8') as f: 
        for line in new_lines: 
            f.write(line) 
    
    if pattern == 'blind' and parameter_dict['use_close_search'] == 'True': 
        # print('plot radar!')
        metric_evaluation(current_path, parameter_dict, summary_path, new_filter_mod, prob_dict)  

    return new_filter_mod



# 数据集级别指标评价 
# Identification efficiency： mod_psm /  all_psm, Modification uniformity：max_mod_psm / mod_psm, Position selectivity
def metric_evaluation(current_path, parameter_dict, summary_path, new_filter_mod, prob_dict): 
    # print(current_path) 
    blind_res = os.path.join(parameter_dict['output_path'], 'blind') 
    blind_res = os.path.join(blind_res, 'pFind.summary') 
    
    with open(blind_res, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    total_psm = int(lines[1].split(':')[1])
    
    # summary_path = os.path.join(current_path, 'pChem.summary') 
    with open(summary_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    # print(prob_dict)
    select_pos_summary = 0 
    psm_summary = 0 
    max_psm = 0  
    mod_name = ''
    if len(lines) < 2:
        print('please change the range of diff!') 
        return 
    
    for i in range(1, len(lines)): 
        line = lines[i].split('\t') 
        # light_mod = line[1]
        mod_name = line[1][12:]
        prob_list = prob_dict[mod_name]
        cur_psm = int(line[5]) 
        psm_summary += cur_psm 
        # line3 = line[3].split('|')
        if i == 1: 
            top1_pos = prob_list[-1][0]
            # print(top1_pos)
            max_psm = cur_psm 
        for pair in prob_list:
            if pair[0] == top1_pos:
                select_pos_summary += float(pair[1])*cur_psm 
        #cur_pos = line[3] 
        #if cur_pos == top1_pos: 
        #    select_pos_summary += float(line3[1])*cur_psm 
    
    #print(psm_summary/total_psm)
    #print(max_psm/psm_summary)
    #print(select_pos_summary/psm_summary)
    x = [psm_summary/total_psm*100.0, max_psm/psm_summary*100.0, select_pos_summary/psm_summary*100.0] 

    # x = update_pos_selectivity(x, parameter_dict, new_filter_mod) 

    metric_path = os.path.join(current_path, 'pChem.metric')
    metric2summary(metric_path, x) 
    # print('metric', x)
    radar_plot(x, current_path, mod_name)



def update_pos_selectivity(x, parameter_dict, new_filter_mod): 
    #print('new filter', new_filter_mod) 
    #print(position_dict) 
    # parameter_dict['side_position']
    psm_res = os.path.join(parameter_dict['output_path'], 'blind') 
    psm_res = os.path.join(psm_res, 'pFind-Filtered.spectra') 
    mass_diff_list = [] 
    for mod in new_filter_mod:
        mass_diff_list.append('PFIND_DELTA_' + mod)
        mass_diff_list.append('PFIND_DELTA_' + str(int(mod)+6))
    
    with open(psm_res, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    spectra_num = len(lines) 
    
    mod_position_dict = {} 
    for mod in mass_diff_list:
        mod_position_dict[mod] = [] 
    
    if parameter_dict['side_position'] == 'False': 
        side_flag = False
    else:
        side_flag = True
        mod_position_dict['N-SIDE'] = [] 
        mod_position_dict['C-SIDE'] = [] 

    total_psm_num = 0
    for i in range(1, spectra_num):
        if len(lines[i]) < 4:
            break 
        sequence = lines[i].split('\t')[5]
        mod_list = lines[i].split('\t')[10].split(';')[:-1]
        for mod in mod_list:
            pos, mod_name = mod.split(',') 
            mod_name = mod_name.split('.')[0]
            if mod_name in mass_diff_list:
                pos = int(pos) 
                total_psm_num += 1 
                if pos == 0 or pos == 1: 
                    mod_position_dict[mod_name].append(sequence[0]) 
                    if side_flag == True:
                        mod_position_dict[mod_name].append('N-SIDE')
                elif pos >= len(sequence):
                    mod_position_dict[mod_name].append(sequence[-1]) 
                    if side_flag == True:
                        mod_position_dict[mod_name].append('C-SIDE')
                else:
                    mod_position_dict[mod_name].append(sequence[pos-1])
    
    top_num = 0 
    max_pos = 'C'
    for i in range(len(mass_diff_list)):
        mod_name = mass_diff_list[i] 
        
        counter = Counter(mod_position_dict[mod_name]) 
        if i == 0: 
            max_pos = counter.most_common()[0][0]
        top_num += counter[max_pos]
        
    x[2] = float(top_num/ total_psm_num) * 100
    return x 





def update_identification_efficiency(current_path, parameter_dict): 
    close_res = os.path.join(parameter_dict['output_path'], 'close') 
    close_res = os.path.join(close_res, 'pFind.summary') 
    
    with open(close_res, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    total_psm = int(lines[1].split(':')[1])
    
    summary_path = os.path.join(current_path, 'pChem-close.summary') 
    with open(summary_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    psm_summary = 0 

    if len(lines) < 2:
        print('please change the range of diff!') 
        return 
    
    for i in range(1, len(lines)): 
        line = lines[i].split('\t') 
        cur_psm = int(line[5]) 
        psm_summary += cur_psm 
    
    metric_path = os.path.join(current_path, 'pChem.metric') 
    x = metric_file_read(metric_path) 
    x[0] = psm_summary/total_psm*100.0 
    print(x)
    radar_plot(x, current_path)



def metric_file_read(metric_path): 
    with open(metric_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    x = []
    for i in range(1,4): 
        x.append(float(lines[i].split('\t')[1][:-1])) 
    return x 



# 将指标结果写入metric文件 
def metric2summary(path, x): 
    # print(path)
    lines = []
    lines.append('Performance Metrics: \n')
    lines.append('Profiling efficiency \t' + str(round(x[0],2)) + '\n') 
    lines.append('PDM homogeneity \t' + str(round(x[1],2)) + '\n') 
    lines.append('Residue selectivity \t' + str(round(x[2],2)) + '\n') 
    with open(path, 'w', encoding='utf-8') as f:
        for line in lines:
            f.write(line) 




# 绘制雷达图 
def radar_plot(x, current_path, mod_name=None):
    # matplotlib.rcParams['font.family']="SimHei"
    radar_labels = np.array(['Profiling efficiency','PDM homogeneity','Residue selectivity'])
    data = np.array(x)
    # 是否需要显示修饰 
    if mod_name is not None: 
        data_labels =(mod_name,)

    angles = np.linspace(0, 2*np.pi, 3, endpoint=False)
    graph = plt.figure(facecolor = "white") 
    
    plt.subplot(111, polar = True)
    plt.plot(angles, data,'o-',linewidth=1, alpha=0.2)
    plt.fill(angles, data, alpha=0.25)
    plt.thetagrids(angles*180/np.pi, radar_labels)
    #plt.figtext(0.52, 0.95, '霍兰德人格分析', ha='center', size=20) 
    if mod_name is not None: 
        legend = plt.legend(data_labels, loc = (0.94, 0.80), labelspacing = 0.1)
        plt.setp(legend.get_texts(), fontsize='large')
    plt.grid(True)
    plt.ylim(0, 100)
    # plt.show()
    pdf_path = os.path.join(current_path, 'radar.pdf') 
    #plt.savefig(pdf_path) 
    #plt.show()
    pp = PdfPages(pdf_path) 
    pp.savefig(graph)
    pp.close()
    plt.close()


if __name__ == "__main__": 
    P = 'NAHSATTWSGQYVGGAEAR' 
    mass = 0.0 
    for a in P:
        mass += amino_acid_dict[a]
    mass += proton_mass + h2o_mass
    print(mass)