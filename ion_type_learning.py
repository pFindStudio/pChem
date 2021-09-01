# 中性丢失鉴定 

from utils import parameter_file_read 
import os 
from collections import Counter 
import numpy as np 
import matplotlib.pyplot as plt 
from parameter import element_dict, amino_acid_dict, common_dict_create 

h2o_mass = element_dict["H"]*2 + element_dict["O"]*1
proton_mass = element_dict['Pm']

# 存放谱图的类
class MassSpectrum:
    def __init__(self, charge, pepmass, peak_list):
        self.charge = charge 
        self.pepmass = pepmass
        self.peak_list = peak_list 


# 从mgf文件中读取谱图
def mgf_read(mgf_path):
    mass_spectra_dict = {}
    with open(mgf_path, 'r') as f:
        lines = f.readlines()
    # print('reading mgf data.')
    i = 0  
    while i < len(lines): 
        if 'BEGIN' in lines[i]: 
            i += 1
            spectrum_name = lines[i].split('=')[1].strip()
            i += 1
            spectrum_charge = int(lines[i].split('=')[1][0])
            i += 2 
            spectrum_pepmass = float(lines[i].split('=')[1])
            spectrum_peak_list = []
            while i < len(lines):
                i += 1 
                if 'END' in lines[i]:
                    break 
                spectrum_peak_list.append(lines[i].split()) 
            # print(spectrum_peak_list)
        spectrum = MassSpectrum(spectrum_charge, spectrum_pepmass, spectrum_peak_list)
        mass_spectra_dict[spectrum_name] = spectrum 
        i += 1 
        # break 
    # print('The number of spectra: ', len(mass_spectra_dict.keys()))
    return mass_spectra_dict


# 读取盲搜的结果 
def blind_res_read(blind_res_path):
    # print(blind_res_path) 
    with open(blind_res_path, 'r') as f:
        lines = f.readlines() 
    return lines[1:]


# 筛选出含有指定修饰的PSM
def PSM_filter(blind_res, modification): 
    filtered_res = []
    for line in blind_res: 
        if modification in line:
            filtered_res.append(line) 
    return filtered_res 


# 给定修饰的位置和质量，修改质量数组
def mass_vector_modify(pos, mass, mass_vector):
    if pos == 0:
        mass_vector[0] += mass 
    elif pos == len(mass_vector)+1:
        mass_vector[len(mass_vector)] += mass
    else:
        mass_vector[pos-1] += mass 
    return mass_vector 


# 生成质量数组，每个位置对应质量 
def mass_vector_generation(peptide_sequence, mod_list, modification, accurate_mass, common_modification_dict):
    mass_vector = []
    
    # 生成原始的质量数组
    for amino_acid in peptide_sequence:
        mass_vector.append(amino_acid_dict[amino_acid])
    # print(mass_vector)

    pos_mod_list = mod_list.split(';')[:-1] 
    for pos_mod in pos_mod_list: 
        pos, mod_name = pos_mod.split(',')
        if modification in mod_name:
            mod_pos = int(pos)
            mass_vector = mass_vector_modify(mod_pos, accurate_mass, mass_vector)
        else:
            mass_vector = mass_vector_modify(int(pos), common_modification_dict[mod_name], mass_vector)
    return mass_vector, mod_pos 


# 位置校准 
def position_correct(mod_pos, sequence_len):
    if mod_pos == 0:
        return 0
    elif mod_pos == sequence_len + 1:
        return sequence_len - 1
    else:
        return mod_pos - 1


# 生成b,y离子谱峰的数组 
def b_y_vector_generation(mass_vector, mod_pos):
    mod_peak_list = []
    sequence_len = len(mass_vector) 
    mass_sum_vector = [mass_vector[0]] 
    
    for i in range(1, sequence_len):
        mass_sum_vector.append(mass_vector[i] + mass_sum_vector[i-1])
    # print(mass_vector)
    # print(mass_sum_vector)

    mod_pos = position_correct(mod_pos, sequence_len) 

    # 生成b离子谱峰
    mod_peak_list = [mass + proton_mass for mass in mass_sum_vector[mod_pos:]]
    # 生成y离子谱峰 
    mod_peak_list += [mass_sum_vector[sequence_len-1] - mass + h2o_mass + proton_mass for mass in mass_sum_vector[:mod_pos]]

    return mod_peak_list 


# 统计质量偏差
def ion_diff_sum(mod_peak_list, peak_list):
    ion_diff_counter = Counter() 
    ion_diff_list = [] 
    weight_ion_diff_list = []
    for mod_peak in mod_peak_list: 
        # 质量偏差保留2位小数 
        ion_diff_fine = [round(mod_peak - float(peak[0]), 6) for peak in peak_list] 
        ion_diff_list += ion_diff_fine 
        weight_ion_diff_fine = [[round(mod_peak - float(peak[0]), 6), float(peak[1])] for peak in peak_list] 
        weight_ion_diff_list += weight_ion_diff_fine
        ion_diff = [round(mass, 2) for mass in ion_diff_fine]
        ion_diff_counter.update(ion_diff) 
    return ion_diff_counter, ion_diff_list, weight_ion_diff_list  



# 统计中性丢失的数目 
def ion_type_compute(filtered_res, modification, accurate_mass, common_modification_dict, mass_spectra_dict): 
    total_ion_diff_counter = Counter() 
    total_ion_diff_list = []
    total_weight_ion_diff_list = []
    #segment = int(len(filtered_res) / 10)
    #i = 0 
    for line in filtered_res:
        #if i % segment == 0:
        #    print('finished ', i / segment * 10,  'percetage')
        line_split = line.split('\t')
        spectrum_name, peptide_sequence, mod_list = line_split[0], line_split[5], line_split[10]
        # print(spectrum_name, peptide_sequence, mod_list)
        mass_vector, mod_pos = mass_vector_generation(peptide_sequence, mod_list, modification, accurate_mass, common_modification_dict)
        # print(mass_vector, mod_pos)

        # 生成b/y离子有关的谱峰 
        mod_peak_list = b_y_vector_generation(mass_vector, mod_pos)
        
        # print(mod_peak_list) 
        if spectrum_name in mass_spectra_dict.keys():
            peak_list = mass_spectra_dict[spectrum_name].peak_list
        else:
            continue
        ion_diff_counter, ion_diff_fine, weight_ion_diff = ion_diff_sum(mod_peak_list, peak_list)
        total_ion_diff_counter.update(ion_diff_counter) 
        total_ion_diff_list += ion_diff_fine 
        total_weight_ion_diff_list += weight_ion_diff 
        # i += 1
        # break  
    # print(total_ion_diff_counter.most_common()[:10])
    return total_ion_diff_counter, total_ion_diff_list, total_weight_ion_diff_list


# 绘制频率曲线图 
def freq_line_plot(total_ion_diff_counter): 
    list_len = total_ion_diff_counter.most_common()[0][1] 
    freq_list = [0] * (list_len + 1) 
    x = [i for i in range(0, list_len+1)]
    for _, v in total_ion_diff_counter.items(): 
        freq_list[v] += 1 
        #freq_list.append(v)
    # 画直方图太笼统 
    #data = np.array(freq_list)
    #plt.hist(data,bins=20)
    plt.plot(x[1:], freq_list[1:])
    plt.xlabel('occur times')
    plt.ylabel('frequency')
    plt.show()


# 精准质量计算 
def accurate_ion_mass_computation(coarse_mass, total_ion_diff_list): 
    mass_sum = 0.0 
    mass_num = 0 
    for ion_diff in total_ion_diff_list:
        if(abs(ion_diff - coarse_mass) < 0.05):
            mass_sum += ion_diff 
            mass_num += 1
    return mass_sum / mass_num 


# 谱峰加权的精准质量计算 
def weight_accurate_ion_mass_computation(coarse_mass, total_weight_ion_diff_list):
    mass_sum = 0.0 
    mass_num = 0.0 
    for ion_diff, peak in total_weight_ion_diff_list:
        if(abs(ion_diff - coarse_mass) < 0.01):
            mass_sum += ion_diff * peak 
            mass_num += peak 
    return mass_sum / mass_num 


# 对离群点进行检测和分析 
def freq_analysis(total_ion_diff_counter): 
    counter_list = []
    for k, v in total_ion_diff_counter.items(): 
        counter_list.append([k, v]) 
    counter_list = sorted(counter_list, key=lambda x: x[1], reverse=True)
    print(counter_list[:10])
    print(len(counter_list))
    counter_list = counter_cluster(counter_list) 
    counter_list = sorted(counter_list, key=lambda x: x[1], reverse=True)
    print(counter_list[:10])
    value_counter_list = [p[1] for p in counter_list]
    arr_mean = np.mean(value_counter_list) 
    # 方差 np.var
    arr_var = np.std(value_counter_list, ddof=1)
    print(arr_mean, arr_var) 


# 对近似点进行聚类 
def counter_cluster(counter_list):
    new_counter_list = [] 
    for pair in counter_list:
        if pair[1] > 1:
            flag = False 
            for cur_pair in new_counter_list: 
                if abs(pair[0]-cur_pair[0]) < 0.02: 
                    cur_pair[1] += pair[1]
                    flag = True
                    break 
            if flag == False: 
                new_counter_list.append(pair) 
        #else:
        #    new_counter_list.append(pair)
    return new_counter_list 


# 特征离子确定 (使用了全集，不正确)
def old_feature_peak_determine(mass_spectra_dict): 
    position_list = []
    for key in mass_spectra_dict:
        cur_position_list = [round(float(p[0]),2) for p in mass_spectra_dict[key].peak_list]
        position_list += cur_position_list
        
    position_counter = Counter(position_list) 
    print(position_counter.most_common()[:10])


def feature_peak_determine(mass_spectra_dict, filtered_res): 
    # 保留6位小数的
    fine_position_list = []
    coarse_position_list = [] 
    # 保存谱峰的精确质量 
    peak_dict = {} 
    for line in filtered_res: 
        line_split = line.split('\t')
        spectrum_name, peptide_sequence, mod_list = line_split[0], line_split[5], line_split[10]
        if spectrum_name in mass_spectra_dict.keys():
            peak_list = mass_spectra_dict[spectrum_name].peak_list
        else:
            continue 
        fine_position_list += [float(p[0]) for p in peak_list] 
        coarse_position_list += [round(float(p[0]),2) for p in peak_list]
    position_counter = Counter(coarse_position_list) 
    filtered_position_list = [p[0] for p in position_counter.most_common()[:300]] 
    for position in filtered_position_list: 
        peak_dict[position] = [] 
    for position in fine_position_list: 
        peak = round(position, 2)
        if peak in peak_dict.keys(): 
            peak_dict[peak].append(position) 
    for key in peak_dict.keys(): 
        peak_dict[key] = np.mean(peak_dict[key])
    #print(peak_dict)

    return filtered_position_list, peak_dict 


def ppm_calculate(a, b, mass_diff_diff): 
    return abs(abs(b-a)-mass_diff_diff)/(mass_diff_diff+0.000001)*1000000 


def feature_pair_find(position_list, peak_dict, mass_diff): 
    # print(position_list) 
    coarse_mass_diff = round(mass_diff, 2)
    pair_list = [] 
    for light_mass in position_list[0]:
        for heavy_mass in position_list[1]:
            if round(heavy_mass - light_mass, 2) == coarse_mass_diff: 
                ppm_error = ppm_calculate(peak_dict[0][light_mass], peak_dict[1][heavy_mass], mass_diff)
                pair_list.append([light_mass, heavy_mass, int(ppm_error)]) 
    print(pair_list)
    return pair_list


# 离子类型学习 
def ion_type_determine(current_path, modification_list, modification_dict, parameter_dict): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    # parameter_dict = parameter_file_read(pchem_cfg_path) 
    # print(parameter_dict)

    # 质谱数据读取 
    mass_spectra_dict = {} 
    for msms_path in parameter_dict['msms_path']:
        mgf_path = msms_path.split('=')[1].split('.')[0] + '.mgf'
        cur_mass_spectra_dict = mgf_read(mgf_path) 
        mass_spectra_dict.update(cur_mass_spectra_dict)
    # print('The number of spectra: ', len(mass_spectra_dict.keys())) 
    #feature_peak_determine(mass_spectra_dict) 
    

    # 读取盲搜得到的结果 
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_res_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    blind_res = blind_res_read(blind_res_path) 

    # 读取常见修饰的列表 
    common_modification_dict = common_dict_create(current_path)
    # print(common_modification_dict)
    position_list = [] 
    peak_dict = [] 

    # 修饰->中性丢失 
    mod2ion = {}
    exist_ion_flag = True 
    ion_list = []

    # 筛选有效的PSM 
    for modification in modification_list:
        # pfind-filtered line 
        mod = modification.split('_')[2] 
        int_mod = int(mod)
        mod2ion[mod] = [] 
        t_ion_list = []
        filtered_res = PSM_filter(blind_res, modification) 

        # 确定报告离子 
        t_position_list, t_peak_dict = feature_peak_determine(mass_spectra_dict, filtered_res)
        peak_dict.append(t_peak_dict)
        position_list.append(t_position_list) 

        
        # total_ion_diff_list: 所有质量差组成的列表 
        total_ion_diff_counter, total_ion_diff_list, total_weight_ion_diff_list = ion_type_compute(filtered_res, modification, modification_dict[modification], common_modification_dict, mass_spectra_dict) 
        # 画频率图 
        #freq_line_plot(total_ion_diff_counter) 
        # freq_analysis(total_ion_diff_counter)

        # 判断是否存在中性丢失 
        if int(total_ion_diff_counter.most_common()[0][0]) == 0:
            exist_ion_flag = False 

        repeat_list = []
        for ion_mass, _ in total_ion_diff_counter.most_common()[:10]:
            accurate_ion_mass = accurate_ion_mass_computation(ion_mass, total_ion_diff_list) 
            # weight_accurate_ion_mass = weight_accurate_ion_mass_computation(ion_mass, total_weight_ion_diff_list)
            # print('average: ', accurate_ion_mass)
            # print('weight: ', weight_accurate_ion_mass) 
            int_ion_mass = int(accurate_ion_mass) 
            if int_ion_mass not in repeat_list and int_mod != int_ion_mass: 
                mod2ion[mod].append(accurate_ion_mass) 
                t_ion_list.append(accurate_ion_mass)
                repeat_list.append(int_ion_mass)
        ion_list.append(t_ion_list)
        # print(filtered_res)
    # 特征离子发现
    #print('Feature ion results')
    #pair_list = feature_pair_find(position_list, peak_dict, parameter_dict['mass_of_diff_diff']) 
    return mod2ion, ion_list, exist_ion_flag


# 利用轻重标记筛选中性丢失 
def ion_filter(ion_list, mass_diff): 
    
    refine_ion_list = [[], []] 
    for light_mass in ion_list[0]:
        for heavy_mass in ion_list[1]: 
            diff = int(heavy_mass - light_mass) 
            if diff == mass_diff: 
                refine_ion_list[0].append(round(light_mass, 6)) 
                refine_ion_list[1].append(round(heavy_mass, 6)) 
            if len(refine_ion_list[0]) >= 3:
                return refine_ion_list
    return refine_ion_list 



if __name__ == "__main__": 
    current_path = os.getcwd() 
    # 需要输入待确认的未知
    modification_list = ['PFIND_DELTA_333', 'PFIND_DELTA_339']
    modification_dict = {'PFIND_DELTA_333': 333.168167, 'PFIND_DELTA_339':339.187878}
    ion_type_determine(current_path, modification_list, modification_dict) 
    
