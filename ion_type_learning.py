# 中性丢失鉴定 

from utils import parameter_file_read 
import os 
from collections import Counter 
import numpy as np 
import matplotlib.pyplot as plt 
from mass_diff_correction import element_dict, amino_acid_dict, common_dict_create 


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
    print('read mgf data...')
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
    segment = int(len(filtered_res) / 10)
    i = 0 
    for line in filtered_res:
        if i % segment == 0:
            print('finished ', i / segment * 10,  'percetage')
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
        i += 1
        # break  
    print(total_ion_diff_counter.most_common()[:20])
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



# 离子类型学习 
def ion_type_determine(current_path, modification_list, modification_dict): 
    pchem_cfg_path = os.path.join(current_path, 'pChem.cfg')
    parameter_dict = parameter_file_read(pchem_cfg_path) 
    # print(parameter_dict)

    # 质谱数据读取 
    mass_spectra_dict = {} 
    for msms_path in parameter_dict['msms_path']:
        mgf_path = msms_path.split('=')[1].split('.')[0] + '.mgf'
        cur_mass_spectra_dict = mgf_read(mgf_path) 
        mass_spectra_dict.update(cur_mass_spectra_dict)
    print('The number of spectra: ', len(mass_spectra_dict.keys())) 

    # 读取盲搜得到的结果 
    blind_path = os.path.join(parameter_dict['output_path'], 'blind')
    blind_res_path = os.path.join(blind_path, 'pFind-Filtered.spectra')
    blind_res = blind_res_read(blind_res_path) 

    # 读取常见修饰的列表 
    common_modification_dict = common_dict_create(current_path)
    # print(common_modification_dict)

    # 筛选有效的PSM 
    for modification in modification_list:
        filtered_res = PSM_filter(blind_res, modification) 
        total_ion_diff_counter, total_ion_diff_list, total_weight_ion_diff_list = ion_type_compute(filtered_res, modification, modification_dict[modification], common_modification_dict, mass_spectra_dict) 
        # 画频率图 
        # freq_line_plot(total_ion_diff_counter) 
        
        for ion_mass, _ in total_ion_diff_counter.most_common()[:8]:
            accurate_ion_mass = accurate_ion_mass_computation(ion_mass, total_ion_diff_list) 
            weight_accurate_ion_mass = weight_accurate_ion_mass_computation(ion_mass, total_weight_ion_diff_list)
            print(accurate_ion_mass)
            print(weight_accurate_ion_mass)
            
        # print(filtered_res)


if __name__ == "__main__": 
    current_path = os.getcwd() 
    modification_list = ['PFIND_DELTA_252', 'PFIND_DELTA_258']
    modification_dict = {'PFIND_DELTA_252': 252.121858, 'PFIND_DELTA_258':258.141955}
    ion_type_determine(current_path, modification_list, modification_dict) 
    
