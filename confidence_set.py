from collections import Counter
import os 



# 将鉴定得到的未知修饰质量写入盲搜/限定式结果文件
def accurate_mass_for_result_file(file_path, mass_diff_dict, parameter_dict, mode='blind'): 
    with open(file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    # print(lines) 

    revise_write = False 
    new_lines = []
    for line in lines: 
        if 'Frequency' in line: 
            revise_write = True 
            new_lines.append(line[:-1] + '\t' + 'Accurate Mass' + '\n')
            continue 
        if revise_write == True: 
            query = line.split()[0] 
            if mode == 'close':
                query = query.split('.')[0] 
            if mode == 'close' and parameter_dict['isotope_labeling'] == 'False':
                query = query.split('_')[2]
            if query in mass_diff_dict.keys(): 
                new_lines.append(line_rewrite(line, mass_diff_dict[query]))
            else:
                new_lines.append(line)
            continue
        if revise_write == True and '----------' in line: 
            revise_write = False 
        new_lines.append(line) 
    #for line in new_lines: 
    #    print(line)

    with open(file_path, 'w', encoding='utf') as f:
        for line in new_lines:
            f.write(line)



def line_rewrite(line, mass): 
    return line[:-1] + '\t' + str(mass) + '\n'



# 计算平均剩余残差平方 
import numpy as np
def residual_compute(target_mass, mass_list): 
    # print(np.mean(mass_list), np.std(mass_list))
    sample_mean = np.mean(mass_list)
    predict = 0.0
    gt = 0.0
    for mass in mass_list: 
        predict += (mass - target_mass) * (mass - target_mass) 
        gt += (mass - sample_mean) * (mass - sample_mean) 
    if predict == 0.0:
        r_value = 0.05
    else:
        r_value = float(gt / predict) 
    r_value = max(r_value, 1.0 - r_value)  
    if r_value < 0.8: 
        r_value = 0.8 + r_value / 10.0 
    return r_value


# 调用sklearn计算r方
'''
from sklearn.metrics import r2_score
def r2_score_compute(target_mass_list, mass_list): 
    return r2_score(target_mass_list, mass_list, multioutput='uniform_average')
'''

# 计算位点分布的可信度 
from scipy.stats import binom_test 
from collections import Counter 
def position_test(position_counter, prior_distribution): 
    total_num = 0
    for v in position_counter.values():
        total_num += v 
    #n = len(position_counter)
    for k, v in position_counter.most_common(): 
        #p = binom_test(v, total_num, p=1/n) 
        p = binom_test(v, total_num, p=prior_distribution[k]) 
        print(k, v, p) 
        # n -= 1
        total_num -= v 
        



# 统计位点分布的先验频率
def prior_distribution_compute(file_path, mod): 

    with open(file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    prior_distribution = {'N-SIDE':0, 'C-SIDE':0}
    num = 0 
    psm = 0
    for line in lines: 
        if mod in line: 
            psm += 1
            sequence = line.split('\t')[5] 
            for aa in sequence: 
                if aa in prior_distribution.keys(): 
                    prior_distribution[aa] += 1 
                else:
                    prior_distribution[aa] = 1
                num += 1  
            prior_distribution['N-SIDE'] += 1
            prior_distribution['C-SIDE'] += 1 
    #print(prior_distribution) 
    if num == 0:
        num += 1
    for key in prior_distribution.keys(): 
        prior_distribution[key] = float(prior_distribution[key] / num) 
    #print(num)
    #print('psm: ', psm)
    #print(prior_distribution)
    return prior_distribution




# N-SIDE的统计方式需要考虑
def total_trail_compute(position_list):
    total_num = 0 
    for pair in position_list:
        #if 'N-CIDE' == pair[0]:
        #    continue 
        total_num += pair[1] 
    return total_num 



# 对每个修饰的所有报告位点计算置信度
# 计算p-value需要信息： mod name、position_list、parameter_dict、 pattern-> 去不同搜索结果里面统计先验频率 
def p_value_for_mod(mod_name, position_list, parameter_dict, pattern): 
    
    # 计算先验频率
    # source_path = os.path.join(parameter_dict['output_path'], 'source') 
    pattern_path = os.path.join(parameter_dict['output_path'], pattern) 
    pfind_summary_path = os.path.join(pattern_path, 'pFind-Filtered.spectra')

    prior_distribution = prior_distribution_compute(pfind_summary_path, 'PFIND_DELTA_' + mod_name) 
    total_trail = total_trail_compute(position_list) 

    p_value_dict = {} 
    for pair in position_list: 
        if pair[0] in prior_distribution.keys():
            p_value_dict[pair[0]] = format(binom_test(pair[1], total_trail, p=prior_distribution[pair[0]], alternative='greater'), '.4f')
        else:
            p_value_dict[pair[0]] = 1.0
    return p_value_dict 
    




if __name__ == "__main__": 

    blind_mass_dict = {'PFIND_DELTA_387.18': 387.175225, 'PFIND_DELTA_393.20': 393.196465, 
    'PFIND_DELTA_325.17': 325.166226, 'PFIND_DELTA_471.24': 471.234108, 'PFIND_DELTA_477.26': 477.257341, 'PFIND_DELTA_371.18': 371.178761, 'PFIND_DELTA_369.17': 369.165961, 'PFIND_DELTA_409.20': 409.19209, 'PFIND_DELTA_375.19': 375.187079, 'PFIND_DELTA_261.11': 261.105591, 
    'PFIND_DELTA_331.18': 331.182578, 'PFIND_DELTA_255.09': 255.086557, 'PFIND_DELTA_403.18': 403.171312, 'PFIND_DELTA_305.11': 305.114061, 'PFIND_DELTA_323.15': 323.146595, 
    'PFIND_DELTA_329.16': 329.157907, 'PFIND_DELTA_330.16': 330.15341, 'PFIND_DELTA_377.20': 377.198503, 'PFIND_DELTA_445.20': 445.197768, 'PFIND_DELTA_336.18': 336.172428, 'PFIND_DELTA_401.19': 401.188336, 'PFIND_DELTA_394.18': 394.178088, 'PFIND_DELTA_395.17': 395.168935, 
    'PFIND_DELTA_399.22': 399.221779, 'PFIND_DELTA_451.23': 0.0, 'PFIND_DELTA_425.27': 425.261016, 'PFIND_DELTA_419.17': 419.167726, 'PFIND_DELTA_388.14': 388.139371, 
    'PFIND_DELTA_389.16': 389.151762, 'PFIND_DELTA_299.16': 299.154141}

    close_mass_dict =  {'PFIND_DELTA_387': 387.175356, 'PFIND_DELTA_393': 393.195355, 
    'PFIND_DELTA_325': 325.16519, 'PFIND_DELTA_331': 331.184624, 'PFIND_DELTA_471': 471.234131, 
    'PFIND_DELTA_477': 477.25614, 'PFIND_DELTA_371': 371.180052, 'PFIND_DELTA_377': 377.200902} 

    blind_summary_path = 'results/source/blind/pFind.summary' 
    close_summary_path = 'results/source/close/pFind.summary' 

    # accurate_mass_for_result_file(blind_summary_path, blind_mass_dict) 
    # accurate_mass_for_result_file(close_summary_path, close_mass_dict)
    
    IPM_mod_static_dict = {'PFIND_DELTA_252.12': Counter({'C': 2783, 'N-SIDE': 304, 'A': 91, 'P': 7, 'G': 6, 'L': 5, 'S': 5, 'M': 4, 'N': 3, 'Q': 3, 'T': 3, 'H': 2, 'C-SIDE': 1, 'F': 1, 'I': 1, 'D': 1, 'E': 1}), 
                        'PFIND_DELTA_258.14': Counter({'C': 2150, 'N-SIDE': 243, 'A': 72, 'P': 8, 'G': 7, 'C-SIDE': 5, 'V': 3, 'S': 3, 'R': 3, 'W': 2, 'T': 2, 'Q': 2, 'N': 2, 'F': 2, 'K': 1, 'H': 1, 'E': 1, 'M': 1, 'D': 1}),} 
    
    ALK_mod_static_dict = {'PFIND_DELTA_387.18': Counter({'C': 87, 'N-SIDE': 5, 'A': 3, 'E': 2, 'T': 1, 'G': 1, 'D': 1, 'F': 1, 'K': 1, 'N': 1, 'M': 1}), 
                        'PFIND_DELTA_393.20': Counter({'C': 48, 'N-SIDE': 5, 'W': 2, 'F': 1, 'V': 1, 'G': 1})}
    new_IPM_mod_static_dict = {'PFIND_DELTA_252.12': Counter({'C': 2756, 'N-SIDE': 167, 'A': 12, 'P': 1}), 
                        'PFIND_DELTA_258.14': Counter({'C': 2150, 'N-SIDE': 243, 'A': 72, 'P': 8, 'G': 7, 'C-SIDE': 5, 'V': 3, 'S': 3, 'R': 3, 'W': 2, 'T': 2, 'Q': 2, 'N': 2, 'F': 2, 'K': 1, 'H': 1, 'E': 1, 'M': 1, 'D': 1}),} 
    Diazo_static_dict = {'PFIND_DELTA_267.16': Counter({'E': 101, 'Y': 49, 'C': 46, 'N-SIDE': 30, 'A': 10, 'D': 8, 'H': 7, 'S': 4, 'G': 4, 'T': 3, 'N': 3, 'V': 3, 'M': 2, 'I': 2, 'L': 1, 'F': 1, 'Q': 1, 'P': 1}),}
    #for key in ALK_mod_static_dict.keys(): 
    #    print(key)
    #    position_test(ALK_mod_static_dict[key]) 
    
    blind_psm_res_path = 'results/diazo_new/source/blind/pFind-Filtered.spectra' 
    mod = 'PFIND_DELTA_267.16' 
    prior_distribution = prior_distribution_compute(blind_psm_res_path, mod) 
    position_test(Diazo_static_dict[mod], prior_distribution) 

    mass_list = [387.175302, 387.174946, 387.177513, 387.176867, 387.175535, 387.174119, 387.177558, 387.17416, 387.175803, 387.179019, 387.172765, 387.182742, 387.172077, 387.17845, 387.176522, 387.174005, 387.177563, 387.179487, 387.173134, 387.177723, 387.176265, 387.172865, 387.182995, 387.172778, 387.17891, 387.180439, 387.178613, 387.172788, 387.173217, 387.17622, 387.175217, 387.176924, 387.175817, 387.180149, 387.176449, 387.17569, 387.181579, 387.177501, 387.173593, 387.175542, 387.172815, 387.17932, 387.181172, 387.178474, 387.179223, 387.176096, 387.176307, 387.175256, 387.17769, 387.176268, 387.178428, 387.17709, 387.178226, 387.178464, 387.174615, 387.173003, 387.178307, 387.175803, 387.17508, 387.174872, 387.178746, 387.174259, 387.171746, 387.172536, 387.174429, 387.181688, 387.177593, 387.179464, 387.179066, 387.175107, 387.172969, 387.173353, 387.17193, 387.181708, 387.170776, 387.170155, 387.175783, 387.172327, 387.174484, 387.176246, 387.176897, 387.176668, 387.175057, 387.175174, 387.177178, 387.172155, 387.173713, 387.183545, 387.174813, 387.177579, 387.188771, 387.173412, 387.172199, 387.179957, 387.180594, 387.177088, 387.171014, 387.174384, 387.172217, 387.175892, 387.176847, 387.176071, 387.175387, 387.175234]
    target_mass = 387.175385 
    print(residual_compute(target_mass, mass_list))





