import numpy as np 
import os 
from scipy.stats import norm 
import matplotlib.pyplot as plt 
from collections import Counter 

from parameter import element_dict, amino_acid_dict, common_dict_create, h2o_mass, proton_mass 


def mass_diff_list_compute(lines, mass, common_dict=None):
    mass_diff_list = [] 
    no_modi_num = 0 
    for line in lines: 
        if 'PFIND_DELTA' in line: 
            no_modi_num += 1
        if mass not in line:
            continue 
        line = line.split('\t') 
        mod_list = line[10].split(';')[:-1] 
        #if len(mod_list) > 1:
        #    continue
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
        
        if mod_mass > 252.110 and mod_mass < 252.130:  
            mass_diff_list.append(round(mod_mass,4))
        
    print(len(lines), no_modi_num)
    return mass_diff_list 


def no_mod_diff_list_compute(lines, common_dict=None):
    no_mod_mass_diff_list = [] 
    for line in lines: 
        if 'PFIND_DELTA' in line: 
            continue

        line = line.split('\t') 
        mod_list = line[10].split(';')[:-1] 
        if len(mod_list) >= 1:
            continue
        parent_mass = float(line[2])
        sequence = line[5]
        amino_mass = 0.0
        for a in sequence: 
            if a in amino_acid_dict.keys():
                amino_mass += amino_acid_dict[a]
        mod_mass = parent_mass - amino_mass - proton_mass - h2o_mass
        
        if len(mod_list) >= 1:
            for mod in mod_list:
                
                mod = mod.split(',')[1]
                mod_mass -= common_dict[mod]
        
        no_mod_mass_diff_list.append(mod_mass)
    
    print(np.mean(no_mod_mass_diff_list))


def gaussion_fit(mass_diff_list): 
    return norm.fit(mass_diff_list)


def hist_plot(data): 
    plt.hist(data, bins=100)
    plt.show()


def accurate_mass_fit(blind_path, mass_list): 
    current_path = os.getcwd() 
    common_dict = common_dict_create(current_path)
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    # origin_lines = lines 
    lines = lines[1:] 
    for mass in mass_list: 
        mass_diff_list = mass_diff_list_compute(lines, mass, common_dict) 
        hist_plot(mass_diff_list) 
        
        mean, _ = gaussion_fit(mass_diff_list)
        print(mean)
    no_mod_diff_list_compute(lines, common_dict)


def mass_read(blind_path, target_mass): 
    mass_diff_diff = 6.020132
    with open(blind_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    ori_data = []  
    for line in lines: 
        if 'PFIND_DELTA' in line: 
            mass = float(line.split('\t')[-1])
            # if mass >= target_mass - 0.01 and mass <= target_mass + 0.01: 
            ori_data.append(mass) 
    data = []
    for mass in ori_data: 
        if mass >= target_mass - 0.01 and mass <= target_mass + 0.01: 
            data.append(mass)
        #if index % 2 == 0:
        #    if mass - mass_diff_diff >= target_mass - 0.01 and mass - mass_diff_diff <= target_mass + 0.01: 
        #        data.append(mass- mass_diff_diff) 
        #if index % 2 == 1:
        #    if mass + mass_diff_diff >= target_mass - 0.01 and mass + mass_diff_diff <= target_mass + 0.01: 
         #       data.append(mass + mass_diff_diff) 
    inc_number = len(data)
    '''
    r = [target_mass - 0.01, target_mass + 0.01] 
    mu0, sigma0 = -1, -1 
    times = 1 
    while True:
        a = np.array(data)
        mu, sigma = np.mean(a), np.std(a)
        if mu0 == mu or times == 3:
            break
        mu0, sigma0 = mu, sigma
        min_dist = min(abs(mu - r[0]), abs(r[1] - mu))
        new_data = []
        for v in ori_data:
            if v >= mu - 3 * sigma and v <= mu + 3 * sigma and abs(v - mu) <= min_dist:
                new_data.append(v)
        data = new_data 
        times += 1
    return np.mean(data), len(ori_data)



    '''
    mu0, sigma0 = np.mean(data), np.std(data)
    times = 1 

    while True: 
        if len(data) <= 100:
            break
        a = [] 
        for mass in ori_data:
            if mass >= mu0 - 0.01 and mass <= mu0 + 0.01: 
                a.append(mass) 
            #if index % 2 == 0:
            #    if mass - mass_diff_diff >= mu0 - 0.01 and mass - mass_diff_diff <= mu0 + 0.01: 
            #        data.append(mass- mass_diff_diff) 
            #if index % 2 == 1:
            #    if mass + mass_diff_diff >= mu0 - 0.01 and mass + mass_diff_diff <= mu0 + 0.01: 
            #        data.append(mass + mass_diff_diff) 
        mu, sigma = np.mean(a), np.std(a) 
        data = a 
        if mu == mu0 or times == 3: 
            break 
        mu0, sigma0 = mu, sigma 
        times += 1 
    return mu0, inc_number
    

    
    

def ppm_calculate(a, target): 
    return abs(a - target)/(target+0.000001)*1000000 


if __name__ == "__main__":  
    
    close_path = 'D:/pChem/pChem_new/0.05Da/QE_Plus_YangJing_WYne_O_TCP_50per_20180823/blind/pFind-Filtered.spectra' 
    blind_path = 'D:/pChem/pChem_new/results/ALK-iter/blind/pFind-Filtered.spectra' 
    # mass_list = ['PFIND_DELTA_252.1'] 
    # accurate_mass_fit(blind_path, mass_list)  
    mass_list = [291.12, 297.14]
    dream_list = [291.12191, 297.142042]
    new_list = []
    for target_mass in mass_list: 
        new_list.append(mass_read(close_path, target_mass)[0]) 

    for i in range(len(mass_list)): 
        print(new_list[i])
        print(ppm_calculate(new_list[i], dream_list[i]))


    '''
    print(ppm_calculate(252.123487, 252.1222))
    print(ppm_calculate(258.143589, 258.142332))
    
    print(ppm_calculate(279.157048,279.1583))
    print(ppm_calculate(285.177051,285.178432))
    '''

