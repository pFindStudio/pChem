import os 
from utils import parameter_file_read 

def ppm_calculate(a, b, mass_diff_diff): 
    return abs(abs(b-a)-mass_diff_diff)/(mass_diff_diff+0.000001)*1000000  


def mgf_path_find(parameter_dict): 
    pparse_path = os.path.join(parameter_dict['output_path'], 'pparse') 
    file_list = os.listdir(pparse_path) 
    mgf_list = []
    for file_name in file_list: 
        if '.mgf' in file_name: 
            mgf_list.append(os.path.join(pparse_path, file_name)) 
    return mgf_list  
    


def toy_data_generate(parameter_dict): 

    # 从pfind-spectra
    blind_res_path = os.path.join(parameter_dict['output_path'], 'blind') 
    blind_file_path = os.path.join(blind_res_path, 'pFind-Filtered.spectra')
    
    with open(blind_file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    
    print('total psm number: ', len(lines) - 1) 
    
    
    id2mass = {} 
    mass_list = []
    for line in lines[1:]: 
        line = line.split('\t') 
        id, mass = line[0], float(line[2])
        # print(id, mass)
        id2mass[id] = mass 
        mass_list.append(mass) 
    
    filtered_list = [] 
    for id in id2mass.keys(): 
        mass = id2mass[id] 
        flag = False 
        for p_mass in mass_list: 
            if ppm_calculate(mass, p_mass, parameter_dict['mass_of_diff_diff']) < 50: 
                flag = True 
                break
        
        if flag == True: 
            filtered_list.append(id) 
    
    mgf_path_list = mgf_path_find(parameter_dict) 
    print(mgf_path_list) 

    refined_mgf_line = []
    for mgf_path in mgf_path_list: 
        with open(mgf_path, 'r', encoding='utf-8') as f:
            lines = f.readlines() 
        num_line = len(lines)
        i = 0 
        while i < num_line: 
            if 'BEGIN IONS' in lines[i]: 
                left_idx = i 
                spectrum_name = lines[i+1].strip().split('=')[1] 
                while 'END IONS' not in lines[i]:
                    i += 1 
                if spectrum_name in filtered_list: 
                    refined_mgf_line += lines[left_idx:i+1] 
            i += 1 
    print(len(filtered_list))
    with open('toy_data.mgf', 'w', encoding='utf-8') as f: 
        for line in refined_mgf_line: 
            f.write(line)




if __name__ == "__main__":  
    cur_path = os.getcwd() 
    cfg_path = os.path.join(cur_path, 'pChem.cfg') 
    parameter_dict = parameter_file_read(cfg_path)
    print(parameter_dict) 
    toy_data_generate(parameter_dict)


