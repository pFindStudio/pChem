# incorporate psite for amino acid position selection filter
import os
from collections import Counter 
from utils import parameter_file_read, parameter_modify 
# D:\pSite\pPredictAA\x64\Release


def psite_file_generate(file_path, target_mod, current_path): 
    with open(file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    
    idx = 1 
    new_lines = []
    spect2pos = {}
    for line in lines[1:]: 
        if target_mod not in line:
            continue 
        line = line.split() 
        spectrum, sequence, mod = line[0], line[5], line[10]  
        mod = mod.split(',')
        if len(mod) > 2:
            continue 
        pos = int(mod[0]) 
        if pos == 0: 
            sequence = 'm' + sequence 
            spect2pos[spectrum] = [sequence[0], 'N-SIDE']
        elif pos >= len(sequence):
            sequence = sequence + 'm' 
            spect2pos[spectrum] = [sequence[len(sequence)-1], 'C-SIDE']
        else:
            sequence = sequence[:pos] + 'm' + sequence[pos:] 
            spect2pos[spectrum] = [sequence[pos - 1]]

        new_lines.append('S' + str(idx) + '\t' + spectrum + '\n') 
        new_lines.append('P1\t' + sequence + '\t0\n') 
        new_lines.append('\n') 
        idx += 1 

    with open(os.path.join(current_path, 'psite.txt'), 'w', encoding='utf-8') as f: 
        for line in new_lines: 
            f.write(line) 
    # return spect2pos 


def spectra_name_list_generate(summary_file_path, mod): 
    spectra_name_list = [] 
    with open(summary_file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    mod = 'PFIND_DELTA_' + mod 
    for line in lines: 
        if mod not in line: 
            continue 
        else:
            spectra_name_list.append(line.split('\t')[0]) 
    return spectra_name_list



# 读取psite的结果文件
def psite_result_read(file_path, blind_summary_file_path, mod): 
    with open(file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    spec2score = {}
    #i = 0
    #total = 0
    spectra_name_list = spectra_name_list_generate(blind_summary_file_path, mod)
    for line in lines: 
        line = line.split('\t') 
        spectrum_name = line[0]
        if spectrum_name not in spectra_name_list:
            continue
        score = float(line[3])
        spec = line[0] 
        # 用于测试不同的psite阈值 
        #if score < 5.0:
        #    continue
        #print(score, spec2pos[spec]) 
        #if spec2pos[spec][0] == 'C':
        #    i += 1  
        #total += 1 
    #print(float(i/total)) 
        spec2score[spec] = score 
    
    return spec2score



def psite_cfg_write(psite_template_path, current_path, source_path): 
    with open(psite_template_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 

    pparse_path = os.path.join(source_path, 'pParse') 
    mgf_name = '.mgf'
    for name in os.listdir(pparse_path): 
        if '.mgf' in name: 
            mgf_name = name  
            break 
    
    mgf_path = os.path.join(pparse_path, mgf_name) 
    psite_input_path = os.path.join(current_path, 'psite.txt')  

    new_lines = [] 
    for line in lines: 
        if 'mgfPath' in line: 
            line = parameter_modify(line, mgf_path) 
        if 'resultPath' in line: 
            line = parameter_modify(line, psite_input_path)
        new_lines.append(line) 
    

    with open(psite_template_path, 'w', encoding='utf-8') as f: 
        for line in new_lines: 
            f.write(line)



# 统计新的位点频率 拷贝自函数mass_static
def position_static(mod, spec_name_list, blind_summary_file_path, side_position='True'): 
    mod_position_list = [] 
    if side_position == 'True': 
        side_flag = True 
    else: 
        side_flag = False 
    with open(blind_summary_file_path, 'r', encoding='utf-8') as f:
        lines = f.readlines() 
    spectra_num = len(lines)


    for i in range(1, spectra_num):
        if len(lines[i]) < 4:
            break
        spec_name = lines[i].split('\t')[0] 
        if spec_name not in spec_name_list: 
            continue 
        
        sequence = lines[i].split('\t')[5]
        mod_list = lines[i].split('\t')[10].split(';')[:-1]
        for mod_line in mod_list: 
            pos, mod_name = mod_line.split(',')
            if mod in mod_name: 
                pos = int(pos)
                if pos == 0 or pos == 1:
                    mod_position_list.append(sequence[0]) 
                    if side_flag == True:
                        mod_position_list.append('N-SIDE')
                elif pos >= len(sequence):
                    mod_position_list.append(sequence[-1]) 
                    if side_flag == True:
                        mod_position_list.append('C-SIDE')
                else:
                    mod_position_list.append(sequence[pos-1]) 
    
    return Counter(mod_position_list)


def local_list_combine(local_list, position_counter): 
    new_local_list = [] 
    cut_flag = False 
    cut_num = int(local_list[0][1]/3) 
    for pos, num in local_list: 
        if position_counter[pos] <= 1 and len(new_local_list) > 1: 
            break 
        if num < cut_num: 
            cut_flag = True
        if cut_flag == False: 
            new_local_list.append([pos, max(num, position_counter[pos])]) 
        else: 
            if position_counter[pos] == 0: 
                new_local_list.append([pos, num]) 
            else:
                new_local_list.append([pos, min(num, position_counter[pos])]) 
    return new_local_list


def psite_run(parameter_dict, current_path, mod, pattern='blind', local_list=None): 

    # 1. 生成psite运行的参数和格式文件
    # parameter_dict = parameter_file_read(cfg_path) 

    #if parameter_dict['report_statistical'] == 'False': 
    #    print('No statistical information will be reported!')
    #    return 
    # print(parameter_dict) 
    # source_path = os.path.join(parameter_dict['output_path'], 'source')
    source_path = parameter_dict['output_path']
    blind_summary_file_path = os.path.join(source_path, pattern) 
    blind_summary_file_path = os.path.join(blind_summary_file_path, 'pFind-Filtered.spectra') 

    bin_path = os.path.join(current_path, 'bin') 
    psite_path = os.path.join(bin_path, 'pSite') 
    psite_template_path = os.path.join(psite_path, 'template') 
    psite_template_path = os.path.join(psite_template_path, 'param_pSite.txt') 
    
    # 产生输入文件 只有第一次调用才会运行 
    if parameter_dict['psite_run'] == 'True': 
        psite_file_generate(blind_summary_file_path, 'PFIND_DELTA_', current_path) 
        
        # 2. 生成psite参数文件 
        psite_cfg_write(psite_template_path, current_path, source_path)

        # 3. 运行pSite输出打分结果 
        # psite_exe_path = os.path.join(psite_path, 'pPredictAA.exe') 
        # 使用mingw64编译的psite不会报告UAC错误
        psite_exe_path = os.path.join(psite_path, 'a.exe') 
        cmd = psite_exe_path + ' ' + psite_template_path 
        os.chdir(psite_path)
        receive = os.system(cmd) 
        print(receive) 
        os.chdir(current_path) 

        parameter_dict['psite_run'] = 'False' 

    # 4. 读取结果文件，卡值后返回新的位点分布 
    psite_res_path = os.path.join(psite_path, 'res1.txt') 
    spec2score_dict = psite_result_read(psite_res_path, blind_summary_file_path, mod) 

    spec2score_list = []
    for spec, score in spec2score_dict.items(): 
        spec2score_list.append([spec, score]) 
    
    # 用于卡值，可以换成其他策略
    cut_off_ratio = 10
    spec2score_list.sort(key=lambda s:s[1]) 
    cut_off_num = int(len(spec2score_list) * (1 - cut_off_ratio / 100.0)) 
    spec_name_list = []
    for i in range(cut_off_num): 
        spec_name_list.append(spec2score_list[i][0]) 
    
    position_counter = position_static(mod, spec_name_list, blind_summary_file_path, parameter_dict['side_position'])
    
    return local_list_combine(local_list, position_counter)



if __name__ == "__main__":  
    '''
    blind_pfind_summary_path = 'results/source/blind/pFind-Filtered.spectra' 
    mod = 'PFIND_DELTA_387.18' 
    spec2pos = psite_file_generate(blind_pfind_summary_path, mod) 
    psite_res_path = 'results/source/pSite/res1.txt' 
    spec2score = psite_result_read(psite_res_path)
    '''
    current_path = os.getcwd() 
    cfg_path = os.path.join(current_path, 'pChem.cfg')
    mod = '387' 
    parameter_dict = parameter_file_read(cfg_path) 
    local_list = [('C', 106), ('N-SIDE', 5), ('A', 3), ('G', 2), ('F', 2), ('T', 1), ('D', 1), ('E', 1), ('K', 1), ('N', 1), ('W', 1)]
    print(psite_run(parameter_dict, current_path, mod, 'blind', local_list))


