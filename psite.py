
import os 
# D:\pSite\pPredictAA\x64\Release

def psite_file_generate(file_path, target_mod): 
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

    with open('psite.txt', 'w', encoding='utf-8') as f: 
        for line in new_lines: 
            f.write(line) 
    return spect2pos 


# 读取psite的结果文件
def psite_result_read(file_path): 
    with open(file_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    spec2score = {}
    #i = 0
    #total = 0
    for line in lines: 
        line = line.split('\t') 
        
        score = float(line[4])
        spec = line[0] 
        #if score < 5.0:
        #    continue
        #print(score, spec2pos[spec]) 
        #if spec2pos[spec][0] == 'C':
        #    i += 1  
        #total += 1 
    #print(float(i/total)) 
        spec2score[spec] = score 
    
    return spec2score


def psite_run():
    return 1


if __name__ == "__main__":  
    blind_pfind_summary_path = 'results/source/blind/pFind-Filtered.spectra' 
    mod = 'PFIND_DELTA_387.18' 
    spec2pos = psite_file_generate(blind_pfind_summary_path, mod) 
    psite_res_path = 'results/source/pSite/res1.txt' 
    spec2score = psite_result_read(psite_res_path) 
    current_path = os.getcwd() 
    cfg_path = os.path.join(current_path, 'pChem.cfg')
    psite_run(cfg_path, current_path) 


