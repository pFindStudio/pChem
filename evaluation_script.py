import os 
import matplotlib.pyplot as plt 
import numpy as np 
from pchem import run 
from curve_fit import mass_read, ppm_calculate 


wangchu_list = ['20181026_1TO1_1.raw', '20181026_1TO2_1.raw', '20181026_1TO5_1.raw', '20181026_1TO10_1.raw'] 
wangchu_mass_list = [[521.30743, 527.32124],
                    [521.30743, 527.32124],
                    [521.30743, 527.32124],
                    [521.30743, 527.32124],]


raw_list=['QE_Plus_YJ_FL_50per_20190501_F1_R1.raw', 'QE_Plus_YJ_HJX_NPM_PH8_50per_20190820_F1_R2.raw',
            'QEy_Plus_YJ_HJX_ENE8_50per_20191002_F1_R1.raw', 'HFX_YangJing_HJX_VSF_R1_20210319.raw',
            'QE_Plus_YJ_HJX_3_PPMS-P8_50per_20190822_F1_R1.raw', 'HFX_YangJing_HeJiXiang_STP-1_20210515_R1.raw', 
            'HFX_YangJing_HJX_NHS100_20210607_R1.raw', 'HFX_YangJing_HeJiXiang_Dia_1_20210401.raw',  
            'QE_Plus_YangJing_APEX_L2_Elution_TCP_100per_20180815.raw', 'QE_Plus_TCP_HRP_Tyr_20200113_F1_R1.raw',
            'SulfenQ_H2O2_LH_C1.raw', 
            'QE_Plus_YangJing_FL_SOH-55_50per_20170618.raw',
            'QE_Plus_YangJing_WYne_C_TCP_50per_20180823.raw', 'QE_Plus_YangJing_WYne_N_TCP_50per_20180823.raw',
            'QE_Plus_YangJing_WYne_O_TCP_50per_20180823.raw',  'QE_Plus_YangJing_FL_ALK_50per_20170531.raw',
            'RKO_aHNE_PC0h_BR1_T1.raw','JY_RKO_aONE_PC0h_20150419_BR1_T1.raw']


gt_mass_list = [[252.12224, 258.142372],
                [292.11716, 298.137292, 310.12772, 316.147852],
                [279.15829, 285.178422],
                [315.12528, 321.145412],
                [227.07285, 233.092982],
                [237.11134, 243.131472], 
                [251.12699, 257.147122],
                [267.15829, 273.178422],
                [372.17976, 378.199892, 209.08004, 215.100172],
                [372.17976, 378.199892, 388.17467, 394.194802],
                [333.16886, 339.206525],
                [418.13109, 424.151222],
                [265.14264, 271.162772, 511.20193, 517.222062],
                [252.12224, 258.142372],
                [267.12191, 273.142042, 291.12191, 297.142042],
                [387.1754, 393.195532, 471.23291, 477.253042],
                [311.18451, 317.204642, 309.16886, 315.188992],
                [346.16411, 352.184242, 289.14264, 295.162772, 307.15321, 313.173342, 578.22846, 584.248592]]

prob_list = ['IPM', 'NPM', 'ENE', 'VSF', 'PPMS', 'STP', 'NHS', 'Diazirine', 'APEX', 'HRP', 
            'Dyn-2', 'BTD', 'WyneC', 'WyneN', 'WyneO', 'DiaAlk', 'Ac4ManNAz', 'aHNE', 'aONE'] 


assert len(raw_list) == len(gt_mass_list)


def config_write(cfg_path, raw_name): 
    with open(cfg_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    output_path = os.path.join('D:\\pChem\\pChem_new\\results', raw_name[:-4]) 
    msms_path = os.path.join('D:\\pChem\\pchem_data', raw_name) 
    new_lines = []
    for line in lines:
        if 'output_path' in line: 
            new_lines.append('output_path='+output_path+'\n') 
        elif 'msmspath1' in line: 
            new_lines.append('msmspath1='+msms_path+'\n')
        else:
            new_lines.append(line) 
    with open(cfg_path, 'w', encoding='utf-8') as f: 
        for line in new_lines:
            f.write(line) 
    return output_path  



def dataset_list_evaluation():
    current_path = os.getcwd() 
    cfg_path = os.path.join(current_path, 'pChem.cfg') 

    computed_mass_list = [] 
    ppm_list = [] 
    slam_ppm_list = []
    slam_number = []

    for i in range(len(raw_list)):
        raw_name = raw_list[i] 
        output_path = config_write(cfg_path, raw_name)
        run()
        blind_res_path = os.path.join(os.path.join(output_path, 'blind'), 'pFind-Filtered.spectra')
        t_computed_list = [] 
        t_ppm_list = [] 
        j = 0
        for target_mass in gt_mass_list[i]: 
            t_comp_mass, t_number = mass_read(blind_res_path, round(target_mass, 2)) 
            t_ppm = ppm_calculate(t_comp_mass, target_mass) 
            print(t_comp_mass)
            print(t_ppm)
            slam_ppm_list.append(t_ppm)
            slam_number.append(t_number)
            t_computed_list.append(t_comp_mass) 
            t_ppm_list.append(t_ppm) 
            j += 1
            
        computed_mass_list.append(t_computed_list)
        ppm_list.append(t_ppm_list)
        # break 
    print(ppm_list)
    print(computed_mass_list) 
    

    write_lines = [] 
    for i in range(len(ppm_list)): 
        line = prob_list[i] + '\t'
        for j in range(len(ppm_list[i])):
            line += str(computed_mass_list[i][j]) + '\t'
            line += str(ppm_list[i][j]) + '\t' 
        line += '\n' 
        write_lines.append(line) 

    with open('evaluation_res', 'w', encoding='utf-8') as f: 
        for line in write_lines:
            f.write(line) 
    
    # scatter_plot(slam_number, slam_ppm_list)
    #new_ppm_list = [] 
    #for i in range(len(slam_ppm_list)): 
    #    if slam_number[i] >= 500: 
    #        new_ppm_list.append(slam_ppm_list[i])
    print(slam_ppm_list)
    print(np.mean(slam_ppm_list))


def scatter_plot(x, y, x_label=None, y_label=None): 
    colors='#DC143C' 
    plt.scatter(x,y,c=colors) 
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.show()



def close_dataset_list_evaluation():
    current_path = os.getcwd() 
    cfg_path = os.path.join(current_path, 'pChem.cfg') 

    computed_mass_list = [] 
    ppm_list = [] 
    slam_ppm_list = []
    slam_number = []

    for i in range(len(raw_list)):
        raw_name = raw_list[i] 
        output_path = config_write(cfg_path, raw_name)
        run()
        blind_res_path = os.path.join(os.path.join(output_path, 'source\\blind'), 'pFind-Filtered.spectra')
        t_computed_list = [] 
        t_ppm_list = [] 
        j = 0
        for target_mass in gt_mass_list[i]: 
            t_comp_mass, t_number = mass_read(blind_res_path, round(target_mass, 3)) 
            t_ppm = ppm_calculate(t_comp_mass, target_mass) 
            print(t_comp_mass)
            print(t_ppm)
            slam_ppm_list.append(t_ppm)
            slam_number.append(t_number)
            t_computed_list.append(t_comp_mass) 
            t_ppm_list.append(t_ppm) 
            j += 1
            
        computed_mass_list.append(t_computed_list)
        ppm_list.append(t_ppm_list)
        # break 
    print(ppm_list)
    print(computed_mass_list) 
    

    write_lines = [] 
    for i in range(len(ppm_list)): 
        line = prob_list[i] + '\t'
        for j in range(len(ppm_list[i])):
            line += str(computed_mass_list[i][j]) + '\t'
            line += str(ppm_list[i][j]) + '\t' 
        line += '\n' 
        write_lines.append(line) 

    with open('evaluation_res', 'w', encoding='utf-8') as f: 
        for line in write_lines:
            f.write(line) 
    
    # scatter_plot(slam_number, slam_ppm_list)
    #new_ppm_list = [] 
    #for i in range(len(slam_ppm_list)): 
    #    if slam_number[i] >= 500: 
    #        new_ppm_list.append(slam_ppm_list[i])
    print(slam_ppm_list)
    print(np.mean(slam_ppm_list))


if __name__ == "__main__":  
    close_dataset_list_evaluation()



