import os 

raw1_list=['QE_Plus_YJ_FL_50per_20190501_F1_R1.raw', 'QE_Plus_YJ_HJX_NPM_PH8_50per_20190820_F1_R2.raw',
            'QE_Plus_YJ_HJX_ENE8_50per_20191002_F1_R1.raw', 'HFX_YangJing_HJX_VSF_R1_20210319.raw',
            'QE_Plus_YJ_HJX_3_PPMS-P8_50per_20190822_F1_R1.raw',] 
raw_list= ['HFX_YangJing_HeJiXiang_STP-1_20210515_R1.raw', 
            'HFX_YangJing_HJX_NHS100_20210607_R1.raw','HFX_YangJing_HeJiXiang_Dia_1_20210401.raw',  
            'QE_Plus_YangJing_APEX_L2_Elution_TCP_100per_20180815.raw', 'QE_Plus_TCP_HRP_Tyr_20200113_F1_R1.raw',
            'SulfenQ_H2O2_LH_C1.raw', 'QE_Plus_YangJing_FL_SOH-55_50per_20170618.raw',
            'QE_Plus_YangJing_WYne_C_TCP_50per_20180823.raw', 'QE_Plus_YangJing_WYne_N_TCP_50per_20180823.raw',
            'QE_Plus_YangJing_WYne_O_TCP_50per_20180823.raw',  'QE_Plus_YangJing_FL_ALK_50per_20170531.raw',
            'RKO_aHNE_PC0h_BR1_T1.raw','JY_RKO_aONE_PC0h_20150419_BR1_T1.raw']


prob_list = ['IPM', 'NPM', 'ENE', 'VSF', 'PPMS', 'STP', 'NHS', 'Diazirine', 'APEX', 'HRP', 
            'Dyn-2', 'BTD', 'WyneC', 'WyneN', 'WyneO', 'DiaAlk', 'Ac4ManNAz', 'aHNE', 'aONE'] 


def parameter_cfg_read(cfg_path): 
    with open(cfg_path, 'r', encoding='utf-8') as f: 
        lines = f.readlines() 
    return lines 


def run_pchem(): 
    cfg_path = 'pChem_null.cfg'  
    cfg_lines = parameter_cfg_read(cfg_path)
    
    for raw in raw_list: 
        new_lines = [] 
        for line in cfg_lines: 
            t_line = line
            if 'output_path=' in line: 
                t_line = t_line[:-1] + raw[:-4] + '\n' 
            elif 'msmspath1=' in line:
                t_line = t_line[:-1] + raw + '\n' 
            new_lines.append(t_line) 
        print(new_lines) 
        with open('pChem.cfg', 'w', encoding='utf-8') as f: 
            for l in new_lines: 
                f.write(l) 
        cmd = 'pChem.exe'
        receive = os.system(cmd) 
        print(receive) 



if __name__ == "__main__":  
    run_pchem()