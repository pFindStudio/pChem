import os

element_dict={
    "C": 12.0000000,
    "H": 1.00782503207,
    "Pm": 1.00727647012,
    "N": 14.0030740048,
    "O": 15.99491461956,
    "S": 31.972071
}

amino_acid_dict={
    "A" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "B" : element_dict["C"]*0,
    "C" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "D" : element_dict["C"]*4 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "E" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*3 + element_dict["S"]*0,
    "F" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "G" : element_dict["C"]*2 + element_dict["H"]*3 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "H" : element_dict["C"]*6 + element_dict["H"]*7 + element_dict["N"]*3 + element_dict["O"]*1 + element_dict["S"]*0,
    "I" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "J" : element_dict["C"]*0,
    "K" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "L" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "M" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*1,
    "N" : element_dict["C"]*4 + element_dict["H"]*6 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "O" : element_dict["C"]*0,
    "P" : element_dict["C"]*5 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Q" : element_dict["C"]*5 + element_dict["H"]*8 + element_dict["N"]*2 + element_dict["O"]*2 + element_dict["S"]*0,
    "R" : element_dict["C"]*6 + element_dict["H"]*12 + element_dict["N"]*4 + element_dict["O"]*1 + element_dict["S"]*0,
    "S" : element_dict["C"]*3 + element_dict["H"]*5 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "T" : element_dict["C"]*4 + element_dict["H"]*7 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "U" : element_dict["C"]*0,
    "V" : element_dict["C"]*5 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "W" : element_dict["C"]*11 + element_dict["H"]*10 + element_dict["N"]*2 + element_dict["O"]*1 + element_dict["S"]*0,
    "X" : element_dict["C"]*6 + element_dict["H"]*11 + element_dict["N"]*1 + element_dict["O"]*1 + element_dict["S"]*0,
    "Y" : element_dict["C"]*9 + element_dict["H"]*9 + element_dict["N"]*1 + element_dict["O"]*2 + element_dict["S"]*0,
    "Z" : element_dict["C"]*0, 
}

h2o_mass = element_dict["H"]*2 + element_dict["O"]*1
proton_mass = element_dict['Pm']


# 读取选择的常见修饰，返回dict
def common_dict_create(current_path):
    modification_path = os.path.join(current_path, 'modification-null.ini')
    with open(modification_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    i = 1
    common_dict = {} 
    while i < len(lines):
        if len(lines[i]) < 2:
            break
        mod_name = lines[i].split()[0]
        eq_idx = mod_name.find('=')
        mod_name = mod_name[eq_idx+1:]
        mod_mass = lines[i+1].split()[2]
        common_dict[mod_name] = float(mod_mass)
        i += 2
    return common_dict
