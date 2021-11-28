from argparse import ArgumentParser 
import numpy as np
import matplotlib.pyplot as plt

element_dict={
    "C": 12.0000000,
    "H": 1.0078246,
    "N": 14.0030732,
    "O": 15.9949141,
    "S": 31.972070,
}

result = [] 
absolute_error = []

def range_calculate(mass, error_range):
    error_range = error_range / 1000000.0
    left = mass / ( 1 + error_range) 
    right = mass / (1 - error_range) 
    return left, right


# 会有大量冗余的计算 
def dfs(current_p, res_mass, p, left, right): 
    print(current_p)
    if res_mass > right:
        return 
    if res_mass >= left and res_mass <= right: 
        result.append(current_p) 
        return
    for a in p:
        dfs(current_p+a, res_mass + element_dict[a], p, left, right)
    return 


# 优化后的计算方式 
def fast_dfs(i, number_list, max_number_list, p, left, right, target_mass): 

    if mass_in_range(number_list, p, left, right, target_mass):
        #result.append(number_list)
        return 
    if i >= len(number_list):
        return 
    for k in range(1, max_number_list[i]+1):
        number_list[i] = k 
        fast_dfs(i+1, number_list, max_number_list, p, left, right, target_mass)
    return


# 根据绝对误差计算分子式 
def ab_fast_dfs(i, number_list, max_number_list, p, targetmass, gamma): 

    if i == len(number_list):
        current_mass = 0.0 
        for j in range(len(number_list)):
            current_mass += element_dict[p[j]]*number_list[j] 
        if abs(current_mass - targetmass) < gamma: 
            result.append(number_list.copy()) 
        return 
    if i > len(number_list):
        return 
    for k in range(1, max_number_list[i]+1):
        number_list[i] = k 
        ab_fast_dfs(i+1, number_list, max_number_list, p, targetmass, gamma)
    return



def mass_in_range(number_list, p, left, right, target_mass):
    current_mass = 0.0
    for i in range(len(number_list)):
        current_mass += element_dict[p[i]]*number_list[i]
    # print(number_list, current_mass)
    if current_mass >= left and current_mass <= right:
        result.append(number_list.copy()) 
        # absolute_error.append(target_mass - current_mass)
        return True
    return False 



def main():
    parser = ArgumentParser()
    parser.add_argument("--mass", type=float, default=252.122, help="the mass for composition analysis")
    parser.add_argument("--element_list", type=str, default="CHNO", help="the element used for inference")
    parser.add_argument("--error_range", type=float, default=5.0, help="the error range (ppm)") 
    args = parser.parse_args()
    
    # 检查输入是否正确 
    p = ""
    for a in args.element_list:
        if a in element_dict.keys():
            p += a 
    
    # 计算误差范围
    left, right = range_calculate(args.mass, args.error_range)
    #print(left, right)
    # current_p = "" 
    # dfs(current_p, 0.0, p, left, right)
    # print(result) 

    number_list = [0] * len(p)
    max_number_list = []
    for a in p:
        max_number_list.append(int(args.mass/element_dict[a])+1)
    # print(max_number_list)
    print('begin to calculation the combination....')
    # print(left, right)
    #fast_dfs(0, number_list, max_number_list, p, left, right, args.mass) 
    ab_fast_dfs(0, number_list, max_number_list, p, args.mass)
    #print(p)
    #print(len(result))
    for i in range(len(result)):
        print(result[i])
        #print(round(absolute_error[i],5))


def static_for_18_probs(): 
    
    mass_list = [252.1222, 292.1172, 279.1583, 269.1376, 
                    336.1182, 315.1253, 227.0728, 511.2025,
                    513.1817, 512.1977, 346.1641, 311.1845,
                    333.1689, 418.1311, 387.1754, 238.1192,
                    372.1789, 267.1583, 279.1583, 292.1172]
    ppm_range_list = [x for x in range(1, 21)] 
    p = 'ONCH' 
    global result 
    final_list = []
    for ppm_error in ppm_range_list: 
        t_number = []
        for mass in mass_list: 
            left, right = range_calculate(mass, ppm_error) 
            # print(left, right)
            number_list = [0] * len(p)
            max_number_list = []
            for a in p:
                max_number_list.append(int(mass/element_dict[a])+1)
            # print(max_number_list)
            # print(left, right)
            fast_dfs(0, number_list, max_number_list, p, left, right, mass) 
            t_number.append(len(result))
            print(t_number[-1])
            result = []
        final_list.append(t_number)
        print('-', len(final_list))
    print(len(final_list))
    
    
    fig, ax = plt.subplots()  

    # ax.boxplot(final_list) 
    ax.violinplot(final_list)
    ax.set_xlabel('ppm')
    ax.set_ylabel('candidate number')
    #ax.set_xticklabels(["girl20", "boy20", "girl30", "boy30"])     # 设置x轴刻度标签
    plt.show()
    

def ab_static_for_18_probs(): 
    
    mass_list = [252.1222, 292.1172, 279.1583, 269.1376, 
                    336.1182, 315.1253, 227.0728, 511.2025,
                    513.1817, 512.1977, 346.1641, 311.1845,
                    333.1689, 418.1311, 387.1754, 238.1192,
                    372.1789, 267.1583, 279.1583, 292.1172]
    ppm_range_list = [x/1000 for x in range(1, 21)] 

    p = 'ONCH' 
    labels = []
    for ppm in ppm_range_list: 
        labels.append(str(ppm))
    print(labels)
    global result 
    final_list = []
    for ppm_error in ppm_range_list: 
        t_number = []
        for mass in mass_list: 
            # left, right = range_calculate(mass, ppm_error) 
            # print(left, right)
            number_list = [0] * len(p)
            max_number_list = []
            for a in p:
                max_number_list.append(int(mass/element_dict[a])+1)
            # print(max_number_list)
            # print(left, right) 
            print(ppm_error, mass)
            ab_fast_dfs(0, number_list, max_number_list, p, mass, ppm_error) 
            t_number.append(len(result))
            print(t_number[-1])
            result = []
        final_list.append(t_number)
        print('-', len(final_list))
    print(len(final_list))
    
    
    fig, ax = plt.subplots()  


    ax.violinplot(final_list) 
    #ax.boxplot(final_list, labels=labels, showfliers=False) 
    ax.set_xlabel('absolute error')
    ax.set_ylabel('candidate number')
    #ax.set_xticklabels(["girl20", "boy20", "girl30", "boy30"])     # 设置x轴刻度标签
    plt.show()


if __name__ == "__main__":
    static_for_18_probs()
    #a=[[1,2,3],[1,4,5]]
    #fig, ax = plt.subplots()  
    #ax.violinplot(a)
    #plt.show()

