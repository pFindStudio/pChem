from argparse import ArgumentParser 


element_dict={
    "C": 12.0000000,
    "H": 1.0078246,
    "N": 14.0030732,
    "O": 15.9949141,
    "S": 31.972070,
}

result = []

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
def fast_dfs(i, number_list, max_number_list, p, left, right): 

    if mass_in_range(number_list, p, left, right):
        #result.append(number_list)
        return 
    if i >= len(number_list):
        return 
    for k in range(max_number_list[i]+1):
        number_list[i] = k 
        fast_dfs(i+1, number_list, max_number_list, p, left, right)
    return


def mass_in_range(number_list, p, left, right):
    current_mass = 0.0
    for i in range(len(number_list)):
        current_mass += element_dict[p[i]]*number_list[i]
    # print(number_list, current_mass)
    if current_mass >= left and current_mass <= right:
        result.append(number_list.copy())
        return True
    return False 


def main():
    parser = ArgumentParser()
    parser.add_argument("--mass", type=float, default=100.0030732, help="the mass for composition analysis")
    parser.add_argument("--element_list", type=str, default="CHNOS", help="the element used for inference")
    parser.add_argument("--error_range", type=float, default=20.0, help="the error range (ppm)") 
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
    fast_dfs(0, number_list, max_number_list, p, left, right) 
    print(p)
    for line in result:
        print(line)




if __name__ == "__main__":
    main() 

