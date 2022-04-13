A = ["r", "s", "t", "u", "v", "w", "x", "y", "z"]
B = [ 0,   1,   1,    0,   1,   2,   2,   0,   1]
result_list = [i for _,i in sorted(zip(B,A))]
print(result_list)

A = ['nMT', "E", "PCP", "AT"]
B = [ 'PCP', 'E', 'AT', 'nMET', 'PCP', 'TE']
result_list = [i for _,i in sorted(zip(B,A))]
print(result_list)