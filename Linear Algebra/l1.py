import re
import numpy as np

#把增广矩阵转换成行阶梯形（REF）
def ref(A,b):
    n = len(A)
    # 前向消元
    for i in range(n):
        # 寻找主元
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k
        # 交换最大行至当前行
        A[i], A[maxRow] = A[maxRow], A[i]
        b[i], b[maxRow] = b[maxRow], b[i]
        # 将当前行的首元素归一化，并消去下方所有行的对应元素
        for k in range(i+1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
            b[k] += c * b[i]
    return A, b

#高斯消元法求解线性方程组
def gaussian_elimination(A, b):
    n = len(A)
    # 前向消元
    for i in range(n):
        # 寻找主元
        maxEl = abs(A[i][i])
        maxRow = i
        for k in range(i+1, n):
            if abs(A[k][i]) > maxEl:
                maxEl = abs(A[k][i])
                maxRow = k
        # 交换最大行至当前行
        A[i], A[maxRow] = A[maxRow], A[i]
        b[i], b[maxRow] = b[maxRow], b[i]
        # 将当前行的首元素归一化，并消去下方所有行的对应元素
        for k in range(i+1, n):
            c = -A[k][i] / A[i][i]
            for j in range(i, n):
                if i == j:
                    A[k][j] = 0
                else:
                    A[k][j] += c * A[i][j]
            b[k] += c * b[i]
    # 回代
    x = [0 for _ in range(n)]
    for i in range(n-1, -1, -1):
        x[i] = b[i] / A[i][i]
        for k in range(i-1, -1, -1):
            b[k] -= A[k][i] * x[i]
    return x

#输入方程组
equations = []
T=int(input('How many equations do you have?\n'))
print('Enter the equation like 3x+4y+z=80:')
for i in range(T):
    equations.append(input())

# 步骤1: 收集所有变量名
all_variables = set()
pattern = r"([a-zA-Z]+)"
for eq in equations:
    variables = re.findall(pattern, eq.split('=')[0])
    all_variables.update(variables)
all_variables = sorted(list(all_variables)) # 排序以保持一致性

# 步骤2: 构建方程系数矩阵
coeff_matrix = []
for eq in equations:
    coeff = [0] * len(all_variables)  #初始化系数为0
    for var in all_variables:
        # 查找变量系数，如果没有找到，则系数为0
        match = re.search(rf"(-?\d*){var}", eq)
        if match:
            coeff_val = match.group(1)
            coeff[all_variables.index(var)] = int(coeff_val) if coeff_val not in ['', '-'] else -1 if coeff_val == '-' else 1
    coeff_matrix.append(coeff)

# 步骤3: 构建常数项向量
const_vector = [int(eq.split('=')[1]) for eq in equations]

# 步骤4: 根据行阶梯形REF的秩，判断方程组的解的情况
A = coeff_matrix
b = const_vector
print('系数矩阵：',A)
print('常数项向量：',b)
A, b = ref(A, b)#把增广矩阵转换成行阶梯形（REF）
augmented_matrix=np.column_stack((A,b))#合成增广矩阵
rank1 = np.linalg.matrix_rank(A)#求系数矩阵的秩
rank2=np.linalg.matrix_rank(augmented_matrix)#求增广矩阵的秩
n=len(all_variables)#未知数的个数
if rank1==rank2:
    if rank1==n:
        print('The system has a unique solution')#有唯一解
#步骤5: 调用高斯消元法求解方程组，输出结果
        x = gaussian_elimination(A, b)
        print('解为：',x)
    else:
        print('The system has infinite solutions')#有无穷多解
else:
    print('The system has no solution')#无解