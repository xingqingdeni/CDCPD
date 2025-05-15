from __future__ import print_function
from ruptures import Binseg
from tigramite.toymodels.structural_causal_processes import structural_causal_process
import math
import numpy as np
random_state = np.random.default_rng(seed=2025)
import matplotlib.pyplot as plt
import numpy as np
import ruptures as rpt


# 数据生成1
## Generate some time series from a structural causal process
def lin_f(x): return x
def nonlin_f(x): return (x ** 1/2 + 5 * math.sin(x))


random_state_1 = np.random.default_rng(300)
random_state_2 = np.random.default_rng(400)
random_state_3 = np.random.default_rng(500)
random_state_4 = np.random.default_rng(600)
random_state_5 = np.random.default_rng(700)
random_state_6 = np.random.default_rng(800)
random_state_7 = np.random.default_rng(900)
# 修改后的links_1，每个变量有自身延迟一阶影响且无环，延迟最高阶为4，结构与其他不同
links_1 = {
    0: [((0, -1), 0.6, lin_f), ((0, -4), 0.1, nonlin_f)],  # 增加了变量0的 -4阶依赖，满足延迟最高阶为4
    1: [((1, -1), 0.5, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -3), 0.1, lin_f)],  # 增加对变量0的 -3阶依赖，不超过最高阶4
    2: [((2, -1), 0.4, lin_f), ((1, -1), 0.2, lin_f), ((0, -2), -0.2, lin_f)],  # 包含不同延迟阶数依赖，最高阶4
    3: [((3, -1), 0.3, lin_f), ((2, -1), 0.2, lin_f), ((1, -1), -0.2, lin_f)],  # 加入变量1的 -4阶依赖，符合要求
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.7, lin_f), ((2, -2), 0.1, nonlin_f)],  # 保持对其他变量不同延迟阶依赖，最高阶4
}

# 修改后的links_2，满足要求，与其他结构不同，延迟最高阶为4
links_2 = {
    0: [((0, -1), 0.6, lin_f), ((0, -2), 0.4, nonlin_f), ((0, -4), 0.1, lin_f)],  # 体现不同延迟阶数，最高阶4
    1: [((1, -1), 0.2, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -4), 0.1, nonlin_f)],  # 增加对变量0的 -4阶依赖
    2: [((2, -1), 0.7, lin_f), ((1, -1), 0.5, lin_f), ((0, -3), -0.5, lin_f)],  # 包含不同变量不同延迟阶依赖，最高阶4
    3: [((3, -1), 0.7, lin_f), ((2, -1), 0.2, lin_f), ((1, 0), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.1, lin_f), ((2, -4), 0.3, lin_f)],  # 加入变量2的 -4阶依赖，满足最高阶4要求
}

# 修改后的links_3，确保每个变量含自身延迟一阶影响且无环，结构区别于其他，延迟最高阶为4
links_3 = {
    0: [((0, -1), 0.6, lin_f), ((0, -3), 0.2, nonlin_f), ((0, -4), 0.1, lin_f)],  # 多种延迟阶数组合，最高阶4
    1: [((1, -1), 0.5, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -2), 0.2, nonlin_f)],  # 合理利用不同延迟阶，不超最高阶
    2: [((2, -1), 0.4, lin_f), ((1, -1), 0.2, lin_f), ((0, -3), -0.2, lin_f)],
    3: [((3, -1), 0.3, lin_f), ((2, -1), 0.2, lin_f), ((1, -1), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.7, lin_f), ((2, -4), 0.1, nonlin_f)],  # 加入变量2的 -4阶依赖，符合最高阶要求
}

# 修改后的links_4，符合条件，与其余结构不同，延迟最高阶为4
links_4 = {
    0: [((0, -1), 0.6, lin_f), ((0, -4), 0.1, nonlin_f)],
    1: [((1, -1), 0.2, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -4), 0.1, lin_f)],
    2: [((2, -1), 0.7, lin_f), ((1, -1), 0.5, lin_f), ((0, -4), -0.5, lin_f)],  # 不同延迟阶数组合，最高阶4
    3: [((3, -1), 0.7, lin_f), ((2, -1), 0.2, lin_f), ((1, 0), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.1, lin_f), ((2, -4), 0.3, lin_f)],
}

# 修改后的links_5，维持要求特性，结构与其他不同，延迟最高阶为4
links_5 = {
    0: [((0, -1), 0.6, lin_f), ((1, -2), 0.2, nonlin_f), ((0, -4), 0.1, lin_f)],
    1: [((1, -1), 0.5, lin_f), ((0, -1), 0.3, nonlin_f)],
    2: [((2, -1), 0.4, lin_f), ((1, -1), 0.2, lin_f), ((0, -4), -0.2, lin_f)],
    3: [((3, -1), 0.3, lin_f), ((2, -1), 0.2, lin_f), ((1, -1), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.7, lin_f), ((2, -4), 0.1, nonlin_f)],
}

# 修改后的links_6，满足每个变量有自身延迟一阶影响及无环，与其他结构有差异，延迟最高阶为4
links_6 = {
    0: [((0, -1), 0.6, lin_f),  ((0, -4), 0.2, lin_f)],
    1: [((1, -1), 0.2, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -3), 0.1, lin_f)],
    2: [((2, -1), 0.7, lin_f), ((1, -1), 0.5, lin_f), ((0, -4), -0.5, lin_f)],
    3: [((3, -1), 0.7, lin_f), ((2, -1), 0.2, lin_f), ((1, 0), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.1, lin_f), ((2, -4), 0.3, lin_f)],
}

# 修改后的links_7，保证符合要求，结构不同于其他，延迟最高阶为4
links_7 = {
    0: [((0, -1), 0.6, lin_f), ((0, -4), 0.1, lin_f)],
    1: [((1, -1), 0.5, lin_f), ((0, -1), 0.3, nonlin_f), ((0, -3), 0.2, nonlin_f)],
    2: [((2, -1), 0.4, lin_f), ((1, -1), 0.2, lin_f), ((0, -4), -0.2, lin_f)],
    3: [((3, -1), 0.3, lin_f), ((2, -1), 0.2, lin_f), ((1, -1), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((3, -1), 0.7, lin_f), ((2, -4), 0.1, nonlin_f)],
}
noises = [random_state_1.standard_normal, random_state_2.standard_normal, random_state_3.standard_normal,random_state_4.standard_normal,
          random_state_5.standard_normal]
ens = 5
data_ens = {}
all_data = []
sum = 7000
for j in range(len([links_1, links_2, links_3, links_4, links_5, links_6, links_7])):
    for i in range(ens):
        if j == 0:
            links = links_1
            T = 1000
        elif j == 1:
            links = links_2
            T = 1000
        elif j == 2:
            links = links_3
            T = 1000
        elif j == 3:
            links = links_4
            T = 1000
        elif j == 4:
            links = links_5
            T = 1000
        elif j == 5:
            links = links_6
            T = 1000
        else:
            links = links_7
            T = 1000
        data, nonstat = structural_causal_process(links,
                                              T, noises=noises,transient_fraction=0)
    all_data.append(data)
data = np.concatenate(all_data, axis=0)

# PELT
algo_pelt= rpt.Pelt(model="rbf",min_size=400,jump = 255)#设置模型为L1范数，最小突变点长度为2
algo = algo_pelt.fit(data)
result = algo.predict(pen = 10)
true_cpts = [1000, 2000, 3000, 4000, 5000, 6000, 7000]
bkps = true_cpts
rpt.display(data, bkps, result)
plt.show()

# Dynp
# algo_Dynp= rpt.Dynp(model="l2",min_size=400,jump = 255)
# algo = algo_Dynp.fit(data)
# result = algo.predict(n_bkps = 2)
# true_cpts = [1500, 3000, 4500]
# bkps = true_cpts
# rpt.display(data, bkps, result)
# plt.show()

# BottomUp
algo_bottomUp= rpt.BottomUp(model="rbf",min_size=400,jump = 255)
algo = algo_bottomUp.fit(data)
result = algo.predict(n_bkps = 6,pen = 10)
true_cpts = [1000, 2000, 3000, 4000, 5000, 6000, 7000]
bkps = true_cpts
rpt.display(data, bkps, result)
plt.show()

# Binseg
binseg_obj = Binseg(model="rbf")
algo = binseg_obj.fit(data)
result = algo.predict(n_bkps = 6)
true_cpts = [1000, 2000, 3000, 4000, 5000, 6000, 7000]
bkps = true_cpts
rpt.display(data, bkps, result)
plt.show()

# window;
algo = rpt.Window(model="rbf").fit(data)
result = algo.predict(n_bkps = 6)
true_cpts = [1000, 2000, 3000, 4000, 5000, 6000, 7000]
bkps = true_cpts
rpt.display(data, bkps, result)
plt.show()
