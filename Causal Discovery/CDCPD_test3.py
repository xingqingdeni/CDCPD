from __future__ import print_function
import warnings
import sys
from mosum import bandwidths_default
from mosum.bandwidth import multiscale_grid
from mosum.bootstrap import confint_multiscale_cpts, confint_mosum_cpts
from mosum.classes import multiscale_cpts
from mosum.mosum_test import pValue
import pandas as pd
import ruptures as rpt
from ruptures import Pelt
from analyze import analyze_causal_relations
from tigramite.toymodels.structural_causal_processes import structural_causal_process
from tigramite.data_processing import DataFrame
import math
import numpy as np
import tigramite.data_processing as pp
from matplotlib import pyplot as plt
from tigramite.independence_tests.cmiknn import CMIknn
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr
from mosum import mosum
from causallearn.search.ConstraintBased.CDNOD import cdnod

random_state = np.random.default_rng(seed=2025)

def multiscale_struct_bottomUp(x, G=None, threshold=['critical_value', 'custom'][0],
                        alpha=0.1, threshold_function=None, eta=0.4, do_confint=False, level=0.05, N_reps=1000, distance = None, min_jump = 0):
    """
    Multiscale MOSUM algorithm with bottom-up merging

    Parameters
    ----------
    x : list
        input data
    G : int
        vector of bandwidths; given as either integers less than `len(x)/2`,
         or numbers between `0` and `0.5` describing the moving sum bandwidths relative to `len(x)`
    threshold : Str
        indicates which threshold should be used to determine significance.
        By default, it is chosen from the asymptotic distribution at the given significance level 'alpha`.
        Alternatively it is possible to parse a user-defined function with 'threshold_function'.
    alpha : float
        numeric value for the significance level with '0 <= alpha <= 1';
        use iff 'threshold = "critical_value"'
    threshold_function : function
    eta : float
        a positive numeric value for the minimal mutual distance of changes,
        relative to moving sum bandwidth (iff 'criterion = "eta"')
    do_confint : bool
        flag indicating whether to compute the confidence intervals for change points
    level : float
        use iff 'do_confint = True'; a numeric value ('0 <= level <= 1') with which '100(1-level)%'
        confidence interval is generated
    N_reps : int
        use iff 'do.confint = True'; number of bootstrap replicates to be generated

    Returns
    -------
    multiscale_cpts object containing
    x : list
        input data
    G : int
        bandwidth vector
    threshold, alpha, threshold_function, eta
        input
    cpts : ndarray
        estimated change point
    cpts_info : DataFrame
        information on change points, including detection bandwidths, asymptotic p-values, scaled jump sizes
    pooled_cpts : ndarray
        change point candidates
    do_confint : bool
        input
    ci
        confidence intervals
    """
    n = len(x)
    if G is None:
        G = bandwidths_default(n, G_min=max(20, np.ceil(0.05 * n)))
        grid = multiscale_grid(G, method='concatenate')
    elif type(G) in [int, float]:
        grid = multiscale_grid([G], method='concatenate')
    elif type(G) == 'multiscale_grid_obj':
        if any(G.grid[1] - G.grid[0] != 0): sys.exit("Expecting a grid of symmetric bandwidths")
        grid = G
    elif type(G) == list:
        G.sort()
        grid = multiscale_grid(G, method='concatenate')
    else:
        sys.exit('Expecting a vector of numbers')
    abs_bandwidth = (np.array(grid.grid) >= 1).all()

    if abs_bandwidth:
        GRID_THRESH = max([20, 0.05 * n])
    else:
        GRID_THRESH = 0.05

    if (threshold == 'critical_value') & (min(grid.grid[0]) < GRID_THRESH):
        warnings.warn(
            'Smallest bandwidth in grid is relatively small (in comparison to n), \n increase the smallest bandwidth or use multiscale.localPrune instead')

    if (not threshold == 'critical_value') & (not threshold == 'custom'):
        warnings.warn('threshold must be either \'critical.value\' or \'custom\'')

    if not (alpha >= 0) & (alpha <= 1): sys.exit("alpha out of range")
    if not (eta <= 1) & (eta > 0): sys.exit("eta out of range")
    if not (not do_confint or N_reps > 0): sys.exit()

    # Retreive change point candidates from all bandwidths.
    cpts_complete = []
    bandwidths_complete = []
    pValues_complete = []
    jumps_complete = []

    GG = len(grid.grid[0])
    for i in range(GG):
        G = grid.grid[0][i]
        if threshold == 'critical_value':
            m = mosum(x, G, threshold='critical_value', alpha=alpha, criterion='eta', eta=eta)
        else:
            threshold_val = threshold_function(G, n, alpha)
            m = mosum(x, G, threshold='custom', threshold_custom=threshold_val, alpha=alpha, criterion='eta', eta=eta)
        if not abs_bandwidth:
            G = int(np.floor(G * n))
        if GG >= 2:
            cpts = m.cpts
            cpts_complete = np.append(cpts_complete, cpts)
            bandwidths_complete = np.append(bandwidths_complete, np.full(len(cpts), G))
            pValues_complete = np.append(pValues_complete, pValue(m.stat[cpts], n, G))
            jumps_complete = np.append(jumps_complete, m.stat[cpts] * np.sqrt(2 / G))

    # Merge candidates.
    if GG >= 2:
        points = [0]
        bandwidths = []
        pValues = []
        jumps = []
        cptsInOrder = range(len(cpts_complete))
        for i in cptsInOrder:
            p = cpts_complete[i]
            G = bandwidths_complete[i]
            pVal = pValues_complete[i]
            jmp = jumps_complete[i]
            # print(np.abs(p-points))
            if min(np.concatenate((np.abs(p - points), float("inf")),
                                  axis=None)) >= distance and jmp > min_jump:  # Note: min(empty_list) = Inf
                points = np.append(points, p)
                bandwidths = np.append(bandwidths, G)
                pValues = np.append(pValues, pVal)
                jumps = np.append(jumps, jmp)
        cpts_merged = pd.DataFrame({"cpts": points[1:], "G_left": bandwidths, "G_right": bandwidths,
                                    "p_value": pValues, "jump": jumps})
        cpts = cpts_merged["cpts"][cpts_merged.cpts.argsort().argsort()]
        G = grid.grid[0]
        if not abs_bandwidth:
            G = np.floor(n * G)
        out = multiscale_cpts(x, cpts, cpts_merged, np.sort(np.unique(cpts_complete)), G,
                              alpha, threshold, threshold_function, 'eta', eta,
                              False, None)  # note
        if do_confint:
            out.ci = confint_multiscale_cpts(out, level=level, N_reps=N_reps)
            out.do_confint = True
    else:
        out = m
        if do_confint:
            out.ci = confint_mosum_cpts(out, level=level, N_reps=N_reps)
            out.do_confint = True
    return out
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

# 修改后的links_1，每个变量有自身延迟一阶影响且无环
links_1 = {
    0: [((0, -1), 0.6, lin_f)],
    1: [((1, -1), 0.5, lin_f), ((0, -1), 0.3, nonlin_f)],
    2: [((2, -1), 0.4, lin_f), ((1, -1), 0.2, lin_f), ((0, -2), -0.2, lin_f)],
    3: [((3, -1), 0.3, lin_f), ((2, -1), 0.2, lin_f), ((1, 0), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f), ((2, 0), 0.7, lin_f)],#改了这块
}
# 修改后的links_2，满足要求
links_2 = {
    0: [((0, -1), 0.6, lin_f), ((0, -2), 0.6, nonlin_f)],
    1: [((1, -1), 0.2, lin_f), ((0, -1), 0.3, nonlin_f)],
    2: [((2, -1), 0.7, lin_f), ((1, -1), 0.5, lin_f), ((0, -1), -0.5, lin_f)],
    3: [((3, -1), 0.7, lin_f), ((2, -1), 0.2, lin_f), ((1, 0), -0.2, lin_f)],
    4: [((4, -1), 0.2, lin_f)],
}

noises = [random_state_1.standard_normal, random_state_2.standard_normal, random_state_3.standard_normal,random_state_4.standard_normal,
          random_state_5.standard_normal]
ens = 5
data_ens = {}
all_data = []
sum = 7000
for j in range(7):
    for i in range(ens):
        if j == 0 or j == 2 or j == 4 or j == 6:
            links = links_1
            T = 1000
        else:
            links = links_2
            T = 1000

        data, nonstat = structural_causal_process(links,
                                              T, noises=noises,transient_fraction=0)
    all_data.append(data)
data = np.concatenate(all_data, axis=0)
dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2', r'X3', r'X4'])
# Define different colors for each variable
colors = ['#0000FF', '#FF8C00', '#FF0000', '#00FF00', '#800080']
# Determine a fixed y-axis range that encompasses all data
y_min = np.min(data)
y_max = np.max(data)
y_range = [y_min, y_max]
# Create a figure with subplots for each variable
fig, axs = plt.subplots(5, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
# Plot time series with different colors for each variable
for i, color in enumerate(colors):
    axs[i].plot(data[:, i], color=color, label=r'$X^{}'.format(i))
    axs[i].set_ylabel(dataframe.var_names[i], fontname='Times New Roman')
    axs[i].set_ylim(y_range)  # Set the same y-axis range for all subplots
# Hide x-axis ticks, labels, and spines for the first three subplots
for ax in axs[:-1]:
    ax.set_xticks([])
    ax.set_xticklabels([])
    ax.spines['bottom'].set_visible(False)
# Set xticks, labels, and spine for the last subplot
axs[-1].set_xticks(np.arange(0, sum, step=200))
axs[-1].set_xticklabels([str(int(tick)) for tick in np.arange(0, sum, step=200)], fontname='Times New Roman')
axs[-1].set_xlabel('Time', fontname='Times New Roman')
axs[-1].spines['bottom'].set_visible(True)
# Adjust layout to prevent overlap
plt.tight_layout()
# Set the spines to only show the left frame for all subplots
for ax in axs:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
plt.show()

#变点检测
first_column_data = data[:,2]
first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# 将秩次转换为从1开始的连续整数
first_column_data_ranks = first_column_data_ranks

#变点检测pelt
algo_pelt= Pelt(model="l2",min_size=400,jump = 255)#设置模型为L1范数，最小突变点长度为2
#检测变点
bkps_pelt = algo_pelt.fit(first_column_data)
breakpoints =algo_pelt.predict(pen=100)
print("Using PELT algorithm:")
print("Detected breakpoints using L2 norm:", breakpoints)

# 打印秩次
# print(first_column_data_ranks)
# detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
result_mosum = multiscale_struct_bottomUp(first_column_data_ranks, G = [220,300,400,500], eta = 2/3, alpha=0.001, distance = 800)
# result_mosum = mosum.multiscale_localPrune(first_column_data, G = [200], criterion="eta", eta = 2/3, alpha=0.01)
# summary and print methods
# result_mosum.summary()
cpts = result_mosum.cpts.tolist()
result_mosum.print()

# plot the output
true_cpts = [1000, 2000, 3000, 4000, 5000, 6000, 7000]
bkps = true_cpts
result = cpts
rpt.display(data, bkps, result)
plt.show()


segments = []
start = 0
cpts.append(sum)
cpts = sorted(cpts)
for cpt in cpts:
    segment = data[start:cpt, :]
    segments.append(segment)
    start = cpt
# PCMCIPLUS
ci_test = CMIknn(significance="fixed_thres", verbosity=3)
ci_test = ParCorr()
tau_max = 2
for i, segment in enumerate(segments):
    segment_df = pp.DataFrame(segment, var_names=[r'X0', r'X1', r'X2', r'X3', r'X4'])
    pcmci_parcorr = PCMCI(
        dataframe=segment_df,
        cond_ind_test=ci_test,
        verbosity=1)
    print(f"Results for segment {i + 1}:")
    results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
                                          pc_alpha=0.001,
                                          reset_lagged_links=False,
                                          )
#CDNOD
c_indx = np.arange(1, sum+1).reshape(-1, 1) # 将指定列索引对应的列数据提取出来，形成二维数组，列数为1
# 使用默认参数构建因果图
cg_default = cdnod(data, c_indx)
# 可视化默认参数构建的因果图（展示）
cg_default.draw_pydot_graph()

#模型评估
# 定义一个列表来存储所有阶段的因果关系字典（真实因果关系）
true_stages_causal = []
# 循环生成三个阶段的因果关系字典并添加到列表中（真实因果关系）
for i in range(1, 8):
    if i == 1 or i == 3 or i == 5 or i == 7:
        stage_causal = {
            "X0": ["X0[t - 1]"],
            "X1": ["X1[t - 1]", "X0[t - 1]"],
            "X2": ["X2[t - 1]", "X1[t - 1]", "X0[t - 2]"],
            "X3": ["X3[t - 1]", "X2[t - 1]", "X0[t]"],
            "X4": ["X4[t - 1]", "X3[t - 1]", "X2[t]"],
        }
    else:
        stage_causal = {
            "X0": ["X0[t - 1]", "X0[t - 2]"],
            "X1": ["X1[t - 1]", "X0[t - 1]"],
            "X2": ["X2[t - 1]", "X1[t - 1]", "X0[t - 1]"],
            "X3": ["X3[t - 1]", "X2[t - 1]", "X1[t]"],
            "X4": ["X4[t - 1]"],
        }
    true_stages_causal.append(stage_causal)

tested_causal = """
### 第一阶段
    Variable X0 has 1 link(s):
        (X0 -1): pval = 0.00000 | val =  0.534
    Variable X1 has 2 link(s):
        (X0 -1): pval = 0.00000 | val =  0.531
        (X1 -1): pval = 0.00000 | val =  0.435
    Variable X2 has 2 link(s):
        (X1 -1): pval = 0.00000 | val =  0.233
        (X2 -1): pval = 0.00000 | val =  0.224
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.298
        (X1  0): pval = 0.00000 | val = -0.234
        (X2 -1): pval = 0.00000 | val =  0.157
    Variable X4 has 2 link(s):
        (X2  0): pval = 0.00000 | val =  0.557
        (X4 -1): pval = 0.00000 | val =  0.151
### 第二阶段
    Variable X0 has 3 link(s):
        (X1 -1): pval = 0.00000 | val =  0.647
        (X0 -1): pval = 0.00000 | val =  0.510
        (X0 -2): pval = 0.00005 | val =  0.128
    Variable X1 has 1 link(s):
        (X1 -1): pval = 0.00000 | val =  0.175
    Variable X2 has 3 link(s):
        (X0 -1): pval = 0.00000 | val = -0.626
        (X1 -1): pval = 0.00000 | val =  0.606
        (X2 -1): pval = 0.00000 | val =  0.582
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.591
        (X1  0): pval = 0.00000 | val = -0.304
        (X2 -1): pval = 0.00000 | val =  0.198
    Variable X4 has 1 link(s):
        (X4 -1): pval = 0.00000 | val =  0.158
### 第三阶段
    Variable X0 has 1 link(s):
        (X0 -1): pval = 0.00000 | val =  0.477
    Variable X1 has 2 link(s):
        (X0 -1): pval = 0.00000 | val =  0.525
        (X1 -1): pval = 0.00000 | val =  0.522
    Variable X2 has 3 link(s):
        (X2 -1): pval = 0.00000 | val =  0.325
        (X1 -1): pval = 0.00000 | val =  0.225
        (X0 -1): pval = 0.00007 | val = -0.123
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.266
        (X1  0): pval = 0.00000 | val = -0.259
        (X2 -1): pval = 0.00000 | val =  0.221
    Variable X4 has 2 link(s):
        (X2  0): pval = 0.00000 | val =  0.518
        (X4 -1): pval = 0.00000 | val =  0.260
### 第四阶段
    Variable X0 has 4 link(s):
        (X1 -1): pval = 0.00000 | val =  0.661
        (X0 -1): pval = 0.00000 | val =  0.511
        (X1 -2): pval = 0.00000 | val = -0.189
        (X1  0): pval = 0.00000 | val = -0.152
    Variable X1 has 0 link(s):
    Variable X2 has 3 link(s):
        (X0 -1): pval = 0.00000 | val = -0.633
        (X2 -1): pval = 0.00000 | val =  0.597
        (X1 -1): pval = 0.00000 | val =  0.547
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.546
        (X1  0): pval = 0.00000 | val = -0.232
        (X2 -1): pval = 0.00000 | val =  0.200
    Variable X4 has 1 link(s):
        (X4 -1): pval = 0.00000 | val =  0.233
### 第五阶段
   Variable X0 has 1 link(s):
        (X0 -1): pval = 0.00000 | val =  0.539
    Variable X1 has 4 link(s):
        (X1 -1): pval = 0.00000 | val =  0.454
        (X0 -1): pval = 0.00000 | val =  0.386
        (X2  0): pval = 0.00002 | val =  0.125
        (X0  0): pval = 0.00003 | val = -0.123
    Variable X2 has 4 link(s):
        (X4  0): pval = 0.00000 | val =  0.421
        (X2 -1): pval = 0.00000 | val =  0.397
        (X0 -1): pval = 0.00000 | val = -0.272
        (X1 -1): pval = 0.00000 | val =  0.224
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.390
        (X2 -1): pval = 0.00000 | val =  0.225
        (X1  0): pval = 0.00000 | val = -0.216
    Variable X4 has 1 link(s):
        (X4 -1): pval = 0.00000 | val =  0.306
### 第六阶段
   Variable X0 has 3 link(s):
        (X1 -1): pval = 0.00000 | val =  0.669
        (X0 -1): pval = 0.00000 | val =  0.517
        (X1 -2): pval = 0.00001 | val = -0.155
    Variable X1 has 0 link(s):
    Variable X2 has 3 link(s):
        (X0 -1): pval = 0.00000 | val = -0.664
        (X1 -1): pval = 0.00000 | val =  0.605
        (X2 -1): pval = 0.00000 | val =  0.597
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.575
        (X1  0): pval = 0.00000 | val = -0.313
        (X2 -1): pval = 0.00000 | val =  0.202
    Variable X4 has 0 link(s):
### 第七阶段
     Variable X0 has 1 link(s):
        (X0 -1): pval = 0.00000 | val =  0.484
    Variable X1 has 4 link(s):
        (X1 -1): pval = 0.00000 | val =  0.468
        (X0 -1): pval = 0.00000 | val =  0.457
        (X2  0): pval = 0.00001 | val =  0.131 | unoriented link
        (X0  0): pval = 0.00003 | val = -0.127
    Variable X2 has 5 link(s):
        (X2 -1): pval = 0.00000 | val =  0.471
        (X4  0): pval = 0.00000 | val =  0.454
        (X0 -1): pval = 0.00000 | val = -0.182
        (X4 -1): pval = 0.00000 | val = -0.139
        (X1  0): pval = 0.00001 | val =  0.131 | unoriented link
    Variable X3 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.345
        (X1  0): pval = 0.00000 | val = -0.257
        (X2 -1): pval = 0.00000 | val =  0.186
    Variable X4 has 1 link(s):
        (X4 -1): pval = 0.00000 | val =  0.341
"""
true_counts, tested_counts, same_counts_matrix = analyze_causal_relations(true_stages_causal, tested_causal)
#CDCPD_precision
precision = (9/12)*(6300/7000)+(400/7000)*(7/12)
print(f"precision的值为: {precision}")
recall = (9/10)*(6300/7000)+(400/7000)*(7/10)
print(f"recall的值为: {recall}")
f1_score = 2*(precision * recall)/(precision+recall)
print(f"precision的值为: {f1_score}")











# for col in range(data.shape[1]):
#     plt.plot(data[:, col], label=f'Column {col}')
# plt.xlabel('Index')
# plt.ylabel('Value')
# plt.title('Data of Each Column')
# plt.legend()
# plt.show()
#
# dataFrame = pp.DataFrame(data)
# cond_ind_test = ParCorr()
# # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
#
# pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# results = pcmci.run_pcmci(tau_max=1, pc_alpha = None)
# pcmci.print_significant_links(p_matrix=results['p_matrix'],
#                                          val_matrix=results['val_matrix'],
#                                          alpha_level=0.05)


#var_process数据生成2
# links = {0: [((0, -1), 0.6)],
#          1: [((1, -1), 0.8), ((0, -1), 0.3), ((2, 0), 0.3)],
#          2: [((2, -1), 0.7), ((1, 0), 0.3)]}
# data, _ = toys.var_process(links, T=1000)
# dataFrame = pp.DataFrame(data)
# cond_ind_test = ParCorr()
# pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# results = pcmci.run_pcmci(tau_max=2, pc_alpha=None)
# pcmci.print_significant_links(p_matrix=results['p_matrix'],
#                                          val_matrix=results['val_matrix'],
#                                          alpha_level=0.05)

#数据生成3
# links_coeffs = {0: [((0, -1), 0.7, lin_f)],
#                 1: [((1, -1), 0.7, lin_f), ((0, 0), 0.2, lin_f), ((2, -2), 0.2, lin_f)],
#                 2: [((2, -1), 0.3, lin_f)],
#                 }
# T = 1000     # time series length
# data, _ = toys.structural_causal_process(links_coeffs, T=T, seed=3)
# T, N = data.shape
# dataFrame = pp.DataFrame(data)
# cond_ind_test = ParCorr()
# # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# results = pcmci.run_pcmci(tau_max=2, pc_alpha = 0.2)
# pcmci.print_significant_links(p_matrix=results['p_matrix'],
#                                          val_matrix=results['val_matrix'],
#                                          alpha_level=0.05)

# #数据生成4
# # Example process to play around with
# T = 1500
# data = np.random.default_rng().normal(1, 1, size = (T, 4))
# # data = random_state.standard_normal((T, 4))
# c = 0.8
# for t in range(0, 500):
#     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
#     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
#     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
#     data[t, 3] += 0.8 * data[t-1, 3]
# for t in range(500, 1000):
#     data[t, 0] += 0.8 * data[t-1, 0] + 0.2 * data[t-1, 2]
#     data[t, 1] += 0.5 * data[t-1, 1]
#     data[t, 2] += 0.3 * data[t-1, 2] + 0.4 * data[t-2, 1]
#     data[t, 3] += 0.3 * data[t-1, 3]
# for t in range(1000, T):
#     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
#     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
#     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
#     data[t, 3] += 0.8 * data[t-1, 3]
# dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2', 'Sun'])
# print(data)
# # Define different colors for each variable
# colors = ['blue', 'orange', 'red', 'green']
# # Determine a fixed y-axis range that encompasses all data
# y_min = np.min(data)
# y_max = np.max(data)
# y_range = [y_min, y_max]
# # Create a figure with subplots for each variable
# fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
# # Plot time series with different colors for each variable
# for i, color in enumerate(colors):
#     axs[i].plot(data[:, i], color=color, label=r'$X^{}'.format(i))
#     axs[i].set_ylabel(dataframe.var_names[i], fontname='Times New Roman')
#     axs[i].set_ylim(y_range)  # Set the same y-axis range for all subplots
# # Hide x-axis ticks, labels, and spines for the first three subplots
# for ax in axs[:-1]:
#     ax.set_xticks([])
#     ax.set_xticklabels([])
#     ax.spines['bottom'].set_visible(False)
# # Set xticks, labels, and spine for the last subplot
# axs[-1].set_xticks(np.arange(0, T, step=200))
# axs[-1].set_xticklabels([str(int(tick)) for tick in np.arange(0, T, step=200)], fontname='Times New Roman')
# axs[-1].set_xlabel('Time', fontname='Times New Roman')
# axs[-1].spines['bottom'].set_visible(True)
# # Adjust layout to prevent overlap
# plt.tight_layout()
# # Set the spines to only show the left frame for all subplots
# for ax in axs:
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
# plt.show()
# #变点检测
# first_column_data = data[:, 0]
# first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# # 将秩次转换为从1开始的连续整数
# first_column_data_ranks = first_column_data_ranks
# # 打印秩次
# # print(first_column_data_ranks)
# # detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# # result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
# #result_mosum = mosum.multiscale_bottomUp(first_column_data_ranks, G = [200,250,300,350], eta = 2/3, alpha=0.1)
# result_mosum = mosum.multiscale_localPrune(first_column_data_ranks, G = [20,40,60,80,100], criterion="epsilon", epsilon = 2/3, alpha=0.01)
# # summary and print methods
# # result_mosum.summary()
# cpts = result_mosum.cpts.tolist()
# result_mosum.print()
# # plot the output
# result_mosum.plot()
# plt.show()
#
# segments = []
# start = 0
# cpts.append(T)
# for cpt in cpts:
#     segment = data[start:cpt, :]
#     segments.append(segment)
#     start = cpt
#
# # PCMCIPLUS
# ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# ci_test = ParCorr()
# tau_max = 2
# link_assumptions = {}
# for j in range(4):
#     if j in [0, 1, 2]:
#         # Directed lagged links
#         link_assumptions[j] = {(var, -lag): '-?>' for var in [0, 1, 2]
#                                for lag in range(1, tau_max + 1)}
#         # Unoriented contemporaneous links
#         link_assumptions[j].update({(var, 0): 'o?o' for var in [0, 1, 2] if var!= j})
#         # Directed lagged and contemporaneous links from the sun (3)
#         link_assumptions[j].update({(var, -lag): '-?>' for var in [3]
#                                     for lag in range(0, tau_max + 1)})
#     else:
#         link_assumptions[j] = {}
#
# all_results = []
# for i, segment in enumerate(segments):
#     segment_df = pp.DataFrame(segment, var_names=[r'X0', r'X1', r'X2', 'Sun'])
#     pcmci_parcorr = PCMCI(
#         dataframe=segment_df,
#         cond_ind_test=ci_test,
#         verbosity=1)
#     print(f"Results for segment {i + 1}:")
#     results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
#                                           pc_alpha=0.001,
#                                           reset_lagged_links=False,
#                                           link_assumptions=link_assumptions)



# # PCMCIPLUS
# ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# ci_test = ParCorr()
# tau_max = 2
# link_assumptions = {} #{}
# for j in range(4):
#     if j in [0, 1, 2]:
#         # Directed lagged links
#         link_assumptions[j] = {(var, -lag): '-?>' for var in [0, 1, 2]
#                          for lag in range(1, tau_max + 1)}
#         # Unoriented contemporaneous links
#         link_assumptions[j].update({(var, 0): 'o?o' for var in [0, 1, 2] if var != j})
#         # Directed lagged and contemporaneous links from the sun (3)
#         link_assumptions[j].update({(var, -lag): '-?>' for var in [3]
#                          for lag in range(0, tau_max + 1)})
#     else:
#         link_assumptions[j] = {}
#
# for j in link_assumptions:
#     print(link_assumptions[j])
# pcmci_parcorr = PCMCI(
#     dataframe=dataframe,
#     cond_ind_test=ci_test,
#     verbosity=1)
# results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
#                 pc_alpha=0.001,
#                 reset_lagged_links=False,
#                 link_assumptions=link_assumptions
#                 ) #, alpha_level = 0.01)
# print(results['graph'].shape)
# # print(results['graph'][:,3,:])
# print(np.round(results['p_matrix'][:,:,0], 2))
# print(np.round(results['val_matrix'][:,:,0], 2))
# print(results['graph'][:,:,0])
#
#
# algo_pelt= Pelt(model="l2",jump = 255)#设置模型为L1范数，最小突变点长度为2
# #检测变点
# bkps_pelt = algo_pelt.fit(first_column_data)
# breakpoints =algo_pelt.predict(pen=0.01)
# print("Using PELT algorithm:")
# print("Detected breakpoints using L2 norm:", breakpoints)


#
# T = 1000
# data = random_state.standard_normal((T, 4))
# # Simple sun
# data[:, 3] = random_state.standard_normal((T))  # np.sin(np.arange(T)*20/np.pi) + 0.1*random_state.standard_normal((T))
# c = 0.8
# for t in range(1, T):
#     data[t, 0] += 0.4 * data[t - 1, 0] + 0.4 * data[t - 1, 1] + c * data[t - 1, 3]
#     data[t, 1] += 0.5 * data[t - 1, 1] + c * data[t, 3]
#     data[t, 2] += 0.6 * data[t - 1, 2] + 0.3 * data[t - 2, 1]  # + c*data[t-1,3]
# dataframe = pp.DataFrame(data, var_names=[r'X^0', r'X^1', r'X^2', 'Sun'])
# tp.plot_timeseries(dataframe); plt.show()

