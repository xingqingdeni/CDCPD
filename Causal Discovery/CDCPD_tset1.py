from __future__ import print_function
import numpy as np
import pandas as pd
from detecta import detect_cusum
import tigramite.data_processing as pp
from matplotlib import pyplot as plt
import ruptures as rpt
from analyze import analyze_causal_relations, calculate_precision, wbs
from tigramite.independence_tests.cmiknn import CMIknn
from tigramite.pcmci import PCMCI
from tigramite.independence_tests.parcorr import ParCorr
import mosum
from ruptures import Pelt, Binseg
from causallearn.search.ConstraintBased.CDNOD import cdnod
from causallearn.utils.GraphUtils import GraphUtils
# 数据生成1
## Generate some time series from a structural causal process
def lin_f(x): return x
def nonlin_f(x): return (x + 5. * x ** 2 * np.exp(-x ** 2 / 20.))
T = 500
sum = 3*T
data = np.random.default_rng(seed=1111).normal(1, 1, size = (sum, 4))
c = 0.8
for t in range(0, T):
    data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
    data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
    data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
    data[t, 3] += 0.8 * data[t-1, 3]
for t in range(T, 2*T):
    data[t, 0] += 0.8 * data[t-1, 0] + 0.2 * data[t-1, 2]
    data[t, 1] += 0.5 * data[t-1, 1]
    data[t, 2] += 0.3 * data[t-1, 2] + 0.4 * data[t-2, 1]
    data[t, 3] += 0.3 * data[t-1, 3]
for t in range(2*T, sum):
    data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
    data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
    data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
    data[t, 3] += 0.8 * data[t-1, 3]
dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2', 'X3'])
print(data)
# Define different colors for each variable
colors = ['blue', 'orange', 'red', 'green']
# Determine a fixed y-axis range that encompasses all data
y_min = np.min(data)
y_max = np.max(data)
y_range = [y_min, y_max]
# Create a figure with subplots for each variable
fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
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
first_column_data = data[:, 0]
first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# 将秩次转换为从1开始的连续整数
first_column_data_ranks = first_column_data_ranks
# pelt
algo_pelt= Pelt(model="l2",min_size=400,jump = 255)#设置模型为L1范数，最小突变点长度为2
bkps_pelt = algo_pelt.fit(first_column_data)
breakpoints =algo_pelt.predict(pen=10)
print("Using PELT algorithm:")
print("Detected breakpoints:", breakpoints)
# NMCD
algo_NMCD= rpt.Dynp(model="rank")
algo = algo_NMCD.fit(first_column_data)
breakpoints = algo.predict(n_bkps = 2)
print("Using NMCD algorithm:")
print("Detected breakpoints:", breakpoints)
# CUSUM
breakpoints = detect_cusum(first_column_data, threshold=10, drift=0.15, show=False)
print("Using CUSUM algorithm:")
print("Detected breakpoints:", breakpoints)
# WBS
breakpoints = wbs(first_column_data, M=10, n_bkps=3)
print("Using WBS algorithm:")
print("Detected breakpoints:", breakpoints)
# PCMCI
ci_test = CMIknn(significance="fixed_thres", verbosity=3)
ci_test = ParCorr()
tau_max = 1
pcmci_parcorr = PCMCI(
    dataframe=dataframe,
    cond_ind_test=ci_test,
    verbosity=1)
print(f"Results for PCMCI:")
results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
                                          pc_alpha=0.01,
                                          reset_lagged_links=False,
                                          )











# 打印秩次
# print(first_column_data_ranks)
# detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
#result_mosum = mosum.multiscale_bottomUp(first_column_data_ranks, G = [200,250,300,350], eta = 2/3, alpha=0.1)
result_mosum = mosum.multiscale_localPrune(first_column_data_ranks, G = [20,40,60,80,100], criterion="epsilon", epsilon = 2/3, alpha=0.01)
# summary and print methods
# result_mosum.summary()
cpts = result_mosum.cpts.tolist()
result_mosum.print()

# plot the output
true_cpts = [500, 1000, 1500]
bkps = true_cpts
result = cpts
rpt.display(data, bkps, result)
plt.show()

segments = []
start = 0
cpts.append(sum)
for cpt in cpts:
    segment = data[start:cpt, :]
    segments.append(segment)
    start = cpt
# PCMCIPLUS
ci_test = CMIknn(significance="fixed_thres", verbosity=3)
ci_test = ParCorr()
tau_max = 1
for i, segment in enumerate(segments):
    segment_df = pp.DataFrame(segment, var_names=[r'X0', r'X1', r'X2', 'X3'])
    pcmci_parcorr = PCMCI(
        dataframe=segment_df,
        cond_ind_test=ci_test,
        verbosity=1)
    print(f"Results for segment {i + 1}:")
    results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
                                          pc_alpha=0.01,
                                          reset_lagged_links=False,
                                          )
#模型评估
# 定义一个列表来存储所有阶段的因果关系字典（真实因果关系）
true_stages_causal = []
# 循环生成三个阶段的因果关系字典并添加到列表中（真实因果关系）
for i in range(1, 4):
    if i == 1:
        stage_causal = {
            "X0": ["X0[t - 1]", "X1[t - 1]", "X3[t - 1]"],
            "X1": ["X1[t - 1]", "X3[t]"],
            "X2": ["X2[t - 1]", "X1[t]"],
            "X3": ["X3[t - 1]"]
        }
    elif i == 2:
        stage_causal = {
            "X0": ["X0[t - 1]", "X2[t - 1]"],
            "X1": ["X1[t - 1]"],
            "X2": ["X2[t - 1]", "X1[t - 2]"],
            "X3": ["X3[t - 1]"]
        }
    else:
        stage_causal = {
            "X0": ["X0[t - 1]", "X1[t - 1]", "X3[t - 1]"],
            "X1": ["X1[t - 1]", "X3[t]"],
            "X2": ["X2[t - 1]", "X1[t]"],
            "X3": ["X3[t - 1]"]
        }
    true_stages_causal.append(stage_causal)
tested_causal = """
### 第一阶段
    Variable X0 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.519
        (X0 -1): pval = 0.00000 | val =  0.403
        (X1 -1): pval = 0.00000 | val =  0.339
    Variable X1 has 2 link(s):
        (X3  0): pval = 0.00000 | val =  0.574
        (X1 -1): pval = 0.00000 | val =  0.263
    Variable X2 has 2 link(s):
        (X2 -1): pval = 0.00000 | val =  0.564
        (X1  0): pval = 0.00000 | val =  0.345
    Variable X3 has 2 link(s):
        (X1  0): pval = 0.00000 | val =  0.574
        (X3 -1): pval = 0.00000 | val =  0.546
### 第二阶段
Variable X0 has 2 link(s):
        (X0 -1): pval = 0.00000 | val =  0.655
        (X2 -1): pval = 0.00000 | val =  0.253
    Variable X1 has 1 link(s):
        (X1 -1): pval = 0.00000 | val =  0.500
    Variable X2 has 1 link(s):
        (X2 -1): pval = 0.00000 | val =  0.358
    Variable X3 has 1 link(s):
        (X3 -1): pval = 0.00000 | val =  0.325
### 第三阶段
Variable X0 has 3 link(s):
        (X3 -1): pval = 0.00000 | val =  0.507
        (X1 -1): pval = 0.00000 | val =  0.413
        (X0 -1): pval = 0.00000 | val =  0.397
    Variable X1 has 2 link(s):
        (X3  0): pval = 0.00000 | val =  0.624
        (X1 -1): pval = 0.00000 | val =  0.330
    Variable X2 has 2 link(s):
        (X2 -1): pval = 0.00000 | val =  0.433
        (X1  0): pval = 0.00000 | val =  0.320
    Variable X3 has 2 link(s):
        (X1  0): pval = 0.00000 | val =  0.624
        (X3 -1): pval = 0.00000 | val =  0.591
"""
true_counts, tested_counts, same_counts_matrix = analyze_causal_relations(true_stages_causal, tested_causal)
precision = calculate_precision(cpts, true_cpts,sum,same_counts_matrix,tested_counts)
print(f"precision的值为: {precision}")
recall = (500/1500)*(8/8)+(5/1500)*(4/6)+(495/1500)*(5/6)+(3/1500)*(4/8)+(497/1500)*(8/8)
print(f"recall的值为: {recall}")
f1_score = 2*(precision * recall)/(precision+recall)
print(f"f1_score: {f1_score}")