from __future__ import print_function
from mosum import bandwidths_default
from mosum.bandwidth import multiscale_grid
from mosum.bootstrap import confint_multiscale_cpts, confint_mosum_cpts
from mosum.classes import multiscale_cpts
import ruptures as rpt
from mosum.mosum_test import pValue

from tigramite.independence_tests.parcorr import ParCorr
from causallearn.search.ConstraintBased.CDNOD import cdnod
from mosum import mosum
from data_processing import DataFrame
import tigramite.data_processing as pp
from matplotlib import pyplot as plt
from independence_tests.cmiknn import CMIknn
from tigramite.pcmci import PCMCI
import tigramite.plotting as tp
import warnings
import sys
import numpy as np
import pandas as pd
from ruptures import Pelt
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
file_path = r'C:\Users\10362\Desktop\论文\数据集\fmri_data.xlsx'
# 使用NumPy的loadtxt函数读取数据
data = pd.read_excel(file_path)
data = data.to_numpy()
# # last_column = data[:, -1]
data = data[:, :-1]
dataframe = pp.DataFrame(data, var_names=[r"CALC","LIPL","LT","LTRIA","LOPER","LIPS","LDLPFC"])
# tp.plot_timeseries(dataframe);
# plt.show()
first_column_data = data[:,3]
# first_column_data = (first_column_data - np.min(first_column_data)) / (np.max(first_column_data) - np.min(first_column_data))
first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# 将秩次转换为从1开始的连续整数
first_column_data_ranks = first_column_data_ranks

#检测变点pelt
algo_pelt= Pelt(model="l2",min_size=523,jump = 255)#设置模型为L1范数，最小突变点长度为2
bkps_pelt = algo_pelt.fit(first_column_data)
breakpoints =algo_pelt.predict(pen=10)
print("Using PELT algorithm:")
print("Detected breakpoints using L2 norm:", breakpoints)
#rosum变点检测
# result_mosum = mosum.mosum(first_column_data_ranks, var_est_method= 'mosum', G=1000, criterion="eta", eta = 0.9, alpha=0.01, boundary_extension=True)
result_mosum = multiscale_struct_bottomUp(first_column_data_ranks, G = [220,400], eta = 2/3, alpha=0.01,distance = 600)
# import mosum
# result_mosum = mosum.multiscale_localPrune(first_column_data_ranks, G = [50,100], criterion="eta", eta = 2/3, alpha=0.1)
cpts = result_mosum.cpts.tolist()
result_mosum.print()
plt.figure(figsize=(10, 3))
# plot the output
true_cpts = [550, 1445, 1936]
bkps = true_cpts
result = cpts
rpt.display(data, bkps, result)
plt.xlim(0, 1936)
plt.show()

segments = []
start = 0
cpts.append(1936)
for cpt in cpts:
    segment = data[start:cpt, :]
    segments.append(segment)
    start = cpt
# PCMCIPLUS
ci_test = CMIknn(significance="fixed_thres", verbosity=3)
ci_test = ParCorr()
tau_max = 1
for i, segment in enumerate(segments):
    segment_df = pp.DataFrame(segment, var_names=[r"CALC","LIPL","LT","LTRIA","LOPER","LIPS","LDLPFC"])
    pcmci_parcorr = PCMCI(
        dataframe=segment_df,
        cond_ind_test=ci_test,
        verbosity=1)
    print(f"Results for segment {i + 1}:")
    results = pcmci_parcorr.run_pcmciplus(tau_min=0,tau_max=1,
                                          pc_alpha=0.05,
                                          reset_lagged_links=False,
                                          )
    tp.plot_graph(
        val_matrix=results['val_matrix'],
        graph=results['graph'],
        var_names=[r"CALC","LIPL","LT","LTRIA","LOPER","LIPS","LDLPFC"],
        link_colorbar_label='MCI',
        arrow_linewidth=3.0,
        curved_radius=0.3,
        );
    plt.show()

#CDNOD
sum = num_rows = data.shape[0]
c_indx = np.arange(1, sum+1).reshape(-1, 1) # 将指定列索引对应的列数据提取出来，形成二维数组，列数为1
# 使用默认参数构建因果图
cg_default = cdnod(data, c_indx)
# 可视化默认参数构建的因果图（展示）
cg_default.draw_pydot_graph()







# # 数据生成1
# ## Generate some time series from a structural causal process
# def lin_f(x): return x
# def nonlin_f(x): return (x ** 1/2 + 5 * math.sin(x))
#
# random_state_1 = np.random.default_rng(300)
# random_state_2 = np.random.default_rng(400)
# random_state_3 = np.random.default_rng(500)
# links_1 = {0: [((0, -1), 0.6, lin_f)],
#          1: [((1, -1), 0.8, lin_f), ((0, -1), 0.3, nonlin_f)],
#          2: [((2, -1), 0.7, lin_f), ((1, 0), -0.2, lin_f)],
#          }
# links_2 = {0: [((0, -1), 0.6, lin_f),((0, -2), 0.6, nonlin_f)],
#          1: [((1, -1), 0.2, lin_f), ((0, -2), 0.3, nonlin_f)],
#          2: [((2, -1), 0.7, lin_f), ((1, 0), -0.5, lin_f)],
#          }
# links_3 = {0: [((0, -1), 0.6, lin_f)],
#          1: [((1, -1), 0.8, lin_f), ((0, -1), 0.8, nonlin_f)],
#          2: [((2, -1), 0.5, lin_f), ((1, 0), -0.7, nonlin_f)],
#          }
# noises = [random_state_1.standard_normal, random_state_2.standard_normal, random_state_3.standard_normal]
# ens = 3
# data_ens = {}
# all_data = []
# T = 1500
# sum = 3*T
# for j in range(len([links_1, links_2, links_3])):
#     for i in range(ens):
#         if j == 0:
#             links = links_1
#         elif j == 1:
#             links = links_2
#         else:
#             links = links_3
#         data, nonstat = structural_causal_process(links,
#                                               T, noises=noises,transient_fraction=0)
#     all_data.append(data)
# data = np.concatenate(all_data, axis=0)
# dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2'])
# # Define different colors for each variable
# colors = ['blue', 'orange', 'red']
# # Determine a fixed y-axis range that encompasses all data
# y_min = np.min(data)
# y_max = np.max(data)
# y_range = [y_min, y_max]
# # Create a figure with subplots for each variable
# fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
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
# axs[-1].set_xticks(np.arange(0, sum, step=200))
# axs[-1].set_xticklabels([str(int(tick)) for tick in np.arange(0, sum, step=200)], fontname='Times New Roman')
# axs[-1].set_xlabel('Time', fontname='Times New Roman')
# axs[-1].spines['bottom'].set_visible(True)
# # Adjust layout to prevent overlap
# plt.tight_layout()
# # Set the spines to only show the left frame for all subplots
# for ax in axs:
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
# plt.show()
#
#
# #变点检测
# first_column_data = data[:,1]
# first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# # 将秩次转换为从1开始的连续整数
# first_column_data_ranks = first_column_data_ranks
# # 打印秩次
# # print(first_column_data_ranks)
# # detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# # result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
# result_mosum = multiscale_struct_bottomUp(first_column_data_ranks, G = [180,220], eta = 2/3, alpha=0.01)
# # result_mosum = mosum.multiscale_localPrune(first_column_data, G = [200], criterion="eta", eta = 2/3, alpha=0.01)
# # summary and print methods
# # result_mosum.summary()
# cpts = result_mosum.cpts.tolist()
# result_mosum.print()
# # plot the output
# result_mosum.plot()
# plt.show()
# segments = []
# start = 0
# cpts.append(sum)
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
# for j in range(3):
#     if j in [0, 1, 2]:
#         # Directed lagged links
#         link_assumptions[j] = {(var, -lag): '-?>' for var in [0, 1, 2]
#                                for lag in range(1, tau_max + 1)}
#         # Unoriented contemporaneous links
#         link_assumptions[j].update({(var, 0): 'o?o' for var in [0, 1, 2] if var!= j})
#     else:
#         link_assumptions[j] = {}
#
# all_results = []
# for i, segment in enumerate(segments):
#     segment_df = pp.DataFrame(segment, var_names=[r'X0', r'X1', r'X2'])
#     pcmci_parcorr = PCMCI(
#         dataframe=segment_df,
#         cond_ind_test=ci_test,
#         verbosity=1)
#     print(f"Results for segment {i + 1}:")
#     results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
#                                           pc_alpha=0.001,
#                                           reset_lagged_links=False,
#                                           link_assumptions=link_assumptions)
#
#
# #CDNOD
# c_indx = np.arange(1, sum+1).reshape(-1, 1) # 将指定列索引对应的列数据提取出来，形成二维数组，列数为1
# # 使用默认参数构建因果图
# cg_default = cdnod(data, c_indx)
# # 可视化默认参数构建的因果图（展示）
# cg_default.draw_pydot_graph()
#
#
# # for col in range(data.shape[1]):
# #     plt.plot(data[:, col], label=f'Column {col}')
# # plt.xlabel('Index')
# # plt.ylabel('Value')
# # plt.title('Data of Each Column')
# # plt.legend()
# # plt.show()
# #
# # dataFrame = pp.DataFrame(data)
# # cond_ind_test = ParCorr()
# # # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# #
# # pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# # results = pcmci.run_pcmci(tau_max=1, pc_alpha = None)
# # pcmci.print_significant_links(p_matrix=results['p_matrix'],
# #                                          val_matrix=results['val_matrix'],
# #                                          alpha_level=0.05)
#
#
# #var_process数据生成2
# # links = {0: [((0, -1), 0.6)],
# #          1: [((1, -1), 0.8), ((0, -1), 0.3), ((2, 0), 0.3)],
# #          2: [((2, -1), 0.7), ((1, 0), 0.3)]}
# # data, _ = toys.var_process(links, T=1000)
# # dataFrame = pp.DataFrame(data)
# # cond_ind_test = ParCorr()
# # pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# # results = pcmci.run_pcmci(tau_max=2, pc_alpha=None)
# # pcmci.print_significant_links(p_matrix=results['p_matrix'],
# #                                          val_matrix=results['val_matrix'],
# #                                          alpha_level=0.05)
#
# #数据生成3
# # links_coeffs = {0: [((0, -1), 0.7, lin_f)],
# #                 1: [((1, -1), 0.7, lin_f), ((0, 0), 0.2, lin_f), ((2, -2), 0.2, lin_f)],
# #                 2: [((2, -1), 0.3, lin_f)],
# #                 }
# # T = 1000     # time series length
# # data, _ = toys.structural_causal_process(links_coeffs, T=T, seed=3)
# # T, N = data.shape
# # dataFrame = pp.DataFrame(data)
# # cond_ind_test = ParCorr()
# # # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# # pcmci = PCMCI(dataframe = dataFrame, cond_ind_test=cond_ind_test)
# # results = pcmci.run_pcmci(tau_max=2, pc_alpha = 0.2)
# # pcmci.print_significant_links(p_matrix=results['p_matrix'],
# #                                          val_matrix=results['val_matrix'],
# #                                          alpha_level=0.05)
#
# # #数据生成4
# # # Example process to play around with
# # T = 1500
# # data = np.random.default_rng().normal(1, 1, size = (T, 4))
# # # data = random_state.standard_normal((T, 4))
# # c = 0.8
# # for t in range(0, 500):
# #     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
# #     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
# #     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
# #     data[t, 3] += 0.8 * data[t-1, 3]
# # for t in range(500, 1000):
# #     data[t, 0] += 0.8 * data[t-1, 0] + 0.2 * data[t-1, 2]
# #     data[t, 1] += 0.5 * data[t-1, 1]
# #     data[t, 2] += 0.3 * data[t-1, 2] + 0.4 * data[t-2, 1]
# #     data[t, 3] += 0.3 * data[t-1, 3]
# # for t in range(1000, T):
# #     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
# #     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
# #     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
# #     data[t, 3] += 0.8 * data[t-1, 3]
# # dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2', 'Sun'])
# # print(data)
# # # Define different colors for each variable
# # colors = ['blue', 'orange', 'red', 'green']
# # # Determine a fixed y-axis range that encompasses all data
# # y_min = np.min(data)
# # y_max = np.max(data)
# # y_range = [y_min, y_max]
# # # Create a figure with subplots for each variable
# # fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
# # # Plot time series with different colors for each variable
# # for i, color in enumerate(colors):
# #     axs[i].plot(data[:, i], color=color, label=r'$X^{}'.format(i))
# #     axs[i].set_ylabel(dataframe.var_names[i], fontname='Times New Roman')
# #     axs[i].set_ylim(y_range)  # Set the same y-axis range for all subplots
# # # Hide x-axis ticks, labels, and spines for the first three subplots
# # for ax in axs[:-1]:
# #     ax.set_xticks([])
# #     ax.set_xticklabels([])
# #     ax.spines['bottom'].set_visible(False)
# # # Set xticks, labels, and spine for the last subplot
# # axs[-1].set_xticks(np.arange(0, T, step=200))
# # axs[-1].set_xticklabels([str(int(tick)) for tick in np.arange(0, T, step=200)], fontname='Times New Roman')
# # axs[-1].set_xlabel('Time', fontname='Times New Roman')
# # axs[-1].spines['bottom'].set_visible(True)
# # # Adjust layout to prevent overlap
# # plt.tight_layout()
# # # Set the spines to only show the left frame for all subplots
# # for ax in axs:
# #     ax.spines['right'].set_visible(False)
# #     ax.spines['top'].set_visible(False)
# # plt.show()
# # #变点检测
# # first_column_data = data[:, 0]
# # first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# # # 将秩次转换为从1开始的连续整数
# # first_column_data_ranks = first_column_data_ranks
# # # 打印秩次
# # # print(first_column_data_ranks)
# # # detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# # # result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
# # #result_mosum = mosum.multiscale_bottomUp(first_column_data_ranks, G = [200,250,300,350], eta = 2/3, alpha=0.1)
# # result_mosum = mosum.multiscale_localPrune(first_column_data_ranks, G = [20,40,60,80,100], criterion="epsilon", epsilon = 2/3, alpha=0.01)
# # # summary and print methods
# # # result_mosum.summary()
# # cpts = result_mosum.cpts.tolist()
# # result_mosum.print()
# # # plot the output
# # result_mosum.plot()
# # plt.show()
# #
# # segments = []
# # start = 0
# # cpts.append(T)
# # for cpt in cpts:
# #     segment = data[start:cpt, :]
# #     segments.append(segment)
# #     start = cpt
# #
# # # PCMCIPLUS
# # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# # ci_test = ParCorr()
# # tau_max = 2
# # link_assumptions = {}
# # for j in range(4):
# #     if j in [0, 1, 2]:
# #         # Directed lagged links
# #         link_assumptions[j] = {(var, -lag): '-?>' for var in [0, 1, 2]
# #                                for lag in range(1, tau_max + 1)}
# #         # Unoriented contemporaneous links
# #         link_assumptions[j].update({(var, 0): 'o?o' for var in [0, 1, 2] if var!= j})
# #         # Directed lagged and contemporaneous links from the sun (3)
# #         link_assumptions[j].update({(var, -lag): '-?>' for var in [3]
# #                                     for lag in range(0, tau_max + 1)})
# #     else:
# #         link_assumptions[j] = {}
# #
# # all_results = []
# # for i, segment in enumerate(segments):
# #     segment_df = pp.DataFrame(segment, var_names=[r'X0', r'X1', r'X2', 'Sun'])
# #     pcmci_parcorr = PCMCI(
# #         dataframe=segment_df,
# #         cond_ind_test=ci_test,
# #         verbosity=1)
# #     print(f"Results for segment {i + 1}:")
# #     results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
# #                                           pc_alpha=0.001,
# #                                           reset_lagged_links=False,
# #                                           link_assumptions=link_assumptions)
#
#
#
# # # PCMCIPLUS
# # ci_test = CMIknn(significance="fixed_thres", verbosity=3)
# # ci_test = ParCorr()
# # tau_max = 2
# # link_assumptions = {} #{}
# # for j in range(4):
# #     if j in [0, 1, 2]:
# #         # Directed lagged links
# #         link_assumptions[j] = {(var, -lag): '-?>' for var in [0, 1, 2]
# #                          for lag in range(1, tau_max + 1)}
# #         # Unoriented contemporaneous links
# #         link_assumptions[j].update({(var, 0): 'o?o' for var in [0, 1, 2] if var != j})
# #         # Directed lagged and contemporaneous links from the sun (3)
# #         link_assumptions[j].update({(var, -lag): '-?>' for var in [3]
# #                          for lag in range(0, tau_max + 1)})
# #     else:
# #         link_assumptions[j] = {}
# #
# # for j in link_assumptions:
# #     print(link_assumptions[j])
# # pcmci_parcorr = PCMCI(
# #     dataframe=dataframe,
# #     cond_ind_test=ci_test,
# #     verbosity=1)
# # results = pcmci_parcorr.run_pcmciplus(tau_max=tau_max,
# #                 pc_alpha=0.001,
# #                 reset_lagged_links=False,
# #                 link_assumptions=link_assumptions
# #                 ) #, alpha_level = 0.01)
# # print(results['graph'].shape)
# # # print(results['graph'][:,3,:])
# # print(np.round(results['p_matrix'][:,:,0], 2))
# # print(np.round(results['val_matrix'][:,:,0], 2))
# # print(results['graph'][:,:,0])
# #
# #
# # algo_pelt= Pelt(model="l2",jump = 255)#设置模型为L1范数，最小突变点长度为2
# # #检测变点
# # bkps_pelt = algo_pelt.fit(first_column_data)
# # breakpoints =algo_pelt.predict(pen=0.01)
# # print("Using PELT algorithm:")
# # print("Detected breakpoints using L2 norm:", breakpoints)
#
#
# #
# # T = 1000
# # data = random_state.standard_normal((T, 4))
# # # Simple sun
# # data[:, 3] = random_state.standard_normal((T))  # np.sin(np.arange(T)*20/np.pi) + 0.1*random_state.standard_normal((T))
# # c = 0.8
# # for t in range(1, T):
# #     data[t, 0] += 0.4 * data[t - 1, 0] + 0.4 * data[t - 1, 1] + c * data[t - 1, 3]
# #     data[t, 1] += 0.5 * data[t - 1, 1] + c * data[t, 3]
# #     data[t, 2] += 0.6 * data[t - 1, 2] + 0.3 * data[t - 2, 1]  # + c*data[t-1,3]
# # dataframe = pp.DataFrame(data, var_names=[r'X^0', r'X^1', r'X^2', 'Sun'])
# # tp.plot_timeseries(dataframe); plt.show()
#
# from __future__ import print_function

#
#
#
#
# random_state = np.random.default_rng(seed=43)
# def lin_f(x): return x
# def nonlin_f(x): return (x + 5. * x ** 2 * np.exp(-x ** 2 / 20.))
# # T = 1500
# # data = np.random.default_rng().normal(2, 1, size = (T, 4))
# # # data = random_state.standard_normal((T, 4))
# #
# # c = 0.8
# # for t in range(0, 500):
# #     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
# #     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
# #     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
# #     data[t, 3] += 0.8 * data[t-1, 3]
# # for t in range(500, T):
# #     data[t, 0] += 0.8 * data[t-1, 0] + 0.2 * data[t-1, 2]
# #     data[t, 1] += 0.5 * data[t-1, 1]
# #     data[t, 2] += 0.3 * data[t-1, 2] + 0.4 * data[t-2, 1]
# #     data[t, 3] += 0.3 * data[t-1, 3]
# # for t in range(1000, T):
# #     data[t, 0] += 0.4 * data[t-1, 0] + 0.4 * data[t-1, 1] + c * data[t-1, 3]
# #     data[t, 1] += 0.5 * data[t-1, 1] + c * data[t, 3]
# #     data[t, 2] += 0.6 * data[t-1, 2] + 0.3 * data[t, 1]
# #     data[t, 3] += 0.8 * data[t-1, 3]
# # dataframe = pp.DataFrame(data, var_names=[r'X0', r'X1', r'X2', 'Sun'])
# # print(data)
# # first_column_data = data[:, 0]
# # first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# # # 将秩次转换为从1开始的连续整数
# # first_column_data_ranks = first_column_data_ranks + 1
# # # 打印秩次
# # print(first_column_data_ranks)
# # # detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# # # result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
# # # result_mosum = mosum.multiscale_bottomUp(values_array_ranks, G = [3,6,9], eta = 2/3, alpha=0.001)
# # result_mosum = mosum.multiscale_localPrune(values_array_ranks, G = [10,15,20], criterion="epsilon", epsilon = 2/3, alpha=0.01)
# # # summary and print methods
# # result_mosum.summary()
# # result_mosum.print()
# # # plot the output
# # result_mosum.plot()
# # plt.show()
#
#
#
# # 你的数值列表
# values = [
#     0.548054828851267, 0.108065477446078, 0.159644361671702, 0.255741068991730, 0.0236786239653456,
#     0.225293637657200, 0.360152126058572, 0.302793450610223, 0.495199084785648, 0.361699731967659,
#     0.288746451650207, 0.311966515786077, 0.398066072479398, 0.190100144613410, 0.0343870104500,
#     0.120888787377312, 0.381299414529124, 0.374678819599434, 0.241816969401687, -0.0172256761064497,
#     0.151564157366174, 0.0792363215656019, -0.0488441740542104, 0.0422466575491896, -0.141851377132136,
#     0.0990236260736834, 0.0312774703118287, -0.0573892384492988, 0.244326531982757,
#     -0.338361947483559, -0.0518875645446012, 0.148484931644656, -0.205625778374573,
#     -0.132200743912698, -0.105698161317981, 0.178318237822030, 0.0762100896878891,
#     0.292863909978405, 0.111439982845626, 0.363857552863368, 0.259801937872358,
#     0.0729997479668308, 0.361578638163063, 0.221927100919805, 0.593398591603231,
#     0.153229764729903, 0.424412325333526, 0.391922534444257, 0.199037897803507,
#     0.353328093416055, 0.181412232448095, -0.0173180679504521, 0.194648522547737,
#     0.362454238070246, 0.118791754506458, -0.114797127570850, 0.129140845292384,
#     0.0713369042573362, 0.141937344050961, -0.0212332958217078, 0.175484041340127,
#     0.0303057626460393, 0.326388097931126, 0.434145710333401, 0.676943470676115,
#     0.764165017618700, 0.745193772541677, 0.742437947678138, 0.863379249067743,
#     0.641559762191089, 0.537269528284888, 0.818959735122729, 0.802404930901335,
#     0.567150028607477, 0.458062137367003, 0.356093469316563, 0.342062991183548,
#     -0.145548368314035, 0.338407940840552, -0.189641146604775, 0.164527185933008,
#     0.236798423799697, 0.0612897884698173, 0.0198912376490328, 0.317866064235381,
#     0.547434741908840, 0.126175577277498, 0.177386779661916, 0.286445142033722,
#     0.245286714322088, 0.210553742380902, 0.261395707143683, 0.190035125511402,
#     0.00429819546956469, 0.287062935616955, -0.0196296300740212, -0.175108889273972,
#     0.00415029306805184, 0.00905924971815275, 0.159231742976046, -0.0349669308678230,
#     0.185055410855722, -0.107094414260001, 0.238686181128347, 0.127796976293462,
#     0.287899402513304, 0.222471697840347, 0.118100832346589, -0.0499735552227914,
#     -0.117977754233259, -0.251194548114207, -0.221601150443199, -0.10563311396510,
#     0.0284672827743480, -0.277953735172174, 0.241655735231637, 0.417487921485027,
#     0.405447520960845, 0.371941030113109, 0.419741085276209, 0.984872581301156,
#     0.370853506201600, 0.620593604977158, 0.426831039099780, 0.400297190171936,
#     0.569835376834383, 0.146452942205907, 0.194572793837573, -0.0668090821430109,
#     0.177693370607884, 0.531436817348630, 0.238008235863703, 0.406571738662761,
#     0.287391557090998, 0.149962863270810, 0.286041155357509, 0.463829378394337,
#     0.326539735624632, 0.198119920305, 0.202865067866757, 0.424533758140583,
#     0.0846531400180585, 0.345620612981813, 0.279262940617219, 0.361775861022818,
#     0.239942738120075, 0.310776860230718, 0.331381616781768, 0.0280069450000969,
#     -0.0904479990722788, -0.132652751412106, -0.0767028444365349, -0.586837300975076,
#     -0.142927271036645, -0.513036396015626, -0.30018665632012, -0.540986058593683,
#     -0.285015049531817, -0.308043222956777, -0.189049403124973, -0.254871984654968,
#     -0.357914253451910, -0.217031561275893, 0.0551229579445153, -0.135593626874090
# ]
#
# # 将列表转换为NumPy数组
# values_array = np.array(values)
# print(values_array)
# values_array_ranks = np.argsort(np.argsort(values_array))
# # 将秩次转换为从1开始的连续整数
# values_array_ranks = values_array_ranks+1
# # 打印秩次
# print(values_array_ranks)
#
#
# # detect changes，使用提取出来的第一列数据进行变点检测，这里带宽G等参数你可以根据实际情况调整
# # result_mosum = mosum.mosum(first_column_data_ranks, G=400, criterion="eta", eta = 2/3, alpha=0.1, boundary_extension=True)
# result_mosum = mosum.multiscale_bottomUp(values_array_ranks, G = [3,6,9], eta = 2/3, alpha=0.001)
# result_mosum = mosum.multiscale_localPrune(values_array_ranks, G = [10,15,20], criterion="epsilon", epsilon = 2/3, alpha=0.01)
# # summary and print methods
# result_mosum.summary()
# result_mosum.print()
#
# # plot the output
# result_mosum.plot()
# plt.show()
















# # Determine a fixed y-axis range that encompasses all data
# # Create a figure with subplots for each variable
# fig, axs = plt.subplots(4, 1, figsize=(10, 8), sharex=False)  # Changed sharex to False
# # Plot time series with different colors for each variable
# for i, ax in enumerate(axs.flat):  # 使用flat属性来迭代每个子图
#     # 计算当前变量的最小值和最大值
#     y_min = np.min(data[:, i])
#     y_max = np.max(data[:, i])
#     # 绘制当前变量的时间序列
#     ax.set_ylim([y_min, y_max])  # 为当前子图设置y轴范围
# for i, color in enumerate(colors):
#     axs[i].plot(data[:, i], color=color, label=r'$X^{}'.format(i))
#     axs[i].set_ylabel(dataframe.var_names[i], fontname='Times New Roman')
# # Hide x-axis ticks, labels, and spines for the first three subplots
# for ax in axs[:-1]:
#     ax.set_xticks([])
#     ax.set_xticklabels([])
#     ax.spines['bottom'].set_visible(False)
# # Set xticks, labels, and spine for the last subplot
# axs[-1].set_xticks(np.arange(0, 3312, step=500))
# axs[-1].set_xticklabels([str(int(tick)) for tick in np.arange(0, 3312, step=500)], fontname='Times New Roman')
# axs[-1].set_xlabel('Time', fontname='Times New Roman')
# axs[-1].spines['bottom'].set_visible(True)
# # Adjust layout to prevent overlap
# plt.tight_layout()
# # Set the spines to only show the left frame for all subplots
# for ax in axs:
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
# plt.show()



