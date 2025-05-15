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

file_path = r'C:\Users\10362\Desktop\论文\数据集\EEG_Eyes.arff.gz.txt'
# 使用NumPy的loadtxt函数读取数据
data = np.loadtxt(file_path, delimiter=',')  # 假设数据由逗号分隔
data = data[3342:6653]
# last_column = data[:, -1]
data = data[:, :-1]
dataframe = pp.DataFrame(data, var_names=[r'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'])
tp.plot_timeseries(dataframe);
# plt.show()
first_column_data = data[:,6]
first_column_data = (first_column_data - np.min(first_column_data)) / (np.max(first_column_data) - np.min(first_column_data))
first_column_data_ranks = np.argsort(np.argsort(first_column_data))
# 将秩次转换为从1开始的连续整数
first_column_data_ranks = first_column_data_ranks
# 打印秩次
print(first_column_data_ranks)
algo_pelt= Pelt(model="l2",min_size=523,jump = 255)#设置模型为L1范数，最小突变点长度为2
#检测变点
bkps_pelt = algo_pelt.fit(first_column_data)
breakpoints =algo_pelt.predict(pen=10)
print("Using PELT algorithm:")
print("Detected breakpoints using L2 norm:", breakpoints)
# result_mosum = mosum.mosum(first_column_data_ranks, var_est_method= 'mosum', G=1000, criterion="eta", eta = 0.9, alpha=0.01, boundary_extension=True)
result_mosum = multiscale_struct_bottomUp(first_column_data_ranks, G = [230,300,400], eta = 2/3, alpha=0.01, distance = 600, min_jump = 0.65)
# result_mosum = mosum.multiscale_localPrune(first_column_data_ranks, G = [100], criterion="eta", eta = 0.99, alpha=0.001)
cpts = result_mosum.cpts.tolist()
result_mosum.print()

# plot the output
true_cpts = [1012, 1904, 2586, 3311]
bkps = true_cpts
result = cpts
# rpt.display(data, bkps, result)
# plt.xlim(0, 3312)
# plt.show()
sum = 3311
cpts.append(sum)
segments = []
start = 0
for cpt in cpts:
    segment = data[start:cpt, :]
    segments.append(segment)
    start = cpt
# PCMCIPLUS
ci_test = CMIknn(significance="fixed_thres", verbosity=1)
ci_test = ParCorr()
for i, segment in enumerate(segments):
    segment_df = pp.DataFrame(segment, var_names=[r'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'])
    pcmci = PCMCI(
        dataframe=segment_df,
        cond_ind_test=ci_test,
        verbosity=1)
    print(f"Results for segment {i + 1}:")
    results = pcmci.run_pcmci(tau_min=1,tau_max=1,pc_alpha=0.65,alpha_level=0.0001)
    tp.plot_graph(
        val_matrix=results['val_matrix'],
        graph=results['graph'],
        var_names=[r'AF3', 'F7', 'F3', 'FC5', 'T7', 'P7', 'O1', 'O2', 'P8', 'T8', 'FC6', 'F4', 'F8', 'AF4'],
        link_colorbar_label='MCI',
        arrow_linewidth=3.0,
        curved_radius=0.3,
        );
    plt.show()
#CDNOD
c_indx = np.arange(1, sum+1).reshape(-1, 1) # 将指定列索引对应的列数据提取出来，形成二维数组，列数为1
# 使用默认参数构建因果图
cg_default = cdnod(data, c_indx)
# 可视化默认参数构建的因果图（展示）
cg_default.draw_pydot_graph()














