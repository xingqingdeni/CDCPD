import numpy as np
def analyze_causal_relations(true_stages_causal, tested_causal):
    # 用于存储每个阶段处理后的检测因果关系字典
    all_stages_causal = []
    # 按阶段分割检测因果关系文本
    stages_text = tested_causal.split("###")
    stages_text = [s.strip() for s in stages_text if s.strip()]

    for stage_text in stages_text:
        stage_causal = {}
        lines = stage_text.splitlines()
        current_variable = None
        for line in lines:
            line = line.strip()
            if line.startswith("Variable"):
                current_variable = line.split()[1]
                stage_causal[current_variable] = []
            elif line.startswith("(") and current_variable:
                parts = line[1:-1].split(":")[0].split()
                var_name = parts[0]
                time_index = parts[1].strip(")")
                if time_index == "0":
                    time_index = "[t]"
                else:
                    time_index = f"[t - {abs(int(time_index))}]"
                stage_causal[current_variable].append(f"{var_name}{time_index}")
        all_stages_causal.append(stage_causal)

    # 统计每一阶段真实因果关系数量、检测因果关系各阶段的因果关系数量
    true_counts = [sum(len(causes) for _, causes in stage.items()) for stage in true_stages_causal]
    tested_counts = [sum(len(causes) for _, causes in stage.items()) for stage in all_stages_causal]

    # 统计每一阶段真实因果关系与其他各阶段检测因果关系相同的因果关系数量
    num_true_stages = len(true_stages_causal)
    num_tested_stages = len(all_stages_causal)
    same_counts_matrix = [[0] * num_tested_stages for _ in range(num_true_stages)]
    for true_stage_index in range(num_true_stages):
        true_stage = true_stages_causal[true_stage_index]
        for tested_stage_index in range(num_tested_stages):
            tested_stage = all_stages_causal[tested_stage_index]
            same_count = 0
            for var in true_stage:
                if var in tested_stage:
                    true_edges = set(true_stage[var])
                    tested_edges = set(tested_stage[var])
                    same_count += len(true_edges & tested_edges)
            same_counts_matrix[true_stage_index][tested_stage_index] = same_count

    # 打印每一阶段的数量统计信息
    print("真实因果关系各阶段的因果关系数量:")
    for i in range(num_true_stages):
        print(f"第{i + 1}阶段: {true_counts[i]}")
    print()

    print("检测因果关系各阶段的因果关系数量:")
    for i in range(num_tested_stages):
        print(f"第{i + 1}阶段: {tested_counts[i]}")
    print()

    print("各阶段真实因果关系与各阶段检测因果关系相同的数量:")
    for true_stage_index in range(num_true_stages):
        print(f"第{true_stage_index + 1}阶段真实因果关系与各阶段检测因果关系相同的数量：")
        for tested_stage_index in range(num_tested_stages):
            print(
                f"  与第{tested_stage_index + 1}阶段检测因果关系相同数量: {same_counts_matrix[true_stage_index][tested_stage_index]}")
        print()

    return true_counts, tested_counts, same_counts_matrix


def calculate_precision(cpts, true_cpts, sum, same_counts_matrix, tested_counts):
    if cpts[0] > true_cpts[0] and cpts[1] > true_cpts[1]:
        precision = (((((true_cpts[0] / sum) * (same_counts_matrix[0][0] / tested_counts[0])
                        + (abs(true_cpts[0] - cpts[0]) / sum) * (same_counts_matrix[1][0] / tested_counts[0]))
                       + (abs(true_cpts[1] - cpts[0]) / sum) * (same_counts_matrix[1][1] / tested_counts[1]))
                      + (abs(cpts[1] - true_cpts[1]) / sum) * (same_counts_matrix[2][1] / tested_counts[1]))
                     + (abs(true_cpts[2] - cpts[1]) / sum) * (same_counts_matrix[2][2] / tested_counts[2]))
        return precision
    elif cpts[0] > true_cpts[0] and cpts[1] < true_cpts[1]:
        precision = (((((true_cpts[0] / sum) * (same_counts_matrix[0][0] / tested_counts[0])
                        + (abs(true_cpts[0] - cpts[0]) / sum) * (same_counts_matrix[0][0] / tested_counts[0]))
                       + (abs(cpts[0] - cpts[1]) / sum) * (same_counts_matrix[1][1] / tested_counts[1]))
                      + (abs(cpts[1] - true_cpts[1]) / sum) * (same_counts_matrix[1][2] / tested_counts[2]))
                     + (abs(true_cpts[2] - true_cpts[1]) / sum) * (same_counts_matrix[2][2] / tested_counts[2]))
        return precision
    elif cpts[0] < true_cpts[0] and cpts[1] > true_cpts[1]:
        precision = (((((cpts[0] / sum) * (same_counts_matrix[0][0] / tested_counts[0])
                        + (abs(true_cpts[0] - cpts[0]) / sum) * (same_counts_matrix[0][1] / tested_counts[1]))
                       + (abs(true_cpts[1] - true_cpts[0]) / sum) * (same_counts_matrix[1][1] / tested_counts[1]))
                      + (abs(cpts[1] - true_cpts[1]) / sum) * (same_counts_matrix[2][1] / tested_counts[1]))
                     + (abs(true_cpts[2] - cpts[1]) / sum) * (same_counts_matrix[2][2] / tested_counts[2]))
        return precision
    elif cpts[0] < true_cpts[0] and cpts[1] < true_cpts[1]:
        precision = (((((cpts[0] / sum) * (same_counts_matrix[0][0] / tested_counts[0])
                        + (abs(true_cpts[0] - cpts[0]) / sum) * (same_counts_matrix[0][1] / tested_counts[1]))
                       + (abs(true_cpts[0] - cpts[1]) / sum) * (same_counts_matrix[1][1] / tested_counts[1]))
                      + (abs(cpts[1] - true_cpts[1]) / sum) * (same_counts_matrix[1][2] / tested_counts[2]))
                     + (abs(true_cpts[2] - true_cpts[1]) / sum) * (same_counts_matrix[2][2] / tested_counts[2]))
        return precision
    else:
        print("不符合计算precision的条件")
        return None

def convert_links(links):
    stage_causal = {}
    for key, dependencies in links.items():
        var_name = f"X{key}"
        causal_list = []
        for (var_index, delay), _, _ in dependencies:
            causal_list.append(f"X{var_index}[t {'-' if delay < 0 else '+'} {abs(delay)}]" if delay != 0 else f"X{var_index}[t]")
        stage_causal[var_name] = causal_list
    return stage_causal

def wbs(data, M=100, n_bkps=2, min_size=30):
    n = len(data)
    intervals = []
    cps = []

    # Step 1: Sample M random intervals
    for _ in range(M):
        start = np.random.randint(0, n - min_size)
        end = np.random.randint(start + min_size, n)
        intervals.append((start, end))

    # Step 2: For each interval, compute max CUSUM statistic
    stats = []
    for (s, e) in intervals:
        best_t, best_val = -1, 0
        for t in range(s + 1, e):
            mu1 = np.mean(data[s:t])
            mu2 = np.mean(data[t:e])
            val = np.abs(mu1 - mu2)
            if val > best_val:
                best_val = val
                best_t = t
        if best_t != -1:
            stats.append((best_val, best_t))

    # Step 3: Select top change points
    stats.sort(reverse=True)
    cps = sorted(list(set([t for _, t in stats[:n_bkps]])))

    return cps
