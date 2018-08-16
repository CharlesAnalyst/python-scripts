#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy

import os
import pandas as pd
import numpy as np

######################################################################
from scipy.stats import norm

data_dir = "/data3/xs/tissue_m6a/2018.1/fig2/fig2_315/1to8_exprCV"
data_1 = "%s/5utr_cds_stop_3utr_expr_cv.txt" % data_dir
data_2 = "%s/mrna_lincRNA_expr_cv.txt" % data_dir
data_3 = "/data3/xs/tissue_m6a/2018.1/fig2/fig2_315/1to8_exprCV/m6A_modified_1to8.txt"
result_dir = "/data5/galaxy/project/trend_test"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)
###############################################################


def main():
    for in_file in [data_1, data_2, data_3]:
        process_data(in_file)


def process_data(in_file):
    df = pd.read_table(in_file, sep="\t")
    df_sort = df.stack().sort_index(level=1)
    level_2_indexes = list(set(df_sort.index.get_level_values(level=1)))
    #
    index_dict = {}
    for x in level_2_indexes:
        feature_type = str(x).split("-")[0]
        index_dict[feature_type] = index_dict.get(feature_type, []) + [x]
    #
    for feature_type, index_list in index_dict.items():
        df_list = []
        index_list = sorted(index_list)
        for i_index in index_list:
            print(i_index)
            i_df = pd.DataFrame(df_sort[:, i_index], columns=["cv"])
            i_df["category"] = i_index
            df_list.append(i_df)
        df_com = pd.concat(df_list)
        result_file = os.path.join(result_dir, "%s_data.txt" % feature_type)
        df_com.to_csv(result_file, sep="\t", header=False, index=False)


def mk_test(x, alpha=0.05):
    """
    Input:
    x: a vector of data
    alpha: significance level(0.05, default)

    Output:
    trend: tells the trend(increasing, decreasing or no trend)
    h: True( if trend is present) or False( if trend is absence)
    p: pvalue of the significance test
    z: normalized test statistics
"""
    # x = np.random.rand(100)
    # trend, h, p, z = mk_test(x, 0.05)
    n = len(x)

    # calculate S
    s = 0
    for k in range(n - 1):
        for j in range(k + 1, n):
            s += np.sign(x[j] - x[k])

    # calculate the unique data
    unique_x = np.unique(x)
    g = len(unique_x)

    # calculate the var(s)
    if n == g:  # there is no tie
        var_s = (n * (n - 1) * (2 * n + 5)) / 18
    else:  # there are some ties in data
        tp = np.zeros(unique_x.shape)
        for i in range(len(unique_x)):
            tp[i] = sum(unique_x[i] == x)
        var_s = (n * (n - 1) * (2 * n + 5) + np.sum(tp * (tp - 1) * (2 * tp + 5))) / 18

    if s > 0:
        z = (s - 1) / np.sqrt(var_s)
    elif s == 0:
        z = 0
    elif s < 0:
        z = (s + 1) / np.sqrt(var_s)

    # calculate the p_value
    p = 2 * (1 - norm.cdf(abs(z)))  # two tail test
    h = abs(z) > norm.ppf(1 - alpha / 2)

    if (z < 0) and h:
        trend = 'decreasing'
    elif (z > 0) and h:
        trend = 'increasing'
    else:
        trend = 'no trend'

    return trend, h, p, z


if __name__ == '__main__':
    main()
