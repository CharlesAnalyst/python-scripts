#!/usr/bin/python
# -*-coding:utf-8-*-
# @author galaxy


import os
import glob


data_dir = "/data5/galaxy/project/data/input_data/total/input_bed"
result_dir = "/data5/galaxy/project/data/input_data/autosomal/input_bed"
if not os.path.exists(result_dir):
    os.makedirs(result_dir)


def extract_autosomal(in_file):
    result_file = os.path.join(result_dir, os.path.basename(in_file))
    # os.system("awk '$1 ~ /^chr(1?[0-9]|2[0-2]|X|Y)$/' %s > %s" % (in_file, result_file))
    os.system("awk '$1 ~ /^chr(1?[0-9]|2[0-2])$/' %s > %s" % (in_file, result_file))


if __name__ == '__main__':
    file_list = glob.glob("%s/*_discrete.bed" % data_dir)
    for i_file in file_list:
        extract_autosomal(i_file)
