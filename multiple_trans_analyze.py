
from __future__ import print_function
import sys
import re
import commands
import os
import time


def help_info():
	print("===================================================")
	print("                  MultiTrans Usage                 ")
	print("===================================================")
	print("================Analyze================")
	print("** Example for extract common transcripts: ")
	print('     python multiple_trans_analyze.py -t trans1.fa,trans2.fa,trans3.fa -o /out_dir/')
	print("---------------------------------------------------")
	print("     --transcripts_file_list (or -t): list of transcriptome.fa separated by commas")
	print("     --out_dir (or -o) : directory of output")
	print("===================================================")
#####################################################################################################


trans_file_list = []
merge_list = []
trans_str = ''
merge_str = ''
out_dir = ''

for i in range(len(sys.argv)):	
	if sys.argv[i] == "--help" or sys.argv[i] == '-h':
		help_info()
		sys.exit(-1)	
	if sys.argv[i] == "--transcripts_file_list" or sys.argv[i] == '-t':
		trans_str = sys.argv[i+1]
	if sys.argv[i] == "--merge" or sys.argv[i] == '-m':
		merge_str = sys.argv[i+1]
	if sys.argv[i] == "--out_dir" or sys.argv[i] == '-o':
		out_dir = sys.argv[i+1]

trans_file_list = trans_str.split(',')

if len(merge_str) > 0:
    merge_list = merge_str.split(',')
    new_trans = []
    for i in range(len(merge_list)):
        idx = int(merge_list[i])-1
        new_trans = new_trans + [trans_file_list[idx]]
    trans_file_list = new_trans


trans_vec = []
trans_name = []
current_tran = ''
for lines in open(trans_file_list[0]):
    if lines[0] == '>':
        if len(current_tran) > 0:
            trans_vec = trans_vec + [current_tran]
        trans_name = trans_name + [lines[1:len(lines)-1]]
        current_tran = ''
    else:
        current_tran = current_tran + lines[0:len(lines)-1]
trans_vec = trans_vec + [current_tran]
for i in range(len(trans_file_list)-1):
    temp_trans = []
    temp_name = []
    current_name = ''
    current_tran = ''
    for lines in open(trans_file_list[i+1]):
        if lines[0] == '>':
            if len(current_tran) > 0:
                if current_tran in trans_vec:
                    idx = trans_vec.index(current_tran)
                    temp_trans = temp_trans + [current_tran]
                    new_name = trans_name[idx] + ',' + current_name
                    temp_name = temp_name + [new_name]
            current_name = lines[1:len(lines)-1]
            current_tran = ''
        else:
            current_tran = current_tran + lines[0:len(lines)-1]
    if current_tran in trans_vec:
        idx = trans_vec.index(current_tran)
        temp_trans = temp_trans + [current_tran]
        new_name = trans_name[idx] + ',' + current_name
        temp_name = temp_name + [new_name]
    trans_vec = temp_trans
    trans_name = temp_name

out_name = out_dir + 'common_transcripts.fa'
out = open(out_name, 'w')
for i in range(len(trans_vec)):
    print('>'+trans_name[i]+'\n',file=out,end='')
    print(trans_vec[i]+'\n',file=out,end='')

