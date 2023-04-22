
#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import os
import time


def help_info():
	print("===================================================")
	print("                  MultiTrans Usage                 ")
	print("===================================================")
	print("** Example for assembly graph: ")
	#print("---------------------------------------------------")
	print('     python MultiTrans.py -m 1 -g assembly_graph_with_scaffolds.gfa -o /out_dir/ -r1 reads_1.fastq -r2 reads_2.fastq')
	print("---------------------------------------------------")
	print("** Example for splicing graph: ")
	print('     python MultiTrans.py -m 2 -g /RawGraphs_/ -o /out_dir/')
	print("----------------------------------------------------")
	print("** Parameters")
	#print("===================================================")
	print("     --test (or -t): testing the software")
	print("     --model (or -m) <int>: 1 assembly graph, 2 splcing graph")
	print("     --graph_file (or -g) <string>: *.gfa file of assembly graph")
	print("                                or  directory of splicing graph")
	print("     --out_dir (or -o): directory of output")
	print("     --reads1 (or -r1) <string>: paired-end reads")
	print("     --reads2 (or -r2) <string>: paired-end reads")
	print("     --help (or -h): MultiTrans usage")
	print("===================================================")
#####################################################################################################

model='' # 1 assembly graph, 2 splicing graph
u = 0.1
graph_file = ''
graph_sum = 0
out_dir = ''
reads_file1 = ''
reads_file2 = ''
is_test = 0

for i in range(len(sys.argv)):
	if sys.argv[i] == "--model" or sys.argv[i] == "-m":
		model = sys.argv[i+1]	
	if sys.argv[i] == "--help" or sys.argv[i] == '-h':
		help_info()
		sys.exit(-1)		
	if sys.argv[i] == "--graph_file" or sys.argv[i] == '-g':
		graph_file = sys.argv[i+1]
	if sys.argv[i] == "--graph_sum" or sys.argv[i] == '-s':
		graph_sum = int(sys.argv[i+1])
	if sys.argv[i] == "--out_dir" or sys.argv[i] == '-o':
		out_dir = sys.argv[i+1]
	if sys.argv[i] == "--reads1" or sys.argv[i] == '-r1':
		reads_file1 = sys.argv[i+1]
	if sys.argv[i] == "--reads2" or sys.argv[i] == '-r2':
		reads_file2 = sys.argv[i+1]
	if sys.argv[i] == "--test" or sys.argv[i] == '-t':
		is_test = 1	
	if sys.argv[i] == "-u":
		u = float(sys.argv[i+1])


if is_test == 0:
	if graph_file=="" or model == '':
		help_info()
		sys.exit(-1)

	if model != '1' and model != '2':
		help_info()
		sys.exit(-1)
	if model == '2' and graph_file[len(graph_file)-1] != '/':
		help_info()
		sys.exit(-1)
	if out_dir == '':
		out_dir = os.getcwd()
		out_dir = out_dir + '/'
MultiTrans_dir = sys.path[0]

if is_test == 1:
	out_dir = MultiTrans_dir + '/test/'
	command = "mkdir " + out_dir + "MultiTrans_Out_Dir"
	os.system(command)
	out_dir = out_dir + "MultiTrans_Out_Dir/"
	command = "mkdir " + out_dir + "SG"
	os.system(command)  
	graph_file = MultiTrans_dir + '/test/RawGraphs_/'
	if graph_sum == 0:
		for each_file in os.listdir(graph_file):
			if each_file[len(each_file)-9:len(each_file)] == '.rg.debug':
					graph_sum = graph_sum + 1
	command = MultiTrans_dir+'/splicing_graph -g ' + str(graph_sum) + ' -i ' + graph_file + ' -o ' + out_dir
	os.system(command) 
	print('solving mixed integer linear programming...')
	command = 'python '+MultiTrans_dir+'/write_program_splicing.py -in ' + out_dir + 'graph.info -u' + str(u)
	os.system(command) 
	print("Success!")
 	graph_file = MultiTrans_dir + '/test/assembly_graph_with_scaffolds.gfa'
 	reads_file1 =  MultiTrans_dir + '/test/reads_1.fastq'
	reads_file2 =  MultiTrans_dir + '/test/reads_2.fastq'
	command = MultiTrans_dir+'/node_seq -g ' + graph_file + ' -o ' + out_dir
	os.system(command) 
	print('algning reads to non-trivial graphs...')
	command = 'bwa index ' + out_dir + 'nodeSeq.fa'
	os.system(command) 
	command = 'bwa mem ' + out_dir + 'nodeSeq.fa ' + reads_file1 + ' '+ reads_file2 + ' >'+out_dir + 'aln_node.sam' 
	os.system(command) 
	command = MultiTrans_dir+'/assembly_graph -g ' + graph_file + ' -m ' + out_dir + 'aln_node.sam' + ' -o ' + out_dir
	os.system(command) 
	print('solving mixed integer linear programming...')
	command = 'python '+MultiTrans_dir+'/write_program_assembly.py -in ' + out_dir + 'graph.info -u' + str(u) 
	os.system(command) 
	print("\n====== MultiTrans pipeline finished.\n")
	print("The assembled transcripts can be found here: "+out_dir+"MultiTrans.fa\n")
	print("Thank you for using MultiTrans!\n")
	sys.exit(-1)

command = "mkdir " + out_dir + "MultiTrans_Out_Dir"
os.system(command)
out_dir = out_dir + "MultiTrans_Out_Dir/"
command = "mkdir " + out_dir + "SG"
os.system(command)  


if model == '1': #assembly graph
	
	if reads_file1 != '' and reads_file2 != '':
		command = MultiTrans_dir+'/node_seq -g ' + graph_file + ' -o ' + out_dir
		os.system(command) 

		print('algning reads to non-trivial graphs...')
		command = 'bwa index ' + out_dir + 'nodeSeq.fa'
		os.system(command) 

		command = 'bwa mem ' + out_dir + 'nodeSeq.fa ' + reads_file1 + ' '+ reads_file2 + ' >'+out_dir + 'aln_node.sam' 
		os.system(command) 

		command = MultiTrans_dir+'/assembly_graph -g ' + graph_file + ' -m ' + out_dir + 'aln_node.sam' + ' -o ' + out_dir
		os.system(command) 
	else:
		command = MultiTrans_dir+'/assembly_graph -g ' + graph_file + ' -o ' + out_dir
		os.system(command) 
		
	print('solving mixed integer linear programming...')
	command = 'python '+MultiTrans_dir+'/write_program_assembly.py -in ' + out_dir + 'graph.info -u' + str(u) 
	os.system(command)  

if model == '2': #splicing graph

	
	if graph_sum == 0:
		for each_file in os.listdir(graph_file):
			if each_file[len(each_file)-9:len(each_file)] == '.rg.debug':
					graph_sum = graph_sum + 1
	command = MultiTrans_dir+'/splicing_graph -g ' + str(graph_sum) + ' -i ' + graph_file + ' -o ' + out_dir
	os.system(command) 
	print('solving mixed integer linear programming...')
	command = 'python '+MultiTrans_dir+'/write_program_splicing.py -in ' + out_dir + 'graph.info -u' + str(u)
	os.system(command)  

print("\n====== MultiTrans pipeline finished.\n")
print("The assembled transcripts can be found here: "+out_dir+"MultiTrans.fa\n")
print("Thank you for using MultiTrans!\n")







