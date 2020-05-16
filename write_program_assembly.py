
#!/usr/bin/env python

from __future__ import print_function
import sys
import re
import os
import time



def revcmps(sequence):
	cmp_seq = ''
	for i in range(len(sequence)):
		c = sequence[i]
		revc = 'N'
		if c == 'G':
			revc = 'C'
		if c == 'g':
			revc = 'C'
		if c == 'C':
			revc = 'G'
		if c == 'c':
			revc = 'G'
		if c == 'A':
			revc = 'T'
		if c == 'a':
			revc = 'T'
		if c == 'T':
			revc = 'A'
		if c == 't':
			revc = 'A'
		cmp_seq = revc + cmp_seq
	return cmp_seq


def load_par(par_file):
	
	C = [] # adj_maxtrix
	infile = open(par_file)
	line = infile.readline()	
	if line[len(line)-1] == '\n':
		line = line[0:len(line)-1]	
	lines = line.split(',;')
	for i in range(len(lines)):
		if len(lines[i]) < 1:
			continue
		x = lines[i].split(',')
		cc = []
		for j in range(len(x)):
			cc = cc + [float(x[j])]
		C = C + [cc]

	A = [] # path node
	line = infile.readline()
	if line[len(line)-1] == '\n':
		line = line[0:len(line)-1]	
	lines = line.split(',;')
	for i in range(len(lines)):
		if len(lines[i]) < 1:
			continue
		x = lines[i].split(',')
		cc = []
		for j in range(len(x)):
			cc = cc + [int(x[j])]
		A = A + [cc]

	P = [] # path 
	line = infile.readline()
	if line[len(line)-1] == '\n':
		line = line[0:len(line)-1]		
	lines = line.split(',;')
	for i in range(len(lines)):
		if lines[i] == '' or lines[i] == '\n':
			continue
		x = lines[i].split(',')
		cc = []
		for j in range(len(x)):
			cc = cc + [int(x[j])]
		P = P + [cc]

	S = [] # path stand
	line = infile.readline()
	if line[len(line)-1] == '\n':
		line = line[0:len(line)-1]		
	lines = line.split(',;')
	for i in range(len(lines)):
		if lines[i] == '' or lines[i] == '\n':
			continue
		x = lines[i].split(',')
		cc = []
		for j in range(len(x)):
			cc = cc + [int(x[j])]
		S = S + [cc]

	R = [] # pair
	line = infile.readline()
	if line[len(line)-1] == '\n':
		line = line[0:len(line)-1]		
	lines = line.split(',;')
	for i in range(len(lines)):
		if lines[i] == '' or lines[i] == '\n':
			continue
		x = lines[i].split(',')
		cc = []
		for j in range(len(x)):
			cc = cc + [int(x[j])]
		R = R + [cc]

	infile.close()
	return C, A, P, S, R




def load_node(node_file):
	node_vec = []
	for lines in open(node_file):
		if lines[0] != '>':
			node_seq = lines[0:len(lines)-1]
			node_vec = node_vec + [node_seq]
	return node_vec



def mmilp_call(C, A, P, R, w1, w2, u, outfile):

	L = len(A) #path sum
	N = len(C) #node sum
	out = open(outfile,"w")
	print('from __future__ import print_function\n',file=out, end='')
	print('import sys\n',file=out, end='')
	print('import re\n',file=out, end='')
	print('from gurobipy import *\n',file=out,end='')
	print('def mmilp_solver(outfile):\n',file=out,end ='')

	print('\tm=100000000\n',file=out,end='')

	print('\tmodel = Model(\'mmip\')\n',file=out, end='')
	print('\tmodel.Params.TimeLimit = 300\n',file=out, end='')

	########################### var #####################################
	for i in range(L):
		print('\tx'+str(i)+' = model.addVar(vtype=GRB.BINARY)\n',file=out,end='')
	for i in range(L):
		print('\tt'+str(i)+' = model.addVar(lb=0,vtype=GRB.CONTINUOUS)\n',file=out,end='')

	for i in range(N):
		for j in range(i,N):
			if C[i][j] > 0:
				print('\tz'+str(i)+str(j)+' = model.addVar(vtype=GRB.BINARY)\n',file=out,end='')
				#print('\ty'+str(i)+str(j)+' = model.addVar(lb=0,vtype=GRB.CONTINUOUS)\n',file=out,end='')
	########################## obj ######################################
	print('\tmodel.setObjective('+str(w1)+'*(',file=out,end='')

	for i in range(L):
		if i < L-1:
			print('x'+str(i)+'+',file=out,end='')
		else:
			print('x'+str(i)+')+'+str(w2)+'*(',file=out,end='')

	for i in range(N):
		if i < N-1:
			print('z'+str(i)+str(i)+'+',file=out,end='')
		else:
			print('z'+str(i)+str(i)+'), GRB.MINIMIZE)\n',file=out,end='')
				
	########################## path tag constraint ############################
	for i in range(L):
		print('\tmodel.addConstr(x'+str(i)+' <= m*t'+str(i)+')\n',file=out,end='')
		print('\tmodel.addConstr(t'+str(i)+' <= m*x'+str(i)+')\n',file=out,end='')

	######################## node coverage constraint ########################
	for i in range(N):
		temp_str = ''
		tag = 0
		for k in range(L):
			if A[k][i] == 0:
				continue
			if tag == 0:
				temp_str = str(A[k][i])+'*t'+str(k)
				tag = tag + 1
			else:
				temp_str = temp_str + ' + ' + str(A[k][i])+'*t'+str(k) 
			
		if temp_str == '':
			continue
		print('\tmodel.addConstr(m*z'+str(i)+str(i)+' >= '+temp_str+' - '+str((1+u)*C[i][i])+')\n',file=out,end='')
		print('\tmodel.addConstr(m*z'+str(i)+str(i)+' >= '+str((1-u)*C[i][i])+' - ('+temp_str+'))\n',file=out,end='')
				
	######################### node constraint ################################
	for i in range(N):
		temp_str = ''
		tag = 0
		for k in range(L):
			temp = A[k][i]
			if temp > 0:
				if tag == 0:
					temp_str = temp_str + 'x'+str(k)
				else:
					temp_str = temp_str + ' + '+'x'+str(k)
				tag = tag + 1
		if temp_str != '':
			print('\tmodel.addConstr('+temp_str+' >= 1)\n',file=out,end='')

	####################### edge constraint ##################################
	for i in range(N-1):
		for j in range(i+1,N):
			if C[i][j] > 0:
				temp_str = ''
				tag = 0
				for k in range(L):
					temp = A[k][i]*A[k][j]
					if temp == 0:
						continue
					fid = P[k].index(i)
					if fid == 0 and P[k][1] == j:
						if tag == 0:
							temp_str = temp_str + 'x'+str(k)
						else:
							temp_str = temp_str + ' + '+'x'+str(k)
						tag = tag + 1
					if fid > 0:
						if P[k][fid-1] == j:
							if tag == 0:
								temp_str = temp_str + 'x'+str(k)
							else:
								temp_str = temp_str + ' + '+'x'+str(k)
							tag = tag + 1
						if fid+1 < len(P[k]) and P[k][fid+1] == j:
							if tag == 0:
								temp_str = temp_str + 'x'+str(k)
							else:
								temp_str = temp_str + ' + '+'x'+str(k)
							tag = tag + 1
				if temp_str == '':
					continue
				print('\tmodel.addConstr('+temp_str+' >= 1)\n',file=out,end='')

	########################### pair constraint ###################################

	for i in range(len(R)):
		xi = R[i][0]
		yi = R[i][1]
		temp_str = ''
		tag = 0
		for k in range(L):
			temp = A[k][xi]*A[k][yi]
			if temp == 0:
				continue
			if tag == 0:
				temp_str = temp_str + 'x'+str(k)
				tag = 1
			else:
				temp_str = temp_str + ' + '+'x'+str(k)
		if temp_str != '':
			print('\tmodel.addConstr('+temp_str+' >= 1)\n',file=out,end='')

	############################### print ###############################################
	print('\tmodel.optimize()\n',file=out,end='')

	print('\tofid = open(outfile, \'w\')\n',file=out, end='')
	for i in range(L):
		if i < L-1:
			print('\tprint(str(x'+str(i)+'.x)+\',\', file=ofid, end=\'\')\n',file=out,end='')
		else:
			print('\tprint(str(x'+str(i)+'.x)+\'\\n\', file=ofid, end=\'\')\n',file=out,end='')

	for i in range(L):
		if i < L-1:
			print('\tprint(str(t'+str(i)+'.x)+\',\', file=ofid, end=\'\')\n',file=out,end='')
		else:
			print('\tprint(str(t'+str(i)+'.x)+\'\\n\', file=ofid, end=\'\')\n',file=out,end='')


	for i in range(N):
		if i < N-1:
			print('\tprint(str(z'+str(i)+str(i)+'.x)+\',\', file=ofid, end=\'\')\n',file=out,end='')
		else:
			print('\tprint(str(z'+str(i)+str(i)+'.x)+\'\\n\', file=ofid, end=\'\')\n',file=out,end='')

	print('mmilp_solver(\''+out_dir+'mmilp.results\')\n',file=out,end='')
	out.close()
	command = 'python ' + outfile + '>out_dir'+'gurobi.results'
	os.system(command)


#####################################################################################################



in_file = ""
u = 0.1
M = 100000000
w1 = 1.0
w2 = 1.0
w= 0.1
for i in range(len(sys.argv)):
	if sys.argv[i] == "-in":
		in_file = sys.argv[i+1]
	if sys.argv[i] == "-w1":
		w1 = float( sys.argv[i+1])
	if sys.argv[i] == "-w2":
		w2 = float( sys.argv[i+1])	
	if sys.argv[i] == "-u":
		u = float(sys.argv[i+1])
	if sys.argv[i] == "-w":
		w = float(sys.argv[i+1])

if in_file=="":
	sys.exit(-1);

graph_sum = 0
out_dir = ""
node_name = ""
par_name = ""

infile = open(in_file)

line = infile.readline()
line = line[0:len(line)-1]
graph_sum = int(line)

line = infile.readline()
out_dir = line[0:len(line)-1]


line = infile.readline()
node_name = line[0:len(line)-1]

line = infile.readline()
par_name = line[0:len(line)-1]

line = infile.readline()
ave_cov = float(line[0:len(line)-1])

trans_name = out_dir + "MultiTrans.fa"
trans_file = open(trans_name, 'w')

x2 = 0
x5 = 0
x8 = 0
x9 = 0
single_sum = 0;
transcript_sum = 0

for i in range(graph_sum):
	#load node sequence
	node_vec = []
	node_file = out_dir + 'SG/' + str(i) + node_name
	for lines in open(node_file):
		if lines[0] != '>':
			node_seq = lines[0:len(lines)-1]
			node_vec = node_vec + [node_seq]	

	par_file = out_dir + 'SG/' + str(i) + par_name
	X = load_par(par_file)
	C = X[0]
	A = X[1]
	P = X[2]
	S = X[3]
	R = X[4]

	if len(node_vec) == 1:				
		#trim
		#if len(node_vec[0]) < 200:
		#	continue
		#if C[0][0] < 1:
		#	continue
		#print(i)
		#if len(node_vec[0]) < 500 and C[0][0] < ave_cov * 0.5 and  C[0][0] <= 4:
		#	x5 = x5 + 1
		#	continue
		#if len(node_vec[0]) < 800 and C[0][0] < ave_cov*0.1 and C[0][0] <= 2:
		#	x8 = x8 + 1
		#	continue
		if len(node_vec[0]) < 200 or len(node_vec[0])*C[0][0] < 0.1*ave_cov:
			continue

		trans_seq = node_vec[0] + '\n'
		trans_id = '>' + str(transcript_sum) + '_'+str(i) + '_0\n'
		print(trans_id, file=trans_file, end='')
		print(trans_seq, file=trans_file, end='')
		transcript_sum = transcript_sum + 1
		single_sum = single_sum + 1;
		continue

	w1 = 0
	#w1 = 1.0 / len(A)
	w2 = 1.0 

	outfile = out_dir + "mmilp.py"
	mmilp_call(C, A, P, R, w1, w2, u, outfile)
	results_name = out_dir + 'mmilp.results'
	out_file = open(results_name)
	line = out_file.readline()
	line = line[0:len(line)-1]
	results_vec = line.split(',')
	line = out_file.readline()
	line = line[0:len(line)-1]
	cov_vec = line.split(',')
	line = out_file.readline()
	line = line[0:len(line)-1]
	score_vec = line.split(',')
	out_file.close()
	

	satisfied_sum = 0;

	for j in range(len(score_vec)):
		if score_vec[j] == '1.0':
			satisfied_sum = satisfied_sum + 1

	if satisfied_sum > w*len(node_vec) and len(P) < 10:
		x9 = x9 + 1
		for j in range(len(P)):			
			trans_seq = ''
			for k in range(len(P[j])):
				node_id = P[j][k]
				node_str = node_vec[node_id]
				if S[j][k] == 0:
					node_str = revcmps(node_vec[node_id])
				if k > 0:
					trans_seq = trans_seq + node_str[int(C[P[j][k-1]][node_id]):len(node_vec[node_id])]
				else:
					trans_seq = node_str
			if len(trans_seq) > 0:
				trans_seq = trans_seq + '\n'

				trans_id = '>' + str(transcript_sum) + '_'+str(i)+'_'+str(j)+'\n'
				print(trans_id, file= trans_file, end='')
				print(trans_seq, file=trans_file, end='')
				transcript_sum = transcript_sum + 1
		continue;

	x2 = x2 + 1
	
	#output transcripts
	graph_trans_id = 0
	for j in range(len(results_vec)):
		if results_vec[j] == '1.0':
			trans_id = '>' + str(transcript_sum) + '_'+str(i) + '_'+str(graph_trans_id)+'\n'
			trans_seq = ''
			for k in range(len(P[j])):
				node_id = P[j][k]
				node_str = node_vec[node_id]
				if S[j][k] == 0:
					node_str = revcmps(node_vec[node_id])
				if k > 0:
					trans_seq = trans_seq + node_str[int(C[P[j][k-1]][node_id]):len(node_str)]
				else:
					trans_seq = node_str
			trans_seq = trans_seq + '\n'
			#if len(trans_seq) > 200:
			print(trans_id, file=trans_file, end='')
			print(trans_seq, file=trans_file, end='')
			transcript_sum = transcript_sum + 1
			graph_trans_id = graph_trans_id + 1

	#output all transcripts

trans_file.close()
command = "cat " + out_dir + "part.fa" + " >> " + out_dir + "MultiTrans.fa"
os.system(command)  
################ out put all path ###########################


















