# MultiTrans 


### About MultiTrans

MultiTrans is a tool for transcripts assmebly, which formulated the transcripts assembly problem into a mixed integer linear programming. 


### Prerequistites 

C++, Python, and  Gurobi are required to be installed.

To use the paired-end information of assembly graphs, BWA is also required to be installed.


### Quick start

**Assembly graph

python MultiTran.py -m 1 -g /path/to/assembly/graph/assembly_graph_with_scaffolds.gfa

**Splicing graph

python MultiTran.py -m 2 -g /path/to/splicing/graph/RawGraphs_/


### Usage

 --model (or -m) <int>: 1 assembly graph, 2 splcing graph
 
 --graph_file (or -g) <string>: *.gfa file of assembly graph  or  directory of splicing graph
                                
 --out_dir (or -o): directory of output
 
 --reads1 (or -r1) <string>: paired-end reads
 
 --reads2 (or -r2) <string>: paired-end reads
 
 --test (or -t): testing the software
   
 --help (or -h): MultiTrans usage


###### Extracting common transcripts from multiple transcriptome.fa ############

### Example

python multiple_trans_analyze.py -t trans1.fa,trans2.fa,trans3.fa -o /out_dir/

### Usage

 --transcripts_file_list (or -t): list of transcriptome.fa separated by commas
 --out_dir (or -o): directory of output


### Feedback and bug reports

Your comments, bug reports, and suggestions are very welcomed. They will help us to further improve MultiTrans. If you have any troubles running MultiTrans, please contact us (Email: zhaojin@mail.sdu.edu.cn).

