#include "AssemblyGraph.h"
#include <iostream>
#include <ctime>
#include <errno.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include <time.h>
#include <getopt.h>


using namespace std;

static const char *short_options = "g:e:m:o:h";

struct option opts[] = {
	{"graph_file",  required_argument,   0,  'g'},
	{"edge_file",  required_argument,   0,  'e'},
	{"mapped_file",  required_argument,   0,  'm'},
	{"out_dir",  required_argument,   0,   'o'},
	{"help", no_argument,  0,   'h'},
	{0,0,0,0} // terminator
};


string usage() {

	stringstream usage_info;
	usage_info
		<< endl
		<< "===============================================================================" << endl
		<< " MultiTrans Usage " << endl
		<< "===============================================================================" << endl
		<< " ** Required: **" << endl
		<< "  -g <string>: graph file" << endl
		<< "  -e <string>: edge file" << endl
		<< "  -m <string>: mapped file" << endl
		<< "  -o <string>: output directory " << endl
		<< " ** Options: **" <<endl
		<< "  -h : help information. " << endl
		<< "===============================================================================" << endl
		<< endl;

	return usage_info.str();

}


int parse_options(int argc, char* argv[]) {

	int option_index = 0;
	int next_option;
	do {
		next_option = getopt_long(argc, argv, "g:e:m:o:h", opts, &option_index);
		switch (next_option) {
		case -1:
			break;
		case 'g':
			g_graph_file = optarg;
			break;
		case 'e':
			g_edge_file = optarg;
			break;
		case 'm':
			g_mapped_file = optarg;
			break;
		case 'o':
			g_out_dir = optarg;
			break;
		case 'h':
			g_help = true;
			break;
		default:
			cout << usage();
			exit(1);
		}

	} while (next_option != -1);

	if (g_help) {
		cout << usage();
		exit (1);
	}

	return 0;

}



int main(int argc, char* argv[]){

	int parse_ret = parse_options(argc,argv);
	if (parse_ret)
		return parse_ret;

	//string out_dir = "/data/RNA-Seq/human/samllSimulation/1000.mix/";
	//string file_name ="/data/RNA-Seq/human/samllSimulation/1000.mix/assembly_graph_with_scaffolds.gfa";	
	//string file_name ="/data/RNA-Seq/SRR5133163/rnaSpades/assembly_graph_with_scaffolds.gfa";	
	
	Graph graph;
	
	graph.load_graph(g_graph_file, g_mapped_file);
	graph.save_graph(g_out_dir);
	
	return 1;

}

