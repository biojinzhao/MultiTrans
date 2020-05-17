
#ifndef AssemblyGRAPH_H
#define AssemblyGRAPH_H


#include<algorithm>
#include<assert.h>
#include<fstream>
#include<map>
#include<numeric>
#include<iomanip>
#include<set>
#include<string>
#include <stdlib.h>
#include<sstream>
#include<vector>
#include<list>


using namespace std;

extern string g_out_dir;
extern string g_graph_file;
extern string g_mapped_file;
extern string g_edge_file;
extern int g_max_path_sum;
extern bool g_help;


string revcomp (const string& kmer);
bool is_not_num(string s);


class Node {

public:

	string sequence;
	vector<pair<int, int> > parents;
	vector<pair<int, int> > children;
	vector<pair<int, int> > pair_nodes;
	double coverage;
	int reads_sum;
	string s_name;
	Node(): sequence("") {};

	Node(const string& mysequence){
		sequence = mysequence;
	}

	Node(const Node& node){

		sequence = node.sequence;
		parents = node.parents;
		children = node.children;
		coverage = node.coverage;
		s_name = node.s_name;
		reads_sum = node.reads_sum;
		pair_nodes = node.pair_nodes;
	}

	void remove() {
		sequence = "";
		coverage = 0.0;
		parents.clear();
		children.clear();
		s_name = "";
		reads_sum = 0;
		pair_nodes.clear();
	}

	bool add_pair(pair<int, int> p) {

		if(p.first < 0)
			return false;

		if(!pair_nodes.empty()) {

			for (int i = 0; i < pair_nodes.size(); i++) 
				if (pair_nodes[i].first == p.first) {
					pair_nodes[i].second = pair_nodes[i].second + p.second;
					return true;
				}
		}

		this -> pair_nodes.push_back(p);
		return true;

	}

	bool add_child(pair<int, int> child) {

		if(child.first < 0)
			return false;

		if(!children.empty()) {

			for (int i = 0; i < children.size(); i++) 
				if (children[i].first == child.first)
					return false;
		}

		this -> children.push_back(child);
		return true;

	}

	bool add_parent(pair<int, int> parent) {

		if(parent.first < 0)
			return false;
		if (!parents.empty()) {
			for (int i = 0; i < parents.size(); i++)
				if (parents[i].first == parent.first)
					return false;
		}

		this -> parents.push_back(parent);
		return true;

	}


	bool delete_child(int child) {

		if(child < 0)
			return false;

		vector<pair<int, int> >::iterator it = children.begin();

		for( ; it != children.end(); it++) {

			if (it->first == child)
				break;

		}

		if(it != children.end()) {
			children.erase(it);
			return true;

		} else {

			return false;

		}

	}

	bool delete_parent(int parent) {

		if(parent < 0)
			return false;

		vector<pair<int, int> >::iterator it = parents.begin();

		for ( ; it != parents.end(); it++) {

			if (it->first == parent)
				break;

		}

		if (it != parents.end()) {
			parents.erase(it);
			return true;

		} else {

			return false;

		}

	}

};



class Graph {

public:

	vector<Node> node_set;
	size_t node_sum;
	int mid_node_id;
	vector<int> node_order;
	map<int, bool> circle_node;
	Graph();
	void set_parents();
	int add_node(Node& node);
	void add_reverse_nodes();
	void load_graph(const string& graph_path, const string& mapped_file);
	void load_graph(const string& graph_path, const string& edge_file, const string& mapped_file);
	void get_non_trivial_graph_nodes(const string& graph_path, string out_dir);
	void save_graph(string out_dir);
	void add_to_comp(vector<int>& is_out, vector<int>& comp_vec, int current_id, int connect_id);
	bool has_circle(vector<int>& comp_vec);
	bool is_circle(int p, int q, set<int>& checked, vector<int>& c_path);
	bool is_circle(int p, int q, set<int>& checked);
	int divide_graph(vector<int>& comp_vec, vector<int>& is_out);
	int delete_circle(vector<int>& comp_vec);
	void topological_sort(map<int, int>& comp_map, vector<int>& node_order);
	void dfs_visit(int i, map<int, int>& comp_map, vector<int>& node_color, vector<int>& node_order);
	void node_path_sum(map<int, int>& comp_map, vector<int>& sum_vec, int p);
	int total_path_num(map<int, int>& comp_map, vector<int>& node_order);
	void get_all_path(int node, vector<int>& path_node, vector<vector<int> >& path_vec, vector<int>& path_tag, int& bad_sum, vector<int>& is_out);
	void get_all_path_sum(int node, vector<int>& path_node, int& path_sum, vector<int>& path_tag, int& bad_sum, vector<int>& is_out);
	void get_left_heavy_path(std::vector<int>& node_order, int pos, std::vector<int>& path_node);
	void get_right_heavy_path(std::vector<int>& node_order, int pos, std::vector<int>& path_node);
	void add_path_edge(int p, int q, int max_overlap);
	int compute_edge_score(int p, int q, map<pair<int, int>, int>& edge_score, vector<bool>& is_checked);
	void add_parent_node(int p, set<int>& parent_set, vector<bool>& is_checked, int& check_len);
	void add_child_node(int p, set<int>& parent_set, vector<bool>& is_checked, int& check_len);
	void greedy_path(pair<int, int> edge, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec);
	void greedy_left_path(int p, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec);
	void greedy_right_path(int p, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec);

};



#endif
