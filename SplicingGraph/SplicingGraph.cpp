
#include "SplicingGraph.h"
#include<iostream>

using namespace std;

string g_graph_file = "";
string g_out_dir = "";
string g_mapped_file = "";
int g_max_path_sum = 500;
bool g_help = false;
int g_graph_sum = 0;

bool is_not_num(string s) {

	stringstream sin(s);
	double t;
	char p;	
	if (!(sin >> t))
		return true;
	if (sin >> p)
		return true;
	else
		return false;
}	


string revcomp (const string& kmer) {

	string revstring;

	for(int i = kmer.size() -1; i >= 0;i--) {

		char c = kmer[i];
		char revchar;

		switch (c) {

		case 'g':
			revchar = 'c';
			break;

		case 'G':
			revchar = 'C';
			break;

		case 'a':
			revchar = 't';
			break;

		case 'A':
			revchar = 'T';
			break;

		case 't':
			revchar = 'a';
			break;

		case 'T':
			revchar = 'A';
			break;

		case 'c':
			revchar = 'g';
			break;

		case 'C':
			revchar = 'G';
			break;

		default:
			revchar = 'N';
		}

		revstring += revchar;

	}

	return (revstring);

}



Graph::Graph() {

	node_sum = 0;

}






void Graph::set_parents() {

	for (int i = 0; i < node_sum; ++i) {

		if (!node_set[i].parents.empty())
			node_set[i].parents.clear();

	}

	for (size_t i = 0; i < node_sum; ++i) {

		vector<pair<int, double> >::iterator it;

		for (it = node_set[i].children.begin(); it != node_set[i].children.end(); ++it)
			node_set[it->first].add_parent(make_pair(i,it->second));

	}

}



int Graph::add_node(Node& node) {

	node_set.push_back(node);
	return (node_sum++);

}


void Graph::get_all_path_sum(int node, vector<int>& path_node, int& path_sum, vector<int>& path_tag, int& bad_sum, vector<int>& is_out) { 
//cout << "get all path" << endl;
	path_node.push_back(node);
	is_out[node] = -2;
	path_tag[node] = path_tag[node] + 1;
	if (path_tag[node] > 1 && bad_sum < 1)
		bad_sum = bad_sum + 1;
	if (node_set[node].children.size() == 0 || (bad_sum == 1 && path_tag[node] > 1)) {
		path_sum = path_sum +1;
		if (path_sum > g_max_path_sum)
			return;
		path_tag[node] = path_tag[node] - 1;
		path_node.pop_back();

	} else {
		for (int i = 0; i < node_set[node].children.size(); ++i) {
			int c = node_set[node].children[i].first;
			get_all_path_sum(c, path_node, path_sum, path_tag, bad_sum, is_out);
			if (path_sum > g_max_path_sum)
				return;

		}
		if (path_tag[path_node[path_node.size()-1]] > 1) 
			bad_sum = bad_sum - 1;
		path_tag[path_node[path_node.size()-1]] = path_tag[path_node[path_node.size()-1]] - 1;
		path_node.pop_back();
	}
	
//cout << "done" << endl;
}

void Graph::get_all_path(int node, vector<int>& path_node, vector<vector<int> >& path_vec, vector<int>& path_tag, int& bad_sum, vector<int>& is_out) { 
//cout << "get all path" << endl;
	path_node.push_back(node);
	is_out[node] = -2;
	path_tag[node] = path_tag[node] + 1;
	if (path_tag[node] > 1 && bad_sum < 1)
		bad_sum = bad_sum + 1;
	if (node_set[node].children.size() == 0 || (bad_sum == 1 && path_tag[node] > 1)) {
		path_vec.push_back(path_node);
		path_tag[node] = path_tag[node] - 1;
		path_node.pop_back();
	} else {
		for (int i = 0; i < node_set[node].children.size(); ++i) {
			int c = node_set[node].children[i].first;
			get_all_path(c, path_node, path_vec, path_tag, bad_sum, is_out);		
		}
		if (path_tag[path_node[path_node.size()-1]] > 1) 
			bad_sum = bad_sum - 1;
		path_tag[path_node[path_node.size()-1]] = path_tag[path_node[path_node.size()-1]] - 1;	
		path_node.pop_back();
	}
	
//cout << "done" << endl;
}




void Graph::load_all_graph(const string& graph_path) {

	cout << "load graph..." << endl;
	for (int i = 0; i < g_graph_sum; ++i) {
		string path = g_graph_file + "comp"+ to_string(i) + ".rg.debug";
		load_graph(path, i);
	}

}


void Graph::save_graph(string out_dir) {

	cout << "node_sum: " << node_set.size() << endl;
	int xyz = 0;
	for (int i = 0; i < node_set.size(); ++i) 
		xyz = xyz + node_set[i].children.size();
	cout << "edge_sum: " << xyz << endl;

	int total_pair = 0;
	for (int i = 0; i < node_set.size(); ++i) {
		total_pair = total_pair + node_set[i].pair_nodes.size();
	}

	cout << "preprocessing..." << endl;
	double total_ave = 0.0;
	int total_len = 0;
	for (int i = 0; i < node_set.size(); ++i) {
		total_ave = total_ave + node_set[i].coverage*node_set[i].sequence.length();
		total_len = total_len + node_set[i].sequence.length();
	}
	total_ave = total_ave / total_len;
	total_ave = total_ave * total_len / node_set.size();
	map<pair<int, int>, int> edge_score;
	const int current_node_sum = node_set.size();
	vector<int> is_out(current_node_sum, -1);
	vector<int> path_tag(current_node_sum, 0);
	vector<bool> is_checked(current_node_sum, false);
	vector<int> comp_vec;
	int current_id = -1;
	int connect_id = 0;
	int mate_id, p;
	vector<int> path_node;
	int delete_sum = 0;

	string long_file_name = out_dir + "part.fa";
	fstream long_file;
	long_file.open(long_file_name.c_str(), fstream::out);


	while (current_id + 1 < node_set.size()) {

		current_id = current_id + 1;
		if (is_out[current_id] >= 0)
			continue;
		comp_vec.clear();
		bool is_conflict = false;
		add_to_comp(is_out, comp_vec, current_id, connect_id);

		for (int i = 0; i < path_tag.size(); ++i)
			path_tag[i] = 0;
		int path_sum = 0;
		for (int i = 0; i < comp_vec.size(); ++i) {
			p = comp_vec[i];
			if (node_set[p].parents.size() == 0) {
				path_node.clear();
				int bad_sum = 0;
				get_all_path_sum(p, path_node, path_sum, path_tag, bad_sum, is_out);
				if (path_sum > g_max_path_sum)
					break;
			}
		}
		for (int i = 0; i < comp_vec.size(); ++i) {
			if (is_out[comp_vec[i]] >= 0) {
				path_node.clear();
				int bad_sum = 0;
				get_all_path_sum(comp_vec[i], path_node, path_sum, path_tag, bad_sum, is_out);
				if (path_sum > g_max_path_sum)
					break;
			}
		}

		if (path_sum <= g_max_path_sum) {	
			for (int i = 0; i < comp_vec.size(); ++i) 
				is_out[comp_vec[i]] = connect_id;
			vector<vector<int> > path_vec;
			for (int i = 0; i < path_tag.size(); ++i)
				path_tag[i] = 0;
			for (int i = 0; i < comp_vec.size(); ++i) {
				p = comp_vec[i];
				if (node_set[p].parents.size() == 0) {
					path_node.clear();
					int bad_sum = 0;
					get_all_path(p, path_node, path_vec, path_tag, bad_sum, is_out);
				}
			}
			for (int i = 0; i < comp_vec.size(); ++i) {
				if (is_out[comp_vec[i]] >= 0) {
					path_node.clear();
					int bad_sum = 0;
					get_all_path(comp_vec[i], path_node, path_vec, path_tag, bad_sum, is_out);
				}
			}
			//output

			string par_file_name = out_dir + "SG/" + to_string(connect_id)+".par";
			string node_file_name = out_dir + "SG/" + to_string(connect_id) + ".node";
			fstream par_file, node_file;
			node_file.open(node_file_name.c_str(), fstream::out);

			for (int i = 0; i < comp_vec.size(); ++i) {
				node_file << ">" << node_set[comp_vec[i]].coverage << "\t" + node_set[comp_vec[i]].s_name << endl;
				node_file << node_set[comp_vec[i]].sequence << endl;
			}
			node_file.close();		
			const int comp_size = comp_vec.size();
			vector<vector<double> > adj_matrix;
			for (int i = 0; i < comp_size; ++i) {
				vector<double> temp(comp_size, 0.0);
				adj_matrix.push_back(temp);
			}
			map<int, int> comp_map;			
			for (int i = 0; i < comp_vec.size(); ++i) {
				adj_matrix[i][i] = node_set[comp_vec[i]].coverage;
				comp_map[comp_vec[i]] = i;
			}			
			for (int i = 0; i < comp_vec.size(); ++i) {
				for (int j = 0; j < node_set[comp_vec[i]].children.size(); ++j) {
					int c = node_set[comp_vec[i]].children[j].first;					
					adj_matrix[i][comp_map[c]] = node_set[comp_vec[i]].children[j].second;
					adj_matrix[comp_map[c]][i] = node_set[comp_vec[i]].children[j].second;
				}
			}
			vector<vector<int> > out_path;
			for (int i = 0; i < path_vec.size(); ++i) {
				vector<int> temp(comp_size, 0);
				out_path.push_back(temp);
			}

			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					int pv = path_vec[i][j];
					out_path[i][comp_map[pv]] = out_path[i][comp_map[pv]] + 1;
				}
			}
			par_file.open(par_file_name.c_str(), fstream::out);
			
			for (int i = 0; i < comp_size; ++i){
				for (int j = 0; j < comp_size; ++j) {
					par_file << adj_matrix[i][j] << ",";
				}
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < out_path.size(); ++i) {
				for (int j = 0; j < comp_size; ++j) {
					par_file << out_path[i][j] << ",";
				}					
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					int pv = path_vec[i][j];
					par_file << comp_map[pv] << ",";
				}
				par_file << ";";
			}
			par_file << endl;

			for (int i = 0; i < comp_vec.size(); ++i) {
				int xi = comp_vec[i];
				for (int j = 0; j < node_set[comp_vec[i]].pair_nodes.size(); ++j) {
					int yi = node_set[comp_vec[i]].pair_nodes[j].first;
					if (comp_map.find(yi) != comp_map.end())
						par_file << i << "," << comp_map[yi] << ",;";
				}
			}
			par_file << endl;

			for (int i = 0; i < comp_vec.size(); ++i) 
				is_out[comp_vec[i]] = connect_id;

			connect_id = connect_id + 1;
		} else { // path_sum > g_max_path_sum

			delete_sum = delete_sum + 1;
			vector<pair<int, int> > min_candidate;
			pair<int, int> edge;
			int min_score = 0.0;
			for (int i = 0; i < comp_vec.size(); ++i) {
				for (int j = 0; j < node_set[comp_vec[i]].children.size(); ++j) {
					edge.first = comp_vec[i];
					edge.second = node_set[comp_vec[i]].children[j].first;
					int score = compute_edge_score(edge.first, edge.second, edge_score, is_checked);
					if (min_candidate.size() == 0) {
						min_score = score;
						min_candidate.push_back(edge);
					} else {
						if (min_score > score) {
							min_score = score;
							min_candidate.clear();
							min_candidate.push_back(edge);
						}
						if (min_score == score) {
							min_candidate.push_back(edge);
						}
					}
				}
			}
	
			edge = min_candidate[0];
			double min_score2 = node_set[edge.first].coverage + node_set[edge.second].coverage;
			for (int i = 1; i < min_candidate.size(); ++i) {
				int current_score = node_set[min_candidate[i].first].coverage + node_set[min_candidate[i].second].coverage;
				if (current_score < min_score2) {
					min_score2 = current_score;
					edge = min_candidate[i];
				}
			}

			node_set[edge.first].delete_child(edge.second);
			node_set[edge.second].delete_parent(edge.first);

			/*if (min_score2 >= 0.2*total_ave || min_score2 > 4) {
				path_node.clear();
				greedy_path(edge, path_node, is_checked, comp_vec);
				string transcript = node_set[path_node[0]].sequence;
				string trans_name = node_set[path_node[0]].s_name;
				for (int i = 1; i < path_node.size(); ++i) {
					for (int j = 0; j < node_set[path_node[i-1]].children.size(); ++j) {
						if (node_set[path_node[i-1]].children[j].first == path_node[i]) {
							trans_name = trans_name + "," + node_set[path_node[i]].s_name;
							transcript = transcript + node_set[path_node[i]].sequence;
							break;
						}
					}
				}
				if (transcript.length() >= 200)
					long_file << ">" << trans_name << endl << transcript << endl;
			}*/

			path_node.clear();
			greedy_path(edge, path_node, is_checked, comp_vec);
			string transcript = node_set[path_node[0]].sequence;
			string trans_name = node_set[path_node[0]].s_name;
			for (int i = 1; i < path_node.size(); ++i) {
				for (int j = 0; j < node_set[path_node[i-1]].children.size(); ++j) {
					if (node_set[path_node[i-1]].children[j].first == path_node[i]) {
						trans_name = trans_name + "," + node_set[path_node[i]].s_name;
						transcript = transcript + node_set[path_node[i]].sequence;
						break;
					}
				}
			}
			//if (transcript.length() >= 200 && transcript.length()*min_score2 >= 0.2*total_ave)
			if (transcript.length() > 0)
				long_file << ">" << trans_name << endl << transcript << endl;
		

		
			for (int i = 0; i < comp_vec.size(); ++i) 
				is_out[comp_vec[i]] = -1;

			current_id = current_id - 1;
		}
	}

	string info_file_name = out_dir + "graph.info";
	fstream info_file;
	info_file.open(info_file_name.c_str(), fstream::out);
	info_file << connect_id << endl;
	info_file << out_dir << endl;
	info_file << ".node" << endl;
	info_file << ".par" << endl;
	info_file << total_ave << endl;
	info_file.close();
	//cout << "deleted_edge: " << delete_sum <<endl;
	cout << "graph_sum:" << connect_id << endl;


}




void Graph::greedy_path(pair<int, int> edge, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec) {
	vector<int> temp_vec;
	greedy_left_path(edge.first, temp_vec, is_checked, comp_vec);
	for (int i = temp_vec.size()-1; i >= 0; --i) {
		path_node.push_back(temp_vec[i]);
	}
	greedy_right_path(edge.second, path_node, is_checked, comp_vec);
	for (int i = 0; i < path_node.size(); ++i)
		is_checked[path_node[i]] = false;
}


void Graph::greedy_left_path(int p, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec) {
	if (is_checked[p])
		return;
	path_node.push_back(p);
	is_checked[p] = true;
	
	if (node_set[p].parents.size() > 0) {
		int max_p = node_set[p].parents[0].first;
		double max_cov = node_set[max_p].coverage;
		for (int i = 0; i < node_set[p].parents.size(); ++i) {
			int x = node_set[p].parents[i].first;
			if (node_set[x].coverage > max_cov) {
				max_cov = node_set[x].coverage;
				max_p = x;
			}
		}
		greedy_left_path(max_p, path_node, is_checked, comp_vec);
	}
}


void Graph::greedy_right_path(int p, vector<int>& path_node, vector<bool>& is_checked, vector<int>& comp_vec) {
	if (is_checked[p])
		return;
	path_node.push_back(p);
	is_checked[p] = true;
	
	if (node_set[p].children.size() > 0) {
		int max_p = node_set[p].children[0].first;
		double max_cov = node_set[max_p].coverage;
		for (int i = 0; i < node_set[p].children.size(); ++i) {
			int x = node_set[p].children[i].first;
			if (node_set[x].coverage > max_cov) {
				max_cov = node_set[x].coverage;
				max_p = x;
			}
		}
		greedy_right_path(max_p, path_node, is_checked, comp_vec);
	}
}


int Graph::compute_edge_score(int p, int q, map<pair<int, int>, int>& edge_score, vector<bool>& is_checked) {
//cout << "compute score" << p << " " << q << endl;
	pair<int, int> edge;

	edge.first = p;
	edge.second = q;
	if (edge_score.find(edge) != edge_score.end()) {
		return edge_score[edge];
	}

	set<int> parent_set;
	set<int> child_set;
	parent_set.insert(p);
	child_set.insert(q);

	int check_len = 0;
	add_parent_node(p, parent_set, is_checked, check_len);
	check_len = 0;
	add_child_node(q, child_set, is_checked, check_len);

	set<int>::iterator it;
	set<int> mapped_node;
	map<int, int> mapped_node_score;
	for (it = parent_set.begin(); it != parent_set.end(); ++it) {
		for (int i = 0; i < node_set[*it].pair_nodes.size(); ++i) {
			int n = node_set[*it].pair_nodes[i].first;
			int pre_size = mapped_node.size();
			mapped_node.insert(n);
			if (mapped_node.size() > pre_size) {
				mapped_node_score[n] = node_set[*it].pair_nodes[i].second;
			} else {
				mapped_node_score[n] = mapped_node_score[n] + node_set[*it].pair_nodes[i].second;
			}
		}
	}

	vector<int> matched;
    set_intersection(child_set.begin(),child_set.end(),mapped_node.begin(),mapped_node.end(),inserter(matched,matched.end()));

	int score = 0;
	for (int i = 0; i < matched.size(); ++i) {
		score = score + mapped_node_score[matched[i]];
	}
	edge_score[edge] = score;

	for (it = parent_set.begin(); it != parent_set.end(); ++it) 
		is_checked[*it] = false;
	for (it = child_set.begin(); it != child_set.end(); ++it) 
		is_checked[*it] = false;
	return score;
}


void Graph::add_parent_node(int p, set<int>& parent_set, vector<bool>& is_checked, int& check_len) {
	if (node_set[p].parents.size() == 0 || is_checked[p]) {
		return;
	} else {
		is_checked[p] = true;
		for (int i = 0; i < node_set[p].parents.size(); ++i) {
			int parent_id = node_set[p].parents[i].first;
			if (!is_checked[parent_id] && check_len < 500) {
				parent_set.insert(parent_id);
				check_len = check_len + node_set[parent_id].sequence.length();
				int current_len = check_len;
				add_parent_node(parent_id, parent_set, is_checked, current_len);
			}
		}
		
	}
}

void Graph::add_child_node(int p, set<int>& child_set, vector<bool>& is_checked, int& check_len) {
	if (node_set[p].children.size() == 0 || is_checked[p])
		return;
	else {
		is_checked[p] = true;
		for (int i = 0; i < node_set[p].children.size(); ++i) {
			int child_id = node_set[p].children[i].first;
			if (!is_checked[child_id] && check_len < 500) {
				child_set.insert(child_id);
				check_len = check_len +node_set[child_id].sequence.length();
				int current_len = check_len;
				add_child_node(child_id, child_set, is_checked, current_len);
			}
		}
	}
}




void Graph::load_graph(const string& graph_path, int graph_id) {
	ifstream ifs(graph_path);
	vector<string> vec_edges;
	vector<string> vec_pairs;
	vector<double> Node_cov;
	vector<string> Node_seq;

	map<int, int> node_map;
	int node_id = 0;

	double cov;
	string temp;
	vector<pair<int,int> > edge_vec;
	int i,j,k,t,sum_j;
	map<string, int>::iterator it;
	getline(ifs,temp);
	while (1) {
		getline(ifs,temp);
		if (temp=="** Inhibit_edges **")
			break;
		else
			vec_edges.push_back(temp);
	}
	while (1) {
		getline(ifs,temp);
		if (temp=="** Pair_edges **")
			break;
	}
	while (1) {
		getline(ifs,temp);
		if (temp=="** Nodes **")
			break;
		else
			vec_pairs.push_back(temp);
	}
	while (getline(ifs,temp)) {
        if (temp.substr(0,10)=="node id = ") {
			Node node;
            for (i=10;i<temp.size();i++) {
                if (temp.substr(i-6,6)=="cov = "){
                    for (j=0;j<temp.size();j++){
                        if (temp.substr(i+j,2)=="se"){
							node.coverage = atof(temp.substr(i,j-1).c_str());
                            break;
                        }
                    }
                }
            }
            getline(ifs,temp);
			if (temp[temp.length()-1] == '\n')
				node.sequence = temp.substr(0, temp.length()-1);
			else
				node.sequence = temp;
			node.s_name = to_string(graph_id) +"_" + to_string(node_id);
			node_map[node_id] = node_set.size();
			node_set.push_back(node);
			node_id = node_id + 1;
        }
    }

	ifs.close();
	if (vec_edges.size()>=1) {
		for (i=0;i<vec_edges.size();i++) {
			int edge_left = -1;
			int edge_right = -1;
			for (j=0;j<vec_edges[i].size();j++) {
				if (j<vec_edges[i].size()-1 && vec_edges[i][j+1]=='-') {
					if (j==0) {
						edge_left = atoi(vec_edges[i].substr(j,1).c_str());
					} else {
						sum_j=1;
						while (1) {
							if ((j-sum_j==0) || (j-sum_j!=0 && vec_edges[i][j-sum_j]==' '))
								break;
							else
								sum_j++;
						}
						if (j-sum_j==0) {
							edge_left = atoi(vec_edges[i].substr(j-sum_j,sum_j+1).c_str());
						} else {
							edge_left = atoi(vec_edges[i].substr(j-sum_j+1,sum_j).c_str());
						}
					}
				}
				if (j>1 && vec_edges[i][j-1]=='>') {
					sum_j=1;
					while (1) {
						if (vec_edges[i][j+sum_j]==':')
							break;
						else
							sum_j++;
					}
					edge_right = atoi(vec_edges[i].substr(j,sum_j).c_str());
					
				}
				if (j>1 && vec_edges[i][j-1]==':') {
					sum_j=1;
					while (1) {
						if (vec_edges[i][j+sum_j]==',' || vec_edges[i][j+sum_j]==';')
							break;
						else
							sum_j++;
					}
					//cout << edge_left << endl;
					//cout << edge_right <<endl;
					cov = atof(vec_edges[i].substr(j,sum_j).c_str());
					//cout << node_set.size() << endl;
					node_set[node_map[edge_left]].add_child(make_pair(node_map[edge_right], cov));
					node_set[node_map[edge_right]].add_parent(make_pair(node_map[edge_left], cov));	
					edge_vec.push_back(make_pair(edge_left, edge_right));				
				}
			}
		}

		if (vec_pairs.size()>=1) {
			for (i=0;i<vec_pairs.size();i++) {
				for (j=0;j<vec_pairs[i].size();j++) {
					if (vec_pairs[i][j]=='-') {
						int edge_id = atoi(vec_pairs[i].substr(0,j).c_str());
						for (k=j+2;k<vec_pairs[i].size();k++) {
							if (vec_pairs[i][k]=='-') {
								int node_id = atoi(vec_pairs[i].substr(j+2,k-j-2).c_str());
								vector<int> pairs;
								if (edge_vec[edge_id].second != node_id && edge_vec[edge_id].first != node_id) {
									int x1 = node_map[node_id];
									int x2 = node_map[edge_vec[edge_id].first];
									int x3 = node_map[edge_vec[edge_id].second];
									node_set[x1].add_pair(make_pair(x2, 1));
									node_set[x1].add_pair(make_pair(x3, 1));
									node_set[x2].add_pair(make_pair(x1, 1));
									node_set[x2].add_pair(make_pair(x3, 1));
									node_set[x3].add_pair(make_pair(x1, 1));
									node_set[x3].add_pair(make_pair(x2, 1));
								}
						    	break;
							}
						}
						break;
					}
				}
			}
		}//end of if (vec_pairs.size()>=1)

	}

}

bool Graph::has_path(int s, int t) {

	bool is_same = false;
	for (int i = 0; i < node_set[s].children.size(); ++i) {
		if (node_set[s].children[i].first == t) {
			return true;
		} 
	}

	for (int i = 0; i < node_set[s].children.size(); ++i) {
		is_same = has_path(node_set[s].children[i].first, t);
	}
	return is_same;
}

void Graph::save_graph2(string out_dir, vector<string>& trans_vec, int& connect_id, double& total_ave) {

cout << "save graph..." << endl;
	const int current_node_sum = node_set.size();
	vector<int> is_out(current_node_sum, -1);
	vector<int> final_vec;
	int current_id = -1;
	int mate_id;
	int mis_num = 0;
	int c_sum = 0;
	int is_same = 0;
	int p, q;


	double current_ave = 0;
	double current_min = 0;

	int xxx = 0;

	int big_sum = 0;

	while (current_id + 1 < node_set.size()) {
		current_id = current_id + 1;
		if (is_out[current_id] != -1)
			continue;
		bool is_continue = false;
		bool is_circle = false;
		final_vec.clear();
		add_to_comp(is_out, final_vec, current_id, connect_id);
		map<int, int> comp_map;
		for (int i = 0; i < final_vec.size(); ++i)
			comp_map[final_vec[i]] = i;
	   	vector<int> node_order;
		topological_sort(comp_map, node_order);
		int total_p_num = total_path_num(comp_map, node_order);	
		if (total_p_num > 100) {
			big_sum = big_sum + 1;
			int min_node = final_vec[0];
			double min_cov = node_set[final_vec[0]].coverage;
			for (int i = 0; i < final_vec.size(); ++i) {
				if (node_set[final_vec[i]].coverage < min_cov) {
					min_cov = node_set[final_vec[i]].coverage;
					min_node = final_vec[i];
				}
			}

			int is_father = 0;
			if (node_set[min_node].parents.size() > 0) {
				p = node_set[min_node].parents[0].first;
				min_cov = node_set[p].coverage;
				is_father = 1;
			} else {
				if (node_set[min_node].children.size() > 0) {
					p = node_set[min_node].children[0].first;
					min_cov = node_set[p].coverage;
					is_father = 2;
				}
			}
			for (int i = 0; i < node_set[min_node].parents.size(); ++i) {
				q = node_set[min_node].parents[i].first;
				if (node_set[q].coverage < min_cov) {
					p = q;
					min_cov = node_set[q].coverage;
					is_father = 1;
				}
			}
			for (int i = 0; i < node_set[min_node].children.size(); ++i) {
				q = node_set[min_node].children[i].first;
				if (node_set[q].coverage < min_cov) {
					p = q;
					min_cov = node_set[q].coverage;
					is_father = 2;
				}
			}

			if (is_father == 1) {
				vector<int> path_node;
				for (int i = node_order.size()-1; i >=0; --i) {
					if (node_order[i] == p) {
						get_left_heavy_path(node_order, i, path_node);
						path_node.push_back(p);
					} else {
						if (node_order[i] == min_node){
							get_right_heavy_path(node_order, i, path_node);								
							break;
						}
					}
				}
				string transcript = node_set[path_node[0]].sequence;
				for (int i = 1; i < path_node.size(); ++i) {
					for (int j = 0; j < node_set[path_node[i-1]].children.size(); ++j) {
						if (node_set[path_node[i-1]].children[j].first == path_node[i]) {
							transcript = transcript + node_set[path_node[i]].sequence;
							break;
						}
					}
				}
				if (transcript.length() >= 200 && min_cov > 0.001*total_ave / connect_id)
					trans_vec.push_back(transcript);
				node_set[min_node].delete_parent(p);
				node_set[p].delete_child(min_node);
			} 
			if (is_father == 2){
				vector<int> path_node;
				for (int i = node_order.size()-1; i >=0; --i) {
					if (node_order[i] == min_node) {
						get_left_heavy_path(node_order, i, path_node);
						path_node.push_back(min_node);
					} else {
						if (node_order[i] == p){
							get_right_heavy_path(node_order, i, path_node);								
							break;
						}
					}
				}
				
				string transcript = node_set[path_node[0]].sequence;
				for (int i = 1; i < path_node.size(); ++i) {
					for (int j = 0; j < node_set[path_node[i-1]].children.size(); ++j) {
						if (node_set[path_node[i-1]].children[j].first == path_node[i]) {
							transcript = transcript + node_set[path_node[i]].sequence;
							break;
						}
					}
				}

				if (transcript.length() >= 200 && min_cov > 0.001*total_ave / connect_id)
					trans_vec.push_back(transcript);

				node_set[min_node].delete_child(p);
				node_set[p].delete_parent(min_node);
			}
			for (int i = 0; i < final_vec.size(); ++i) {
				is_out[final_vec[i]] = -1;
			}

			connect_id = connect_id - 1;
			current_id = current_id - 1;
		} else {
			/*
			if (final_vec.size() > 1) {
				int total_edge = 0;
				for (int i = 0; i < final_vec.size(); ++i) {
					total_edge = total_edge + node_set[final_vec[i]].children.size();
				}
			}
			*/
			vector<vector<int> > path_vec;
			vector<int> path_node;
			for (int i = 0; i < final_vec.size(); ++i) {
				p = final_vec[i];
				if (node_set[p].parents.size() == 0) {
					path_node.clear();
					get_all_path(p, path_node, path_vec);
				}
			}
			//output graphs
			string par_file_name = out_dir + "SG/" + to_string(connect_id)+".par";
			string node_file_name = out_dir + "SG/" + to_string(connect_id) + ".node";
			fstream par_file, node_file;
			par_file.open(par_file_name.c_str(), fstream::out);
			node_file.open(node_file_name.c_str(), fstream::out);
			const int final_size = final_vec.size();
			for (int i = 0; i <final_size; ++i) {
				node_file << ">" << i << "\t" << node_set[final_vec[i]].coverage << "\t" + node_set[final_vec[i]].s_name << endl;
				node_file << node_set[final_vec[i]].sequence << endl;
			}
			node_file.close();
			vector<vector<double> > adj_matrix;
			for (int i = 0; i < final_size; ++i) {
				vector<double> temp(final_size, 0.0);
				adj_matrix.push_back(temp);
			}
			for (int i = 0; i < final_size; ++i) {
				adj_matrix[i][i] = node_set[final_vec[i]].coverage;
			}

			for (int i = 0; i < final_size; ++i) {
				for (int j = 0; j < node_set[final_vec[i]].children.size(); ++j) {
					int c = node_set[final_vec[i]].children[j].first;
					adj_matrix[i][comp_map[c]] = node_set[final_vec[i]].children[j].second;
					adj_matrix[comp_map[c]][i] = node_set[final_vec[i]].children[j].second;
				}
			}

			vector<vector<int> > out_path;
			for (int i = 0; i < path_vec.size(); ++i) {
				vector<int> temp(final_size, 0);
				out_path.push_back(temp);
			}

			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					out_path[i][comp_map[path_vec[i][j]]] = 1;
				}
			}

			for (int i = 0; i < final_size; ++i){
				for (int j = 0; j < final_size; ++j) {
					par_file << adj_matrix[i][j] << ",";
				}
				par_file << ";";
			}
			par_file << endl;

			map<vector<int>, bool>::iterator p_it;
			for (p_it = pair_edges.begin(); p_it != pair_edges.end(); ++p_it) {
				vector<int> pair_vec = p_it -> first;
				for (int i = 0; i < pair_vec.size(); ++i)
					par_file << pair_vec[i] << ",";
				par_file << ";";
			}
			par_file << endl;

			for (int i = 0; i < out_path.size(); ++i) {
				for (int j = 0; j < final_size; ++j)
					par_file << out_path[i][j] << ",";
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					par_file << comp_map[path_vec[i][j]] << ",";
				}
				par_file << ";";
			}
			par_file << endl;
			current_min = node_set[final_vec[0]].coverage;
			current_ave = 0.0;
			for (int i = 0; i < final_vec.size(); ++i) {
				if (node_set[final_vec[i]].coverage < current_min) 
					current_min = node_set[final_vec[i]].coverage;
				current_ave = current_ave + node_set[final_vec[i]].coverage;
			}

			total_ave = total_ave + current_ave / final_vec.size();

			par_file << current_min;

			par_file.close();
			
		}

		connect_id = connect_id + 1;

	}


}

bool Graph::has_circle(vector<int>& comp_vec) {

	for (int i = 0; i < comp_vec.size(); ++i) {
		for (int j = 0; j < node_set[comp_vec[i]].children.size(); ++j) {
			set<int> checked;
			int c = node_set[comp_vec[i]].children[j].first;
			if (is_circle(comp_vec[i], c, checked))
				return true;

		}
	}

	return false;

}



bool Graph::is_circle(int p, int q, set<int>& checked, vector<int>& c_path) {

	if (node_set[q].children.empty() || checked.find(q) != checked.end()) // case2: contain circle but the initial p,q is not in the circle
		return false;
	else
		checked.insert(q);

	bool flag = false;

	for (int it = 0; it <node_set[q].children.size();++it) {
		if (node_set[q].children[it].first == p) {
			c_path.push_back(q);
			return true;
		} else {
			flag = is_circle(p, node_set[q].children[it].first, checked, c_path);
			if (flag) {
				c_path.push_back(q);
				break;
			}
		}

	}

	return flag;

}



bool Graph::is_circle(int p, int q, set<int>& checked) {


	if (node_set[q].children.empty() || checked.find(q) != checked.end())
		return false;
	else
		checked.insert(q);

	bool flag = false;

	vector<pair<int, double> >::iterator it;
	for (it = node_set[q].children.begin(); it != node_set[q].children.end(); it++) {

		if ((*it).first == p) {
			flag = true;
			break;
		} else {
			flag = is_circle(p, (*it).first, checked);
			if (flag) {
				break;
			}
		}

	}

	return flag;

}


void Graph::add_to_comp(vector<int>& is_out, vector<int>& comp_vec, int current_id, int connect_id) {

	if (is_out[current_id] == connect_id)
		return;
	comp_vec.push_back(current_id);
	is_out[current_id] = connect_id;

	int p;
	for (int i = 0; i < node_set[current_id].parents.size(); ++i) {
		p = node_set[current_id].parents[i].first;
		add_to_comp(is_out, comp_vec, p, connect_id);
	}
	for (int i = 0; i < node_set[current_id].children.size(); ++i) {
		p = node_set[current_id].children[i].first;
		add_to_comp(is_out, comp_vec, p, connect_id);
	}
}



void Graph::topological_sort(map<int, int>& comp_map, vector<int>& node_order) {

	//cout << "topological sort..." << endl;
	const int comp_size = comp_map.size();
	vector<int> node_color(comp_size, 0);
	
	map<int, int>::iterator it;
	for (it = comp_map.begin(); it != comp_map.end(); ++it) {
		int p = it -> first;
		if (node_set[p].parents.size() == 0) {
			dfs_visit(p, comp_map, node_color, node_order);

		}

	}

	//cout << "done" << endl;

}


void Graph::dfs_visit(int i, map<int, int>& comp_map, vector<int>& node_color, vector<int>& node_order) {

	node_color[comp_map[i]] = 1;

	if (node_set[i].children.size() == 0)
		node_order.push_back(i);
	else {
		for (int j  = 0; j < node_set[i].children.size(); ++j) {
			int p = node_set[i].children[j].first;
			if (node_color[comp_map[p]] == 0)
				dfs_visit(p, comp_map, node_color, node_order);
		}
		node_order.push_back(i);
	}

}

void Graph::node_path_sum(map<int, int>& comp_map, vector<int>& sum_vec, int p) {

	int i = comp_map[p];

	if (node_set[p].parents.size() == 0)
		sum_vec[i] = 1;
	else {

		for (int j = 0; j < node_set[p].parents.size(); ++j) {
			int k = node_set[p].parents[j].first;
			sum_vec[i] = sum_vec[i] + sum_vec[comp_map[k]];
			if (sum_vec[i] > 100)
				return;
		}
	}
}


int Graph::total_path_num(map<int, int>& comp_map, vector<int>& node_order) {
//cout << "total path num " << endl;
	int total_num = 0;
	int comp_size = comp_map.size();
	vector<int> sum_vec(comp_size, 0);

	for (int i = node_order.size()-1; i >= 0; --i) {
		node_path_sum(comp_map, sum_vec, node_order[i]);
	}
	
	for (int i = 0; i < node_order.size(); ++i) {
		int p = node_order[i];
		if (node_set[p].children.size() == 0){
			total_num = total_num + sum_vec[comp_map[p]];
		}
	}
//cout << "done" << endl;
	return total_num;

}


void Graph::get_all_path(int node, vector<int>& path_node, vector<vector<int> >& path_vec) { 
//cout << "get all path" << endl;
	path_node.push_back(node);
	if (node_set[node].children.size() == 0) {
		path_vec.push_back(path_node);
	} else {
		for (int i = 0; i < node_set[node].children.size(); ++i) {
			int c = node_set[node].children[i].first;
			get_all_path(c, path_node, path_vec);
		}
	}
	path_node.pop_back();
//cout << "done" << endl;
}


void Graph::get_left_heavy_path(std::vector<int>& node_order, int pos, std::vector<int>& path_node) {
	
	//cout << "get left heavy path" << endl;
	int p = node_order[pos];
	if (node_set[p].parents.size() == 0)
		return;
	const int graph_size = node_set.size();
	std::vector<int> parent_vec(graph_size, -1);
	std::vector<double> max_vec(graph_size, 0.0);
	for (int i = node_order.size()-1; i >= pos; --i) {
		p = node_order[i];
		max_vec[p] = node_set[p].coverage;
	}
	for (int i = node_order.size()-1; i >= pos; --i) {
		p = node_order[i];
		if (node_set[p].parents.size() == 0)
			continue;
		parent_vec[p] = node_set[p].parents[0].first;
		double max_parent = max_vec[parent_vec[p]];
		for (int j = 1; j < node_set[p].parents.size(); ++j) {
			int x = node_set[p].parents[j].first;
			if (max_vec[x] > max_parent) {
				parent_vec[p] = x;
				max_parent = max_vec[x];
			}
		}
		if (max_parent > 0 && max_parent < max_vec[p]) {
			max_vec[p] = max_parent;
		}
	}

	std::vector<int> temp;
	p = node_order[pos];
	while(parent_vec[p] >= 0) {
		temp.push_back(parent_vec[p]);
		p = parent_vec[p];
	}

	for (int i = temp.size()-1; i >= 0; --i) 
		path_node.push_back(temp[i]);

	//cout << "done" << endl;
}


void Graph::get_right_heavy_path(std::vector<int>& node_order, int pos, std::vector<int>& path_node) {

	//cout << "get right heavy path" << endl;
	int p = node_order[pos];
	if (node_set[p].children.size() == 0) {
		path_node.push_back(p);
		return;
	}			

	const int graph_size = node_set.size();
	std::vector<int> parent_vec(graph_size, -1);
	std::vector<double> max_vec(graph_size, 0.0);

	for (int i = pos; i >= 0; --i) {
		p = node_order[i];
		max_vec[p] = node_set[p].coverage;
	}

	for (int i = pos-1; i >= 0; --i) {
		p = node_order[i];
		if (node_set[p].parents.size() == 0) {
			max_vec[p] = 0;
			continue;
		}				
		parent_vec[p] = node_set[p].parents[0].first;
		double max_parent = max_vec[parent_vec[p]];
		for (int j = 1; j < node_set[p].parents.size(); ++j) {
			int x = node_set[p].parents[j].first;
			if (max_vec[x] > max_parent) {
				parent_vec[p] = x;
				max_parent = max_vec[0];
			}
		}
		if (max_parent == 0.0) {
			parent_vec[p] = -1;
		} else if (max_parent < max_vec[p]) {
			max_vec[p] = max_parent;
		}
	}

	p = node_order[0];
	for (int i = pos-1; i >= 0; --i) {
		int x = node_order[i];
		if (node_set[x].children.size() == 0) {
			if (max_vec[p] < max_vec[x])
				p = x;
		}			
	}

	std::vector<int> temp;
	temp.push_back(p);
	while(parent_vec[p] >= 0) {
		temp.push_back(parent_vec[p]);
		p = parent_vec[p];
	}

	for (int i = temp.size()-1; i >= 0; --i) 
		path_node.push_back(temp[i]);	
	
	// << "done" << endl;
}

