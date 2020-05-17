
#include "AssemblyGraph.h"
#include<iostream>

using namespace std;

string g_graph_file = "";
string g_out_dir = "";
string g_mapped_file = "";
string g_edge_file = "";
int  g_max_path_sum = 500;
bool g_help = false;

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

		vector<pair<int, int> >::iterator it;

		for (it = node_set[i].children.begin(); it != node_set[i].children.end(); ++it)
			node_set[it->first].add_parent(make_pair(i,it->second));

	}

}



int Graph::add_node(Node& node) {

	node_set.push_back(node);
	return (node_sum++);

}


void Graph::add_reverse_nodes() {

	mid_node_id = node_set.size();
	for (int i = 0; i < mid_node_id; ++i) {
		string node_seq = node_set[i].sequence;
		Node node;
		node.sequence = revcomp(node_seq);
		node.coverage = node_set[i].coverage;
		node.s_name = node_set[i].s_name.substr(0, node_set[i].s_name.length()-1)+"-";
		node.reads_sum = 0;
		node_set.push_back(node);
	}
	
	node_sum = node_set.size();
}



void Graph::load_graph(const string& graph_path, const string& mapped_file) {
cout << "load graph... " << endl;
	string s;
	string node_name, edge_name;
	map<string, int> node_map;
	double cov;
	fstream in2(graph_path);
	cout << "dealing file: " << graph_path << endl;
	getline(in2,s);
	bool is_first_edge = true;
	int kmer_length = 0;
	int mis_sum = 0;
	int max_overlap_len = 0;
	string pre_end = "";
	string node_seq;
	while (s.length() > 1) {
		
		if (s[0] == 'S') {
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			node_name = s.substr(s1+1, s2-s1-1);
			node_seq = s.substr(s2+1, s3-s2-1);
			Node node;
			node_map[node_name] = node_set.size();
			cov = atof(s.substr(s3+6).c_str());
			node.sequence = node_seq;
			node.coverage = cov;
			node.s_name = node_name + "+";
			node.reads_sum = 0;
			node_set.push_back(node);
			
			
		} else if (s[0] == 'L') {
			if (is_first_edge) {
				is_first_edge = false;
				add_reverse_nodes();
			}
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string::size_type s4 = s.substr(s3+1).find("\t") + s3 + 1;
			string::size_type s5 = s.substr(s4+1).find("\t") + s4 + 1;
			
			string map_info = s.substr(s5 + 1);
			if (map_info.substr(map_info.length()-1) != "M" || is_not_num(map_info.substr(0,map_info.length()-1))) {
				getline(in2,s);
				continue;
			}


			string stand1 = s.substr(s2+1, 1);
			string stand2 = s.substr(s4+1, 1);
			string name1 = s.substr(s1+1, s2-s1-1);
			string name2 = s.substr(s3+1, s4-s3-1);
			int n1 = node_map[name1];
			int n2 = node_map[name2];
			int rev_n1, rev_n2;
			if (stand1 == "-") {
				rev_n1 = n1;
				n1 = n1 + mid_node_id;
			} else {
				rev_n1 = mid_node_id + n1;
			}	
			
			if (stand2 == "-") {
				rev_n2 = n2;
				n2 = n2 + mid_node_id;
			} else {
				rev_n2 = mid_node_id + n2;
			}
			if (n1 == n2) {	
				int cnid = node_map[name1];
				circle_node[cnid] = true;
				mis_sum = mis_sum + 1;
				getline(in2,s);
				continue;
			}
			int overlap_len = atoi(map_info.substr(0,map_info.length()-1).c_str());
			if (overlap_len > max_overlap_len)
				max_overlap_len = overlap_len;
			
			node_set[n1].add_child(make_pair(n2, overlap_len));
			node_set[rev_n2].add_child(make_pair(rev_n1, overlap_len));
			node_set[n2].add_parent(make_pair(n1, overlap_len));
			node_set[rev_n1].add_parent(make_pair(rev_n2, overlap_len));
		}else if (s[0] == 'P') { // check ?????
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string name1 = s.substr(s1+1, s2-s1-1);
			string node_str = s.substr(s2+1, s3-s2-1);
			string::size_type s4 = node_str.find(",");
			
			if (name1[name1.length()-1] != '1') {
				string child_node = "";
				if (s4 == string::npos) {
					child_node = node_str;
				} else {
					child_node = node_str.substr(0,s4);
				}
				int child = node_map[child_node.substr(0, child_node.length()-1)];
				if (child_node[child_node.length()-1] == '-') 
					child = child + mid_node_id;
				int parent = node_map[pre_end.substr(0,pre_end.length()-1)];
				if (pre_end[pre_end.length()-1] == '-')
					parent = parent + mid_node_id;
				//if (node_set[parent].children.size() == 0 && node_set[child].parents.size() == 0)
				add_path_edge(parent, child, max_overlap_len);
				
			}
			pre_end = node_str;
			while (s4 != string::npos) {
				pre_end = pre_end.substr(s4+1);
				s4 = pre_end.find(",");
			}
			
		}

		getline(in2,s);

	}

	kmer_length = max_overlap_len - 1;
	for (int i = 0; i < mid_node_id; ++i) {
		if (node_set[i].sequence.length() <= kmer_length)
			node_set[i].coverage = 1.0*node_set[i].coverage / kmer_length;
		else
			node_set[i].coverage = 1.0*node_set[i].coverage / (node_set[i].sequence.length()-kmer_length);
		node_set[i+mid_node_id].coverage = node_set[i].coverage;		
	}
	cout << "node_sum: " << node_set.size() / 2 << endl;
	int xyz = 0;
	for (int i = 0; i < node_set.size(); ++i) 
		xyz = xyz + node_set[i].children.size();
	cout << "edge_sum: " << xyz / 2 << endl;
	cout << "dealing file: " << mapped_file << endl;
	fstream in3(mapped_file);
	getline(in3,s);
	while(s.length() > 0 && s[0] == '@'){
		getline(in3,s);
	}

	string pre_read = "";
	string read_str, node_str;
	int pre_node = -1;
	string pre_s = "";
	while (s.length() > 9) {
		string::size_type s1 = s.find("\t");
		read_str = s.substr(0, s1);
		string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
		string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
		node_str = s.substr(s2+1, s3-s2-1);
		string::size_type s4 = s.substr(s3+1).find("\t") + s3 + 1;
		string::size_type s5 = s.substr(s4+1).find("\t") + s4 + 1;
		string::size_type s6 = s.substr(s5+1).find("\t") + s5 + 1;
		string mapq = s.substr(s5+1, s6-s5-1);
		int n1 = -1;
		if (node_str != "*" && mapq.find("I") == string::npos && mapq.find("D") == string::npos && mapq.find("N") == string::npos &&
		mapq.find("S")==string::npos && mapq.find("H") == string::npos && mapq.find("P")==string::npos &&
		mapq.find("=")==string::npos && mapq.find("X")==string::npos) {
			n1 = atoi(node_str.c_str());
			int rev_n = 0;
			if (n1 >= mid_node_id) {
				rev_n = n1 - mid_node_id;
			} else {
				rev_n = n1 + mid_node_id;
			}

			node_set[n1].reads_sum = node_set[n1].reads_sum + 1;
			if (pre_read == read_str && pre_node != n1 && n1 != pre_node && pre_node >= 0) {				
				node_set[n1].add_pair(make_pair(pre_node, 1));
				node_set[pre_node].add_pair(make_pair(n1, 1));
				int rev_pre = 0;
				if (pre_node < mid_node_id)
					rev_pre = pre_node + mid_node_id;
				else
					rev_pre = pre_node - mid_node_id;
				node_set[rev_n].add_pair(make_pair(rev_pre, 1));
				node_set[rev_pre].add_pair(make_pair(rev_n, 1));
			}					
		}
		pre_s = s;
		pre_node = n1;	
		pre_read = read_str;	
		getline(in3,s);
	}

	int total_pair = 0;
	for (int i = 0; i < node_set.size(); ++i) {
		total_pair = total_pair + node_set[i].pair_nodes.size();
	}


}


void Graph::load_graph(const string& graph_path, const string& edge_file, const string& mapped_file) {
cout << "load graph... " << endl;
	string s;
	fstream in(edge_file);
	string node_name, edge_name;
	map<string, int> node_map;
	map<string, int> edge_map;
	double cov;
cout << "dealing file: " << edge_file << endl;
	getline(in,s);
	while (s.length() > 1) {
		string::size_type s1 = s.find("EDGE_");
		string::size_type s2 = s.find("_length_");
		string::size_type s3 = s.find("_cov_");
		string::size_type s4 = s.find("\'");
		string::size_type s5 = s.find(";");
		string::size_type s6 = s.find(":");
		string::size_type s7 = s.find(",");
		edge_map[s.substr(1,s.length()-1)] = node_set.size();
		bool is_continue = false;
		if (s4 != string::npos) { //exist reverse
			if (s5 != string::npos && s6 == string::npos && s7 == string::npos) //single edge
				is_continue = true;
			if (s6 != string::npos && s7 == string::npos && s4 < s6) // multi edge ：
				is_continue = true;
			if (s6 == string::npos && s7 != string::npos && s4 < s7) // multi edge ,
				is_continue = true;
			if (s6 != string::npos && s7 != string::npos && s4 < s6 && s6 < s7) // multi edge ：，
				is_continue = true;
			if (s6 != string::npos && s7 != string::npos && s4 < s7 && s7 < s6)
				is_continue = true;
		}

		if (is_continue) {
			getline(in,s);
			while (s.length() > 0 && s[0] != '>') { // the first node is reverse
				getline(in,s);	
			}
			continue;
		}

		node_name = s.substr(s1+5,s2-s1-5) + "+";
		cov = atof(s.substr(s3+5, s5-s3-5).c_str());
		if (s6 == string::npos && s7 == string::npos) {	//single edge		
			Node node;
			node.s_name = node_name;
			node.coverage = cov;
			node_map[node_name.substr(0,node_name.length()-1)] = node_set.size();
			getline(in,s);		
			string trans_seq = "";
			while (s.length() > 0 && s[0] != '>') {
				trans_seq = trans_seq + s;
				getline(in,s);	
			}		
			node.sequence = trans_seq;
			node.reads_sum = 0;
			node_set.push_back(node);
		} else {
			int len = atoi(s.substr(s2+8,s3-s2-8).c_str());
			int total_len = len;
			double total_cov = cov*len;
			if (s6 != string::npos && s7 == string::npos) {
				s = s.substr(s6+1);
				s1 = s.find("EDGE_");
				while (s1 != string::npos) {
					s2 = s.find("_length_");
					s3 = s.find("_cov_");
					s4 = s.find("\'");
					s5 = s.find(";");
					s6 = s.find(":");
					node_name = node_name + "," + s.substr(s1+5,s2-s1-5);
					string stand = "+";
					if (s4 != string::npos) {
						if (s6 == string::npos)
							stand = "-";
						else {
							if (s4 < s6)
								stand = "-";
						}
					}
					node_name = node_name + stand;
					len = atoi(s.substr(s2+8,s3-s2-8).c_str());
					cov = atof(s.substr(s3+5, s5-s3-5).c_str());
					total_cov = total_cov + len*cov;
					total_len = total_len + len;
					if (s6 == string::npos) {
						break;
					} else {
						s = s.substr(s6+1);
						s1 = s.find("EDGE_");
					}
				}
				Node node;
				node.s_name = node_name;
				node.coverage = total_cov*1.0/total_len;
				node_map[node_name] = node_set.size();
				getline(in,s);		
				string trans_seq = "";
				while (s.length() > 0 && s[0] != '>') {
					trans_seq = trans_seq + s;
					getline(in,s);	
				}		
				node.sequence = trans_seq;
				node.reads_sum = 0;
				node_set.push_back(node);
			} // all :
			if (s6 == string::npos && s7 != string::npos) {
				s = s.substr(s7+1);
				s1 = s.find("EDGE_");
				while (s1 != string::npos) {
					s2 = s.find("_length_");
					s3 = s.find("_cov_");
					s4 = s.find("\'");
					s5 = s.find(";");
					s7 = s.find(",");
					node_name = node_name + "," + s.substr(s1+5,s2-s1-5);
					string stand = "+";
					if (s4 != string::npos) {
						if (s7 == string::npos)
							stand = "-";
						else {
							if (s4 < s7)
								stand = "-";
						}
					}
					node_name = node_name + stand;
					len = atoi(s.substr(s2+8,s3-s2-8).c_str());
					cov = atof(s.substr(s3+5, s5-s3-5).c_str());
					total_cov = total_cov + len*cov;
					total_len = total_len + len;
					if (s7 == string::npos) {
						break;
					} else {
						s = s.substr(s7+1);
						s1 = s.find("EDGE_");
					}
				}
				Node node;
				node.s_name = node_name;
				node.coverage = total_cov*1.0/total_len;
				node_map[node_name] = node_set.size();
				getline(in,s);		
				string trans_seq = "";
				while (s.length() > 0 && s[0] != '>') {
					trans_seq = trans_seq + s;
					getline(in,s);	
				}		
				node.sequence = trans_seq;
				node.reads_sum = 0;
				node_set.push_back(node);
			} //all ,
			if (s6 != string::npos && s7 != string::npos) {
				if (s6 < s7) {
					s = s.substr(s7+1);
				} else {
					s = s.substr(s7+1);
				}
				s1 = s.find("EDGE_");
				while (s1 != string::npos) {
					s2 = s.find("_length_");
					s3 = s.find("_cov_");
					s4 = s.find("\'");
					s5 = s.find(";");
					s6 = s.find(":");
					s7 = s.find(",");
					node_name = node_name + "," + s.substr(s1+5,s2-s1-5);
					string stand = "+";
					if (s4 != string::npos) {
						if (s7 == string::npos && s6 == string::npos)
							stand = "-";
						else {
							if (s7 != string::npos && s6 == string::npos && s4 < s7)
								stand = "-";
							if (s7 == string::npos && s6 != string::npos && s4 < s6)
								stand = "-";
							if (s7 != string::npos && s6 != string::npos && s4 < s6 && s6 < s7)
								stand = "-";
							if (s7 != string::npos && s6 != string::npos && s4 < s7 && s7 < s6)
								stand = "-";
						}
					}
					node_name = node_name + stand;
					len = atoi(s.substr(s2+8,s3-s2-8).c_str());
					cov = atof(s.substr(s3+5, s5-s3-5).c_str());
					total_cov = total_cov + len*cov;
					total_len = total_len + len;
					if (s7 == string::npos && s6 == string::npos) {
						break;
					} else {
						if (s7 == string::npos)
							s = s.substr(s6+1);
						if (s6 == string::npos)
							s = s.substr(s7+1);
						if (s6 != string::npos && s7 != string::npos) {
							if (s6 > s7)
								s = s.substr(s7+1);
							else
								s = s.substr(s6+1);
						}
						s1 = s.find("EDGE_");
					}
				}
				Node node;
				node.s_name = node_name;
				node.coverage = total_cov*1.0/total_len;
				node_map[node_name] = node_set.size();
				getline(in,s);		
				string trans_seq = "";
				while (s.length() > 0 && s[0] != '>') {
					trans_seq = trans_seq + s;
					getline(in,s);	
				}		
				node.sequence = trans_seq;
				node.reads_sum = 0;
				node_set.push_back(node);
			} //: and ,
		} //multi edges		
		
	} //while (s.length() > 1)
cout << "edge sum: " << node_set.size() <<endl;
cout << "dealing file: " <<graph_path << endl;
	fstream in2(graph_path);
	getline(in2,s);
	bool is_first_edge = true;
	int kmer_length = 0;
	int mis_sum = 0;
	int max_overlap_len = 0;
	string pre_end = "";
	string node_seq;
	while (s.length() > 1) {
		
		if (s[0] == 'S') {
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			node_name = s.substr(s1+1, s2-s1-1);
			if (node_map.find(node_name) == node_map.end()) {
				node_seq = s.substr(s2+1, s3-s2-1);
				Node node;
				node_map[node_name] = node_set.size();
				cov = atof(s.substr(s3+6).c_str());
				node.sequence = node_seq;
				node.coverage = cov;
				node.s_name = node_name + "+";
				node.reads_sum = 0;
				node_set.push_back(node);
			}
			
		} else if (s[0] == 'L') {
			if (is_first_edge) {
				is_first_edge = false;
				add_reverse_nodes();
			}
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string::size_type s4 = s.substr(s3+1).find("\t") + s3 + 1;
			string::size_type s5 = s.substr(s4+1).find("\t") + s4 + 1;
			
			string map_info = s.substr(s5 + 1);
			if (map_info.substr(map_info.length()-1) != "M" || is_not_num(map_info.substr(0,map_info.length()-1))) {
				getline(in,s);
				continue;
			}


			string stand1 = s.substr(s2+1, 1);
			string stand2 = s.substr(s4+1, 1);
			string name1 = s.substr(s1+1, s2-s1-1);
			string name2 = s.substr(s3+1, s4-s3-1);
			int n1 = node_map[name1];
			int n2 = node_map[name2];
			int rev_n1, rev_n2;
			if (stand1 == "-") {
				rev_n1 = n1;
				n1 = n1 + mid_node_id;
			} else {
				rev_n1 = mid_node_id + n1;
			}	
			
			if (stand2 == "-") {
				rev_n2 = n2;
				n2 = n2 + mid_node_id;
			} else {
				rev_n2 = mid_node_id + n2;
			}
			if (n1 == n2) {	
				int cnid = node_map[name1];
				circle_node[cnid] = true;
				mis_sum = mis_sum + 1;
				getline(in2,s);
				continue;
			}
			int overlap_len = atoi(map_info.substr(0,map_info.length()-1).c_str());
			if (overlap_len > max_overlap_len)
				max_overlap_len = overlap_len;
			
			node_set[n1].add_child(make_pair(n2, overlap_len));
			node_set[rev_n2].add_child(make_pair(rev_n1, overlap_len));
			node_set[n2].add_parent(make_pair(n1, overlap_len));
			node_set[rev_n1].add_parent(make_pair(rev_n2, overlap_len));
		}else if (s[0] == 'P') { // check ?????
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string name1 = s.substr(s1+1, s2-s1-1);
			string node_str = s.substr(s2+1, s3-s2-1);
			string::size_type s4 = node_str.find(",");
			
			if (name1[name1.length()-1] != '1') {
				string child_node = "";
				if (s4 == string::npos) {
					child_node = node_str;
				} else {
					child_node = node_str.substr(0,s4);
				}
				int child = node_map[child_node.substr(0, child_node.length()-1)];
				if (child_node[child_node.length()-1] == '-') 
					child = child + mid_node_id;
				int parent = node_map[pre_end.substr(0,pre_end.length()-1)];
				if (pre_end[pre_end.length()-1] == '-')
					parent = parent + mid_node_id;
				//if (node_set[parent].children.size() == 0 && node_set[child].parents.size() == 0)
				add_path_edge(parent, child, max_overlap_len);
				
			}
			pre_end = node_str;
			while (s4 != string::npos) {
				pre_end = pre_end.substr(s4+1);
				s4 = pre_end.find(",");
			}
			
		}

		getline(in2,s);

	}

	kmer_length = max_overlap_len - 1;
	for (int i = 0; i < mid_node_id; ++i) {
		if (node_set[i].sequence.length() <= kmer_length)
			node_set[i].coverage = 1.0*node_set[i].coverage / kmer_length;
		else
			node_set[i].coverage = 1.0*node_set[i].coverage / (node_set[i].sequence.length()-kmer_length);
		node_set[i+mid_node_id].coverage = node_set[i].coverage;		
	}
	cout << "edge_sum: " << node_set.size() / 2 << endl;
	cout << "dealing file: " << mapped_file << endl;
	fstream in3(mapped_file);
	getline(in3,s);
	while(s.length() > 0 && s[0] == '@') 
		getline(in3,s);

	string pre_read = "";
	string read_str, edge_str;
	int pre_node = -1;

	while (s.length() > 1) {
		//cout << s << endl;
		string::size_type s1 = s.find("\t");
		read_str = s.substr(0, s1);
		string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
		string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
		edge_str = s.substr(s2+1, s3-s2-1);
		int n1 = -1;
		if (edge_str.length() > 1) {
			n1 = edge_map[edge_str];
			node_set[n1].reads_sum = node_set[n1].reads_sum + 1;
			if (pre_read == read_str && pre_node != n1 && n1 >= 0 && pre_node >= 0) {
				node_set[n1].add_pair(make_pair(pre_node, 1));
				node_set[pre_node].add_pair(make_pair(n1, 1));
			}
			pre_read = read_str;
			pre_node = n1;
		}		
		getline(in3,s);
	}

	int total_pair = 0;
	for (int i = 0; i < node_set.size(); ++i) {
		total_pair = total_pair + node_set[i].pair_nodes.size();
	}


}


void Graph::get_non_trivial_graph_nodes(const string& graph_path, string out_dir) {

	string s;
	string node_name, edge_name;
	map<string, int> node_map;
	double cov;
	fstream in2(graph_path);
	getline(in2,s);
	bool is_first_edge = true;
	int kmer_length = 0;
	int mis_sum = 0;
	int max_overlap_len = 0;
	string pre_end = "";
	string node_seq;
	while (s.length() > 1) {
		
		if (s[0] == 'S') {
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			node_name = s.substr(s1+1, s2-s1-1);
			node_seq = s.substr(s2+1, s3-s2-1);
			Node node;
			node_map[node_name] = node_set.size();
			cov = atof(s.substr(s3+6).c_str());
			node.sequence = node_seq;
			node.coverage = cov;
			node.s_name = node_name + "+";
			node.reads_sum = 0;
			node_set.push_back(node);
			
			
		} else if (s[0] == 'L') {
			if (is_first_edge) {
				is_first_edge = false;
				add_reverse_nodes();
			}
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string::size_type s4 = s.substr(s3+1).find("\t") + s3 + 1;
			string::size_type s5 = s.substr(s4+1).find("\t") + s4 + 1;
			
			string map_info = s.substr(s5 + 1);
			if (map_info.substr(map_info.length()-1) != "M" || is_not_num(map_info.substr(0,map_info.length()-1))) {
				getline(in2,s);
				continue;
			}


			string stand1 = s.substr(s2+1, 1);
			string stand2 = s.substr(s4+1, 1);
			string name1 = s.substr(s1+1, s2-s1-1);
			string name2 = s.substr(s3+1, s4-s3-1);
			int n1 = node_map[name1];
			int n2 = node_map[name2];
			int rev_n1, rev_n2;
			if (stand1 == "-") {
				rev_n1 = n1;
				n1 = n1 + mid_node_id;
			} else {
				rev_n1 = mid_node_id + n1;
			}	
			
			if (stand2 == "-") {
				rev_n2 = n2;
				n2 = n2 + mid_node_id;
			} else {
				rev_n2 = mid_node_id + n2;
			}
			if (n1 == n2) {	
				int cnid = node_map[name1];
				circle_node[cnid] = true;
				mis_sum = mis_sum + 1;
				getline(in2,s);
				continue;
			}
			int overlap_len = atoi(map_info.substr(0,map_info.length()-1).c_str());
			if (overlap_len > max_overlap_len)
				max_overlap_len = overlap_len;
			
			node_set[n1].add_child(make_pair(n2, overlap_len));
			node_set[rev_n2].add_child(make_pair(rev_n1, overlap_len));
			node_set[n2].add_parent(make_pair(n1, overlap_len));
			node_set[rev_n1].add_parent(make_pair(rev_n2, overlap_len));
		}else if (s[0] == 'P') { // check ?????
			string::size_type s1 = s.find("\t");
			string::size_type s2 = s.substr(s1+1).find("\t") + s1 + 1;
			string::size_type s3 = s.substr(s2+1).find("\t") + s2 + 1;
			string name1 = s.substr(s1+1, s2-s1-1);
			string node_str = s.substr(s2+1, s3-s2-1);
			string::size_type s4 = node_str.find(",");
			
			if (name1[name1.length()-1] != '1') {
				string child_node = "";
				if (s4 == string::npos) {
					child_node = node_str;
				} else {
					child_node = node_str.substr(0,s4);
				}
				int child = node_map[child_node.substr(0, child_node.length()-1)];
				if (child_node[child_node.length()-1] == '-') 
					child = child + mid_node_id;
				int parent = node_map[pre_end.substr(0,pre_end.length()-1)];
				if (pre_end[pre_end.length()-1] == '-')
					parent = parent + mid_node_id;
				//if (node_set[parent].children.size() == 0 && node_set[child].parents.size() == 0)
				add_path_edge(parent, child, max_overlap_len);
				
			}
			pre_end = node_str;
			while (s4 != string::npos) {
				pre_end = pre_end.substr(s4+1);
				s4 = pre_end.find(",");
			}
			
		}

		getline(in2,s);

	}

	string long_file_name2 = out_dir + "nodeSeq.fa";
	fstream long_file2;
	long_file2.open(long_file_name2.c_str(), fstream::out);
	for (int i = 0; i < node_set.size(); ++i){
		if (node_set[i].parents.size() > 0 || node_set[i].children.size() > 0)
			long_file2 << ">" << i << endl << node_set[i].sequence << endl;
	}

}



void Graph::save_graph(string out_dir) {

cout << "preprocessing..." << endl;

	double total_ave = 0.0;
	int total_len = 0;
	for (int i = 0; i < mid_node_id; ++i) {
		total_ave = total_ave + node_set[i].coverage*node_set[i].sequence.length();
		total_len = total_len + node_set[i].sequence.length();
	}
	total_ave = total_ave / total_len;
	total_ave = total_ave * total_len / mid_node_id;
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
		for (int i = 0; i < comp_vec.size(); ++i) {
			if (comp_vec[i] < mid_node_id) 
				mate_id = comp_vec[i] + mid_node_id;
			else 
				mate_id = comp_vec[i] - mid_node_id;
			if (is_out[mate_id] == connect_id){
				is_conflict = true;
			}
			is_out[mate_id] = connect_id;
		}
//cout <<comp_vec.size() << " " << current_id  << " " << connect_id << endl;
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
//cout << path_sum << endl;
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
			set<int> first_node;
			for (int i = 0; i < comp_vec.size(); ++i) {
				if (comp_vec[i] < mid_node_id)
					first_node.insert(comp_vec[i]);
				else
					first_node.insert(comp_vec[i]-mid_node_id);
			}

			string par_file_name = out_dir + "SG/" + to_string(connect_id)+".par";
			string node_file_name = out_dir + "SG/" + to_string(connect_id) + ".node";
			fstream par_file, node_file;
			node_file.open(node_file_name.c_str(), fstream::out);

			set<int>::iterator it;
			for (it = first_node.begin(); it != first_node.end(); ++it) {
				node_file << ">" << node_set[*it].coverage << "\t" + node_set[*it].s_name << endl;
				node_file << node_set[*it].sequence << endl;
			}
			node_file.close();		
			const int first_size = first_node.size();
			vector<vector<double> > adj_matrix;
			for (int i = 0; i < first_size; ++i) {
				vector<double> temp(first_size, 0.0);
				adj_matrix.push_back(temp);
			}
			map<int, int> comp_map;			
			int idx = 0;
			for (it = first_node.begin(); it != first_node.end(); ++it) {
				adj_matrix[idx][idx] = node_set[*it].coverage;
				comp_map[*it] = idx;
				idx = idx + 1;
			}			
			for (int i = 0; i < comp_vec.size(); ++i) {
				int xi = comp_vec[i];
				if (xi >= mid_node_id)
					xi = comp_vec[i] - mid_node_id;
				for (int j = 0; j < node_set[comp_vec[i]].children.size(); ++j) {
					int c = node_set[comp_vec[i]].children[j].first;					
					int yi = c;					
					if (yi >= mid_node_id)
						yi = yi - mid_node_id;
					adj_matrix[comp_map[xi]][comp_map[yi]] = node_set[comp_vec[i]].children[j].second;
					adj_matrix[comp_map[yi]][comp_map[xi]] = node_set[comp_vec[i]].children[j].second;
				}
			}
			vector<vector<int> > out_path;
			for (int i = 0; i < path_vec.size(); ++i) {
				vector<int> temp(first_size, 0);
				out_path.push_back(temp);
			}

			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					int pv = path_vec[i][j];
					if (pv >= mid_node_id)
						pv = pv - mid_node_id;
					out_path[i][comp_map[pv]] = out_path[i][comp_map[pv]] + 1;
				}
			}
			par_file.open(par_file_name.c_str(), fstream::out);
			
			for (int i = 0; i < first_size; ++i){
				for (int j = 0; j < first_size; ++j) {
					par_file << adj_matrix[i][j] << ",";
				}
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < out_path.size(); ++i) {
				for (int j = 0; j < first_size; ++j) {
					par_file << out_path[i][j] << ",";
				}					
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					int pv = path_vec[i][j];
					if (pv >= mid_node_id)
						pv = pv - mid_node_id;
					par_file << comp_map[pv] << ",";
				}
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < path_vec.size(); ++i) {
				for (int j = 0; j < path_vec[i].size(); ++j) {
					int pv = path_vec[i][j];
					if (pv >= mid_node_id)
						par_file << 0 << ",";
					else
						par_file << 1 << ",";
				}
				par_file << ";";
			}
			par_file << endl;
			for (int i = 0; i < comp_vec.size(); ++i) {
				int xi = comp_vec[i];
				if (xi >= mid_node_id)
					xi = comp_vec[i] - mid_node_id;
				for (int j = 0; j < node_set[comp_vec[i]].pair_nodes.size(); ++j) {
					int yi = node_set[comp_vec[i]].pair_nodes[j].first;
					if (yi >= mid_node_id)
						yi = yi - mid_node_id;
					if (comp_map.find(yi) != comp_map.end())
						par_file << comp_map[xi] << "," << comp_map[yi] << ",;";
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
			int rev_p, rev_q;
			if (edge.first >= mid_node_id) 
				rev_p = edge.first - mid_node_id;
			else 
				rev_p = edge.first + mid_node_id;
			if (edge.second >= mid_node_id)
				rev_q = edge.second - mid_node_id;
			else
				rev_q = edge.second + mid_node_id;
			node_set[rev_q].delete_child(rev_p);
			node_set[rev_p].delete_parent(rev_q);

			/*if (min_score2 >= 0.2*total_ave || min_score2 > 4) {
				path_node.clear();
				greedy_path(edge, path_node, is_checked, comp_vec);
				string transcript = node_set[path_node[0]].sequence;
				string trans_name = node_set[path_node[0]].s_name;
				for (int i = 1; i < path_node.size(); ++i) {
					for (int j = 0; j < node_set[path_node[i-1]].children.size(); ++j) {
						if (node_set[path_node[i-1]].children[j].first == path_node[i]) {
							trans_name = trans_name + "," + node_set[path_node[i]].s_name;
							transcript = transcript + node_set[path_node[i]].sequence.substr(node_set[path_node[i-1]].children[j].second);
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
						transcript = transcript + node_set[path_node[i]].sequence.substr(node_set[path_node[i-1]].children[j].second);
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


int Graph::delete_circle(vector<int>& comp_vec) {
	//cout << "delete_circle" << endl;

	int dele_sum = 0;
	for (int i = 0; i < comp_vec.size(); ++i) {
		int j = 0;
		while (j < node_set[comp_vec[i]].children.size()) {

			set<int> checked;
			int c = node_set[comp_vec[i]].children[j].first;
			vector<int> c_path;
			if (is_circle(comp_vec[i], c, checked, c_path)) {
				c_path.push_back(comp_vec[i]);
				dele_sum = dele_sum + 1;
				if (c_path.size() == 2) {
					node_set[c_path[1]].delete_parent(c_path[0]);
					node_set[c_path[0]].delete_child(c_path[1]);
					continue;
				}
				int min_id = 0;
				double min_cov = node_set[c_path[0]].coverage;
				for (int k = 0; k < c_path.size(); ++k) {
					if (node_set[c_path[k]].coverage < min_cov) {
						min_cov = node_set[c_path[k]].coverage;
						min_id = k;
					}
				}
				if (min_id > 0 && min_id < c_path.size()-1) {
					if (node_set[c_path[min_id-1]].coverage > node_set[c_path[min_id+1]].coverage) {
						node_set[c_path[min_id]].delete_parent(c_path[min_id+1]);
						node_set[c_path[min_id+1]].delete_child(c_path[min_id]);
					} else {
						node_set[c_path[min_id]].delete_child(c_path[min_id-1]);
						node_set[c_path[min_id-1]].delete_parent(c_path[min_id]);
					}
				}
				if (min_id == 0) {
					if (node_set[c_path[c_path.size()-1]].coverage > node_set[c_path[min_id+1]].coverage) {
						node_set[c_path[min_id]].delete_parent(c_path[min_id+1]);
						node_set[c_path[min_id+1]].delete_child(c_path[min_id]);
					} else {
						node_set[c_path[min_id]].delete_child(c_path[c_path.size()-1]);
						node_set[c_path[c_path.size()-1]].delete_parent(c_path[min_id]);
					}
				}
				if (min_id == c_path.size()-1){ 
					if (node_set[c_path[min_id-1]].coverage > node_set[0].coverage) {
						node_set[c_path[min_id]].delete_parent(c_path[0]);
						node_set[c_path[0]].delete_child(c_path[min_id]);
					} else {
						node_set[c_path[min_id]].delete_child(c_path[min_id-1]);
						node_set[c_path[min_id-1]].delete_parent(c_path[min_id]);
					}
				}
			} else {
				j = j + 1;
			}

		}
	}

	//cout << "done" << endl;
	return dele_sum;
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

	vector<pair<int, int> >::iterator it;
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

	/*int mate_id = 0;
	if (current_id < mid_node_id)
		mate_id = current_id + mid_node_id;
	else
		mate_id = current_id - mid_node_id;
	add_to_comp(is_out, comp_vec, mate_id);
	*/
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

int Graph::divide_graph(vector<int>& comp_vec, vector<int>& is_out) {

//cout << "divide graph..." << endl;
	int mate_id;
	for(int i = 0; i < comp_vec.size(); ++i) {
		if (comp_vec[i] < mid_node_id)
			mate_id = comp_vec[i] + mid_node_id;
		else
			mate_id = comp_vec[i] - mid_node_id;
		if (is_out[comp_vec[i]] == -2 || is_out[mate_id] == -2) {
			is_out[comp_vec[i]] = -3;
			continue;
		}
		is_out[comp_vec[i]] = -2;
	}
//cout << "a" << endl;
	const int comp_size = is_out.size();
	vector<int> good_vec(comp_size, 0);
	vector<int> bad_vec(comp_size, 0);
	int p,q;
	for(int i = 0; i < comp_vec.size(); ++i) {
		p = comp_vec[i];
		if (is_out[p] != -2)
			continue;
		for (int j = 0; j < node_set[p].children.size(); ++j) {
			q = node_set[p].children[j].first;
			if (is_out[q] == -2) {
				good_vec[p] = good_vec[p] + 1;
				good_vec[q] = good_vec[q] + 1;
			} else {
				bad_vec[p] = bad_vec[p] + 1;
				if (q < mid_node_id)
					mate_id = q + mid_node_id;
				else
					mate_id = q - mid_node_id;
				bad_vec[q] = bad_vec[q] + 1;
			}
		}
		for (int j = 0; j < node_set[p].parents.size(); ++j) {
			q = node_set[p].parents[j].first;
			if (is_out[q] == -2) {
				good_vec[p] = good_vec[p] + 1;
				good_vec[q] = good_vec[q] + 1;
			} else {
				bad_vec[p] = bad_vec[p] + 1;
				if (q < mid_node_id)
					mate_id = q + mid_node_id;
				else
					mate_id = q - mid_node_id;
				bad_vec[q] = bad_vec[q] + 1;
			}
		}
	}

//cout << "b" << endl;
	int gap = 0;
	for (int i = 1; i < comp_vec.size();++i) {
		if (is_out[comp_vec[i]] == -2) {
			gap = bad_vec[comp_vec[i]] - good_vec[comp_vec[i]];
			p = comp_vec[i];
			break;
		}
	}
	for (int i = 1; i < comp_vec.size();++i) {
		if (is_out[comp_vec[i]] == -2 && (bad_vec[comp_vec[i]] - good_vec[comp_vec[i]] > gap)) {
			gap = bad_vec[comp_vec[i]] - good_vec[comp_vec[i]];
			p = comp_vec[i];
		}
	}
//cout << "c" << endl;
	int bad_sum = 0;
	for (int i = 1; i < comp_vec.size();++i) {
		if (is_out[comp_vec[i]] == -2 ) {
			bad_sum = bad_sum + bad_vec[comp_vec[i]];
		}
	}
	int current_bad = 0;
//cout << "d" << endl;
	while (gap > 0) {
		if (p < mid_node_id)
			mate_id = p + mid_node_id;
		else
			mate_id = p - mid_node_id;
		is_out[p] = -3;
		is_out[mate_id] = -2;
		good_vec[mate_id] = bad_vec[p];
		bad_vec[mate_id] = good_vec[p];
		bad_vec[p] = 0;
		good_vec[p] = 0;
		for (int j = 0; j < node_set[p].children.size(); ++j) {
			q = node_set[p].children[j].first;
			if (is_out[q] == -2) {
				good_vec[q] = good_vec[q] - 1;
				bad_vec[q] = bad_vec[q] + 1;
			} else {
				if (q < mid_node_id)
					mate_id = q + mid_node_id;
				else
					mate_id = q - mid_node_id;
				bad_vec[q] = bad_vec[q] - 1;
				good_vec[q] = good_vec[q] + 1;
			}
		}
		for (int j = 0; j < node_set[p].parents.size(); ++j) {
			q = node_set[p].parents[j].first;
			if (is_out[q] == -2) {
				good_vec[q] = good_vec[q] - 1;
				bad_vec[q] = bad_vec[q] + 1;
			} else {
				if (q < mid_node_id)
					mate_id = q + mid_node_id;
				else
					mate_id = q - mid_node_id;
				bad_vec[q] = bad_vec[q] - 1;
				good_vec[q] = good_vec[q] + 1;
			}
		}

		current_bad = 0;
		for (int i = 1; i < comp_vec.size();++i) {
			if (is_out[comp_vec[i]] == -2 ) {
				current_bad = current_bad + bad_vec[comp_vec[i]];
			}
		}
		if (current_bad >= bad_sum)
			break;
		else
			bad_sum = current_bad;
		gap = 0;
		for (int i = 1; i < comp_vec.size();++i) {
			if (is_out[comp_vec[i]] == -2 ) {
				gap = bad_vec[comp_vec[i]] - good_vec[comp_vec[i]];
				p = comp_vec[i];
				break;
			}
		}
		for (int i = 1; i < comp_vec.size();++i) {
			if (is_out[comp_vec[i]] == -2 && (bad_vec[comp_vec[i]] - good_vec[comp_vec[i]] > gap)) {
				gap = bad_vec[comp_vec[i]] - good_vec[comp_vec[i]];
				p = comp_vec[i];
			}
		}

	}



	for(int i = 0; i < comp_vec.size(); ++i) {
		p = comp_vec[i];
		if (is_out[p] != -2)
			continue;
		vector<int> dele_list;
		for (int j = 0; j < node_set[p].children.size(); ++j) {
			q = node_set[p].children[j].first;
			if (is_out[q] != -2) {
				dele_list.push_back(q);
			}
		}
		for (int j = 0; j < dele_list.size(); ++j) {
			node_set[p].delete_child(dele_list[j]);
			node_set[dele_list[j]].delete_parent(p);
		}
		dele_list.clear();
		for (int j = 0; j < node_set[p].parents.size(); ++j) {
			q = node_set[p].parents[j].first;
			if (is_out[q] != -2) {
				dele_list.push_back(q);
			}
		}
		for (int j = 0; j < dele_list.size(); ++j) {
			node_set[p].delete_parent(dele_list[j]);
			node_set[dele_list[j]].delete_child(p);
		}
	}
//cout << "done" << endl;
	return bad_sum;

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
			if (sum_vec[i] > g_max_path_sum)
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


void Graph::add_path_edge(int p, int q, int max_overlap) {

	int check_len = max_overlap;
	if (node_set[p].sequence.length() < check_len)
		check_len = node_set[p].sequence.length();
	if (node_set[q].sequence.length() < check_len)
		check_len = node_set[q].sequence.length();

	bool is_aligned = false;

	while (check_len > 0) {
		int match = 0;
		int mismatch = 0;
		for (int i = 0; i < check_len; ++i) {
			if (node_set[p].sequence[node_set[p].sequence.length()-check_len+i] == node_set[q].sequence[i]) 
				match++;
			else
				mismatch++;
		}
		if ((float)mismatch / check_len < 0.35) {
			node_set[p].add_child(make_pair(q, check_len));
			node_set[q].add_parent(make_pair(p,check_len));
			int rev_p = 0;
			int rev_q = 0;
			if (p < mid_node_id)
				rev_p = p + mid_node_id;
			else
				rev_p = p - mid_node_id;
			if (q < mid_node_id)
				rev_q = q + mid_node_id;
			else
				rev_q = q - mid_node_id;
			node_set[rev_q].add_child(make_pair(rev_p, check_len));
			node_set[rev_p].add_parent(make_pair(rev_q, check_len));
			break;
		}
		check_len = check_len - 1;
	}

	if (check_len == 0) {
		node_set[p].add_child(make_pair(q,0));
		node_set[q].add_parent(make_pair(p,0));
		int rev_p = 0;
		int rev_q = 0;
		if (p < mid_node_id)
			rev_p = p + mid_node_id;
		else
			rev_p = p - mid_node_id;
		if (q < mid_node_id)
			rev_q = q + mid_node_id;
		else
			rev_q = q - mid_node_id;
		node_set[rev_q].add_child(make_pair(rev_p, check_len));
		node_set[rev_p].add_parent(make_pair(rev_q, check_len));
	}

}
