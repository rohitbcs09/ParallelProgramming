#include <iostream>
#include <vector>
#include <map>
#include <limits>
#include <iostream>                                                                                                                 
#include <vector>                                                               
#include <unordered_map>                                                        
#include <set> 
#include "utils.h"

using namespace std;

int num_edges = 0;
int num_vertices = 0;
int heapSize = 0;

std::set<int> vertices;	

long int mst_weight = 0;	

int mst_edge_count = 0;

unordered_map<int, int> parent;

std::unordered_map<int, std::unordered_map<int, int> > Adj_list;

class UnionFind {
public:
    int setSize;
    std::unordered_map<int, int> rank;
    std::unordered_map<int, int> parent; 

    UnionFind(std::set<int>& vertex) {
    	setSize = vertex.size();
    	for(auto iter = vertex.begin(); iter != vertex.end(); ++iter) {
    		parent[*iter] = *iter;
            rank[*iter] = 1;
    	}
    }
    
    ~UnionFind() {

    }

    int find_set(int v)	{
    	int root = parent[v];
    	while(root != parent[root]) {
    		root = parent[root];
    	}

    	int p = parent[v];
    	while(p != root) {
    		int newpar = parent[p];
    		parent[p] = root;
    		p = newpar;
    	}
        return root;
    }

    void union_set(int x, int y) {
    	int id1 = find_set(x);
    	int id2 = find_set(y);
        if (id1 == id2) {	
        	return;
		}
		if(rank[id1] < rank[id2]) { 
			parent[id1] = id2;
			rank[id2] += rank[id1];
		} 
		else { 
			parent[id2] = id1;
			rank[id1] += rank[id2]; 
		}
		setSize--;
    }

    bool connected(int x, int y) {
        return (find_set(x) == find_set(y));
    }

    int size() {
        return setSize;
    }
};

void fredman_tarjan(UnionFind disjoint_set) {

	if(heapSize == 0) {	
		heapSize = num_edges/num_vertices;	
	}
	else {
		heapSize = pow(2,heapSize);
	}

	std::set<int> visited;

	for(auto vertex = vertices.begin(); vertex != vertices.end(); vertex++) {	

		if(visited.find(disjoint_set.find_set(*vertex)) == visited.end()) {
			
			Fib_heap priority_queue; 
			auto edge_iter = Adj_list[*vertex];
			unordered_map< int , Fib_handle > location;	

			for(auto v = edge_iter.begin(); v != edge_iter.end(); v++) {
				Edge edge;
				edge.v = (*v).first;
				edge.w = (*v).second;
				Fib_handle pos = priority_queue.emplace(edge);
				location[v->first] = pos;
			}

			while(!priority_queue.empty()) {

				auto edge = priority_queue.top();

				if((priority_queue.size() > heapSize) || 
				   (visited.find(edge.v) != visited.end()) || 
				   (disjoint_set.connected(edge.v, *vertex))) {
					break;
				}
				
				visited.insert(*vertex);
				visited.insert(edge.v);

				priority_queue.pop();
				mst_weight += edge.w;	
                ++mst_edge_count;
                
				disjoint_set.union_set(edge.v, *vertex);

				auto v_edge_iter = Adj_list[edge.v];
				for(auto v = v_edge_iter.begin(); v != v_edge_iter.end(); v++) {
					Edge edge;
					edge.v = (*v).first;
					edge.w = (*v).second;
					if( ((*v).first) == (*vertex)) {
						continue;
					}
					if(location.find((*v).first) == location.end())	{
						Fib_handle pos = priority_queue.emplace(edge);
						location[(*v).first] = pos;
					}
					else {
						if(visited.find((*v).first) == visited.end()) {
							Fib_handle pos = location[(*v).first];
							if( (*pos).w> (*v).second ) {
								priority_queue.update(pos, edge);
								priority_queue.update(pos);
							}
						}
					}
				}
			}
		}
	}

	if(disjoint_set.size() == 1) {
		return;
	}

	unordered_map<int, unordered_map<int, int> > new_adj_list;
	std::set<int> newvertices;

    /*
        The contraction is performed as follows:
    	1) Relabel the vertices based on the component they belong to.
    	2) The edges between vertices belonging to the same set are removed. 
    	3) Only the min-cost edge joining two different sets is kept.
    */
	for(auto it = vertices.begin(); it != vertices.end(); it++) {	
		int set_rep = disjoint_set.find_set(*it);
		if(newvertices.find(set_rep) == newvertices.end()) {
			newvertices.insert(set_rep);
		}
	}
	vertices = newvertices;

	unordered_map<int, int> mymap;
	for(auto it = newvertices.begin(); it != newvertices.end(); it++) {
		new_adj_list[*it] = mymap;
	}

	for(auto it = Adj_list.begin(); it != Adj_list.end(); it++) {

		for (auto vecit = ((*it).second).begin(); vecit != ((*it).second).end(); vecit++) {

			if( !disjoint_set.connected( (*it).first,(*vecit).first )) {

				if((new_adj_list[disjoint_set.find_set((*it).first)]).find((*vecit).first) 
                    == (new_adj_list[(*it).first]).end()) {
					new_adj_list[disjoint_set.find_set((*it).first)]
                                [disjoint_set.find_set((*vecit).first)] = (*vecit).second;
				}
				else if(new_adj_list[disjoint_set.find_set((*it).first)]
                        [disjoint_set.find_set( (*vecit).first)] > (*vecit).second) {
					new_adj_list[disjoint_set.find_set((*it).first)]
                                [disjoint_set.find_set( (*vecit).first)] = (*vecit).second;
				}
			}
		}
	}

	Adj_list = new_adj_list;
	
    fredman_tarjan(disjoint_set);
}

int main () {
	time_t start_time,end_time;
	ifstream f;
	f.open("Test/output_graph_2.txt");
	while(1) {
		string vertex;
		f>>vertex;
		if(vertex.compare("#") == 0) {
			break;
        }
		num_vertices++;
		vertices.insert(std::atoi(vertex.c_str()));
	}
	while(!f.eof()) {
	    int v1, v2, weight;
		num_edges++;
		f>>v1>>v2>>weight;
		Adj_list[v1][v2] = weight;
		Adj_list[v2][v1] = weight;
	}
	f.close();
	
	UnionFind disjoint_set(vertices);
    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::cout << "Finding MST .... \n" ;
    fredman_tarjan(disjoint_set);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";

	cout<<"The weight of the MST is "<<mst_weight<<endl;
	cout<<"The number of edges in MST are "<<mst_edge_count<<endl;

	return 1;
}
