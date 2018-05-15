#include "util.h"
#include <chrono>

using namespace std;

vector<string> vertices;	
unordered_map<string, unordered_map<string,int> > adj_lis;
unordered_map<string, string> parent;	
unordered_map< string , bin_handle > pointer;
long long mst_weight = 0;


void prim() {
    set<string> visited;
    bin_heap pq; 
    for(auto vit = vertices.begin(); vit != vertices.end(); ++vit) {
        parent[*vit] = *vit;
        pair<string,int> somepair = make_pair(*vit , INF);
	bin_handle pos = pq.emplace(somepair);
	pointer[*vit] = pos;
    }

    pair<string,int> mypair = make_pair(vertices[0] , 0);
    pq.update(pointer[vertices[0]], mypair);

    while(!pq.empty()) {

        auto it = (pq.top());
	pq.pop();
	mst_weight += it.second;	
	visited.insert(it.first);	
	auto edgeit = adj_lis[it.first];
	bin_handle loophandler;
	for(auto vecit = edgeit.begin(); vecit != edgeit.end(); vecit++) {
            if(visited.find((*vecit).first) == visited.end()) {
	        loophandler = pointer[(*vecit).first];
		if( (*loophandler).second > (*vecit).second ) {
		    parent[(*vecit).first] = it.first;
	            mypair = make_pair( (*vecit).first, (*vecit).second);
	            pq.update(loophandler, mypair);
		    pq.update(loophandler);
		}
	    }
       }	
    }
}


void parse_input_file() {
    string filename("Test/output_graph_2.txt");
    ifstream f;
    f.open(filename);
    string name;

    while(1) {
        f>>name;
        if(name.compare("#") == 0) {
            break;
        }
        vertices.push_back(name);
    }
    string v1,v2;
    int weight;
    while( !f.eof() ) {
        f>>v1>>v2>>weight;
        adj_lis[v1][v2] = weight;
        adj_lis[v2][v1] = weight;
    }
    f.close();
}


void print_mst_edges() {
    std::cout<<"The Edges of the MST are:"<<std::endl;
    for(auto hashit = parent.begin(); hashit != parent.end(); ++hashit) {
        std::cout<<hashit->first<<" "<<hashit->second<<std::endl;
    }
}


int main(int argc, char** argv) { 
    parse_input_file();
    using namespace std::chrono;    
    high_resolution_clock::time_point start_time = high_resolution_clock::now();
    prim();
    high_resolution_clock::time_point end_time = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(end_time - start_time);
    std::cout << "Exectution Time: " << time_span.count() << " seconds.";
    std::cout << std::endl;
    //print_mst_edges();
    return 0;
}

