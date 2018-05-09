#include <bits/stdc++.h>
using namespace std;

class GraphConstructor {
    
    vector<vector<int>> adjacency_list;

    long long num_verts = 0;
    long long num_edges = 0;
    long long pres_edges = 0;
    int gen_file_no = 0;

    public :
	GraphConstructor(long long, long long, int);
	void build_spanning_tree();
	void add_random_edges(); 
 	void write_graph_to_file();

};

GraphConstructor :: GraphConstructor(long long a, long long b, int c) {
    num_verts = a;
    num_edges = b;
    gen_file_no = c;
}


void GraphConstructor :: build_spanning_tree() {
    srand (time(NULL));
    adjacency_list.resize(num_verts);

    for(int i=1; i<num_verts; i++) {
        int random_vertex = rand()%i;
        adjacency_list[random_vertex].push_back(i);
        pres_edges++;
    }
}


void GraphConstructor :: add_random_edges() {
    long long max_edges = num_verts*num_verts;
    long long counter = 0;

    // Two vertices are randomly generated, if they are not connected, they are connected.
    while((counter < max_edges) && (pres_edges < num_edges))
    {
        int v1 = rand()%num_verts;
        int v2 = rand()%num_verts;

        if(v1 > v2) {
            swap(v1,v2);
        }

        if( (!(v1==v2)) && (find(adjacency_list[v1].begin(), adjacency_list[v1].end(), v2) == adjacency_list[v1].end()) )
        {
           // An Edge connecting v1 and v2 does not exist.
           adjacency_list[v1].push_back(v2);
           //adj_lis[v2].push_back(v1);
           pres_edges++;
        }
        counter++;
    }
}

void GraphConstructor :: write_graph_to_file() {

    int max_weight = 50;

    ofstream f;
    string s("Test/output_graph_");
    s += to_string(gen_file_no);
    s += ".txt";
    f.open(s);

    bool first = true;

    cout<<"Writing generated file to - " << s <<endl;

    // Writing the names of the vertices in the file.
    for(int i = 0; i < num_verts; ++i)
    {
        if(first) {
                f<<i;
                first = false;
        }
        else
                f<<endl<<i;
    }

    f<<endl<<"#";

    for(int i = 0; i<num_verts; i++)
    {
        for(int j=0; j< (adjacency_list[i]).size(); j++)
        {
            int weight = (rand() % max_weight) + 1;
		    if (first) {
		    	f<<i<<" "<<adjacency_list[i][j]<<" "<<weight;
		    	first = false;
		    } 
            else {
                f<<endl<<i<<" "<<adjacency_list[i][j]<<" "<<weight;
		    }
        }
    }

    f.close();

}

int main(int argc, char** argv){
 
    long long num_verts = stoll(argv[1]);
    long long num_edges = stoll(argv[2]);
    int gen_file_no = atoi(argv[3]);
     
    GraphConstructor gfc(num_verts, num_edges, gen_file_no);
    gfc.build_spanning_tree();
    gfc.add_random_edges();
    gfc.write_graph_to_file();
    return 0;
}
