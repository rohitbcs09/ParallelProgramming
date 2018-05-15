#include <bits/stdc++.h>
#include <unordered_map>
#include <cilk/cilk.h>
#include <cilk/cilk_api.h>
#include <chrono>

using namespace std;

const int INF = 1000000;

struct Edge {
    int u, v, w;

    Edge(int u, int v, int w) : u(u), v(v), w(w) {}

    bool operator <(const Edge& x) const {
        return w < x.w;
    }
};

void Par_Simulate_Priority_CW_BS(int n, std::vector<Edge> &E, std::vector<int> &R) {
    std::vector<int> B(n+1, 0);
    std::vector<int> l(n+1, 0);
    std::vector<int> h(n+1, 0);

    std::vector<int> lo(n+1, 0);
    std::vector<int> hi(n+1, 0);
    std::vector<int> md(n+1, 0);

    int m = E.size();

    #pragma cilk grainsize = 2048
    cilk_for(int u = 1; u<=n; ++u) {
        l[u] = 0;
        h[u] = m-1;
    }

    for(int k = 0; k < 1 + log2(m); ++k){
	#pragma cilk grainsize = 2048
        cilk_for(int u = 1; u<=n; ++u) {
            B[u] = 0;
            lo[u] = l[u];
            hi[u] = h[u];
        }
	#pragma cilk grainsize = 2048
        cilk_for(int i = 0; i< m; ++i) {
            int u = E[i].u;
            md[u] = (lo[u] + hi[u])/2;
            if( i>=lo[u] && i<=md[u]) {
                B[u] = 1;
            }
        }
	#pragma cilk grainsize = 2048
        cilk_for(int i = 0; i< m; ++i) {
            int u = E[i].u;
            md[u] = (lo[u] + hi[u])/2;
            if( B[u] == 1 && i>=lo[u] && i<=md[u]) {
                h[u] = md[u];
            }
            else if(B[u] == 0 && i>=md[u] &&
                    i<=hi[u] ) {
                l[u] = md[u] + 1;
            }
        }
	#pragma cilk grainsize = 2048
        cilk_for(int i = 0; i<m; ++i) {
            int u = E[i].u;
            if(i==l[u]) {
                R[u] = i;
            }
        }
    }
}

int edmonds(vector<Edge>& edgeList, int V, int R) {
    
    // Parallel QuickSort to sort the weights of graph
    std::sort(edgeList.begin(), edgeList.end());

    std::vector<int> Rank(V+1, 0);

    vector<Edge> minInEdge(V+1, Edge(0, 0, INF));

    //for (Edge e : edgeList) {
    /*//#pragma cilk grainsize = 40
    for (int i = 0; i< edgeList.size(); ++i) {
        Edge e = edgeList[i]; 
        minInEdge[e.v] = min(minInEdge[e.v], e);
    }*/

    minInEdge[R] = Edge(0, R, 0);
    
    Par_Simulate_Priority_CW_BS(V, edgeList, Rank);

    for(int i = 0; i < edgeList.size(); ++i) {
        Edge e(edgeList[i].u, edgeList[i].v, edgeList[i].w);
        if(Rank[e.u] == i) {
           minInEdge[e.u] = e; //min(minInEdge[e.u], e); 
        }
    }

    // assign vertices to their cyclic group
    vector<int> group(V+1, 0);
    vector<bool> visited(V+1, false), isCycleGroup(V+1, false);
    int cnt = 1;
    for (int i = 1; i < V+1; i++) {
        if (visited[i])
            continue;

        int node = i; 
        vector<int> path;
        while (node != 0 && !visited[node]) {
            visited[node] = true;
            path.push_back(node);
            node = minInEdge[node].u;
        }

        bool isCycle = false;
        for (int v : path) {
            group[v] = cnt;
            if (v == node)
                isCycleGroup[cnt] = isCycle = true;
            if (!isCycle)
                cnt++;
        }
        if (isCycle)
            cnt++;
    }

    // when there are no cycles
    if (cnt == V+1) {
        int result = 0;
        for(Edge e : minInEdge) {
            if(e.w != INF) {
                //std::cout << e.u << " " << e.v << " " << e.w << "\n";
                result += e.w;
            }
        }
        return result;
    }

    int result = 0;
    for (Edge e : minInEdge)
        if (isCycleGroup[group[e.v]]) {
            if(e.w != INF) {
                //std::cout << e.u << " " << e.v << " " << e.w << "\n";
                result += e.w;
            }
        }

    // form new graph with groups
    vector<Edge> n_edgeList;
    for (Edge e : edgeList) {
        int u = group[e.u], v = group[e.v], w = e.w;
        if (u == v)
            continue;
        else if (isCycleGroup[v])
            n_edgeList.push_back(Edge(u, v, w - minInEdge[e.v].w));
        else
            n_edgeList.push_back(Edge(u, v, w));
    }

    return result + edmonds(n_edgeList, cnt, group[R]);
}

int main(int argc, char** argv) {

    __cilkrts_set_param("nworkers", "68");

    int V, E, R=1; 
    vector<Edge> edgeList;
    std::unordered_map<int, bool> vertex;

    std::string line;
    std::ifstream infile("Test/Directed_Weighted_Graph.in");
    std::getline(infile, line);
    std::istringstream iss(line);
    while(std::getline(infile, line)) {
        int u, v, w;
        std::istringstream iss(line);
        iss >> u >> v >> w; 
        edgeList.push_back(Edge(u, v, w));
        vertex.insert(std::pair<int, bool>(u, true));
        vertex.insert(std::pair<int, bool>(v, true));
    }

    E = edgeList.size();
    V = vertex.size();

    using namespace std::chrono;
    high_resolution_clock::time_point t1 = high_resolution_clock::now();

    std::cout << "Finding MST .... \n" ;
    int result = edmonds(edgeList, V, R);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";
    //printf("%d\n", result);
    return 0;
}
