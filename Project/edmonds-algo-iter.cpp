#include <bits/stdc++.h>
#include <unordered_map>

using namespace std;

const int INF = 1000000;

struct Edge {
    int u, v, w;

    Edge(int u, int v, int w) : u(u), v(v), w(w) {}

    // order edges by weight
    bool operator <(const Edge& x) const {
        return w < x.w;
    }
};

int edmonds(vector<Edge>& edgeList, int V, int R) {
    
    // determine min cost of edge entering each vertex
    vector<Edge> minInEdge(V+1, Edge(-1, -1, INF));

    for (Edge e : edgeList) {
        minInEdge[e.v] = min(minInEdge[e.v], e);
    } 
    
    minInEdge[R] = Edge(-1, R, 0);

    // assign vertices to their cyclic group
    vector<int> group(V+1, 0);
    vector<bool> visited(V+1, false), isCycleGroup(V+1, false);
    int cnt = 1;
    for (int i = 1; i < V+1; i++) {
        if (visited[i])
            continue;

        int node = i; vector<int> path;
        while (node != -1 && !visited[node]) {
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
        for (Edge e : minInEdge) {
            if(e.w != INF) {
                std::cout << e.u << " " << e.v << " " << e.w << "\n";
                result += e.w;
            }
        }
        return result;
    }

    int result = 0;
    for (Edge e : minInEdge) {
        if (isCycleGroup[group[e.v]]) {
            std::cout << e.u << " " << e.v << " " << e.w << "\n";
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
    std::sort(edgeList.begin(), edgeList.end());
    int result = edmonds(edgeList, V, R);

    high_resolution_clock::time_point t2 = high_resolution_clock::now();
    duration<double> time_span = duration_cast<duration<double>>(t2 - t1);
    std::cout << "Running Time: " << time_span.count() << " seconds.\n";
    //printf("%d\n", result);
    return 0;
}
