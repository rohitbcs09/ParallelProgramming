#include<bits/stdc++.h>
using namespace std;
 
#define RUN 5
#define MAX_VERTICES 50 
#define MAX_EDGES 100 
#define MAXWEIGHT 100
 
int main() {
    set<pair<int, int> > container;
    set<pair<int, int> >::iterator it;
 
    freopen("Test/Directed_Weighted_Graph.in", "w", stdout);
 
    srand(time(NULL));
 
    int NUM;    
    int NUMEDGE;
 
    for (int i=1; i<=RUN; i++) {

        NUM = 1 + rand() % MAX_VERTICES;
        NUMEDGE = 1 + rand() % MAX_EDGES;
 
        while (NUMEDGE > NUM*(NUM-1)/2) {
            NUMEDGE = 1 + rand() % MAX_EDGES;
        }
 
        for (int j=1; j<=NUMEDGE; j++) {

            int a = 1 + rand() % NUM;
            int b = 1 + rand() % NUM;
            pair<int, int> p = make_pair(a, b);
 
            // Search for a random "new" edge every time
            // Note - In a tree the edge (a, b) is same
            // as the edge (b, a)
            while (container.find(p) != container.end()) {
                a = 1 + rand() % NUM;
                b = 1 + rand() % NUM;
                p = make_pair(a, b);
            }
            container.insert(p);
        }
 
        for (it=container.begin(); it!=container.end(); ++it)
        {
            int wt = 1 + rand() % MAXWEIGHT;
            printf("%d %d %d\n", it->first, it->second, wt);
        }
 
        container.clear();
        printf("\n");
 
    }
 
    fclose(stdout);
    return(0);
}
