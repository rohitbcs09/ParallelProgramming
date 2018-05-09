#include <bits/stdc++.h>                                                                                                                 
#include <boost/heap/fibonacci_heap.hpp>                                        
#include <boost/heap/d_ary_heap.hpp>                                            
                                                                                
#define INF numeric_limits<int>::max()                                          
                                                                                
using namespace std;                                                            

typedef struct edge {
    int v;
    int w;
} Edge;
                                                                                
struct Comparator {
    bool operator() (const Edge& a, const Edge& b) const {
    	return a.w > b.w;
    }
};
                                                                       
typedef boost::heap::d_ary_heap< Edge, boost::heap::arity<2>, boost::heap::compare< Comparator >,
                                 boost::heap::mutable_<true> > Bin_heap;

typedef boost::heap::d_ary_heap< Edge, boost::heap::arity<2>, boost::heap::compare< Comparator >,
                                 boost::heap::mutable_<true> >::handle_type Bin_handle;
    
typedef boost::heap::fibonacci_heap< Edge, boost::heap::compare< Comparator >,
                                     boost::heap::mutable_<true> > Fib_heap;

typedef boost::heap::fibonacci_heap< Edge, boost::heap::compare< Comparator >,
                                     boost::heap::mutable_<true> >::handle_type Fib_handle;
