CXX=icc
CXXFLAGS=-std=c++11 -I boost_1_65_1/

all:
	$(CXX) $(CXXFLAGS) -o graph_generator graph_generator.cpp
	$(CXX) $(CXXFLAGS) -o dir_graph_generator dir_graph_generator.cpp
	$(CXX) $(CXXFLAGS) -o fredman-tarjan fredman-tarjan.cpp
	$(CXX) $(CXXFLAGS) -o edmonds-algo-iter edmonds-algo-iter.cpp 
	$(CXX) $(CXXFLAGS) -o edmonds-algo-par edmonds-algo-par.cpp

clean:
	rm -rf graph_generator
	rm -rf fredman_tarjan
	rm -rf edmonds-algo-iter 
	rm -rf edmonds-algo-par
