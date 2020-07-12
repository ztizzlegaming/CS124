# Programming Assignment 1 - Minimum Spanning Trees
## Harvard University, CS 124
### Qiang Fei, Jordan Turley

We developed code to generate random graphs and use the Prim and Kruskal algorithms to find the average minimum spanning tree weight.

### Compile

Use make:

    make
    
Or compile manually:

    g++ -std=c++0x randmst.cpp -o randmst -O3

### Run

Run for a given number of points in the graph, for a given number of trials, for a given dimension:

    ./randmst 0 [numpoints] [numtrials] [dimension]

For dimension = 0, the edge weights are Uniform(0, 1). For dimension = 2, 3, or 4, the points are generated in n-dimensional space and the edge weights are calculated as the Euclidean distance between the points.
