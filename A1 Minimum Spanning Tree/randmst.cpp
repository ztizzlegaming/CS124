//CS124 Programming Assignemnt 1
//Qiang Fei and Jordan Turley
//
//Compile: g++ randmst.cpp -o randmst -O3
//Run:     ./randmst 0 numpoints numtrials dimensions
//
//To use graph from class, run: ./randmst 1 0 0 0

#include <iostream>
#include <vector>
#include <random>
#include <algorithm>
#include <math.h>

using namespace std;

const int DEBUG = 1;

//Holds the two endpoint vertices and the weight of an edge
struct Edge {
	int v1;
	int v2;
	double weight;
};

//Holds the two arrays (p, rank) needed for a disjoint set
struct DisjointSet {
	int *p;
	int *rank;
};

//Holds one element stored in a heap, which contains a vertex and a value/weight
struct HeapElem {
	int vertex;
	double value;
};

//Holds the elements of a heap
struct Heap {
	vector<struct HeapElem> heap;
};

//Holds the elements of a minimum spanning tree we need, tree weight and
//maximum edge weight
struct MinimumSpanningTree {
	double treeWeight;
	double maxEdgeWeight;
};

/**
 * Creates a new disjoint set struct with a given size
 */
struct DisjointSet makeDisjointSet(int size) {
	struct DisjointSet disjointSet;
	disjointSet.p = new int[size];
	disjointSet.rank = new int[size];
	return disjointSet;
}

/**
 * Frees up the memory used by the disjoint set. Should be called after you
 * are finished using the disjoint set.
 */
void freeDisjointSet(struct DisjointSet disjointSet) {
	delete[] disjointSet.p;
	delete[] disjointSet.rank;
}

/**
 * Makes a new set in a given disjoint set for a point x
 */
void ds_makeSet(int x, struct DisjointSet disjointSet) {
	disjointSet.p[x] = x;
	disjointSet.rank[x] = 0;
}

/**
 * Links x and y in a given disjoint set
 */
void ds_link(int x, int y, struct DisjointSet disjointSet) {
	if (disjointSet.rank[x] > disjointSet.rank[y]) {
		int t = x;
		x = y;
		y = t;
	}
	if (disjointSet.rank[x] == disjointSet.rank[y]) {
		disjointSet.rank[y] = disjointSet.rank[y] + 1;
	}
	disjointSet.p[x] = y;
}

/**
 * Finds the parent set for a given point x in a given disjoint set
 */
int ds_find(int x, struct DisjointSet disjointSet) {
	if (x != disjointSet.p[x]) {
		disjointSet.p[x] = ds_find(disjointSet.p[x], disjointSet);
	}
	return disjointSet.p[x];
}

/**
 * Takes the union of the sets containing x and y in a given disjoint set
 */
void ds_union(int x, int y, struct DisjointSet disjointSet) {
	ds_link(ds_find(x, disjointSet), ds_find(y, disjointSet), disjointSet);
}

/**
 * Initializes a new heap struct
 */
struct Heap makeHeap() {
	struct Heap heap;
	return heap;
}

/**
 * Initializes a new HeapElem with a given vertex and value
 */
struct HeapElem makeHeapElem(int vertex, double value) {
	struct HeapElem elem;
	elem.vertex = vertex;
	elem.value = value;
	return elem;
}

/**
 * Calculates the parent element index of a given heap element index
 */
int heap_parent(int i) {
	return i / 2;
}

/**
 * Calculates the left child index of a given heap element index
 */
int heap_left(int i) {
	return 2 * i;
}

/**
 * Calculates the right child index of a given heap element index
 */
int heap_right(int i) {
	return 2 * i + 1;
}

/**
 * Minifies a heap, guaranteeing the heap property that the parent element of
 * any element is smaller than the child element. Called after removing the
 * minimum element from the heap.
 */
void heap_minHeapify(int n, struct Heap *heap) {
	int l = heap_left(n);
	int r = heap_right(n);
	int smallest;
	if (l < heap->heap.size() && heap->heap[l].value != -1 && heap->heap[l].value < heap->heap[n].value) {
		smallest = l;
	} else {
		smallest = n;
	}

	if (r < heap->heap.size() && heap->heap[r].value != -1 && heap->heap[r].value < heap->heap[smallest].value) {
		smallest = r;
	}

	if (smallest != n) {
		struct HeapElem t = heap->heap[n];
		heap->heap[n] = heap->heap[smallest];
		heap->heap[smallest] = t;
		heap_minHeapify(smallest, heap);
	}
}

/**
 * Fetches the minimum element from the heap and removes it.
 */
struct HeapElem heap_removeMin(struct Heap *heap) {
	struct HeapElem min = heap->heap[0];
	heap->heap[0] = heap->heap[heap->heap.size() - 1];
	heap->heap.erase(heap->heap.begin() + heap->heap.size() - 1);
	heap_minHeapify(0, heap);
	return min;
}

/**
 * Inserts a new element into the heap and moves elements to ensure the heap
 * property still holds.
 */
void heap_insert(struct HeapElem heapElem, struct Heap *heap) {
	heap->heap.push_back(heapElem);
	int n = heap->heap.size() - 1;
	while (n != 0 && heap->heap[heap_parent(n)].value > heap->heap[n].value) {
		struct HeapElem t = heap->heap[heap_parent(n)];
		heap->heap[heap_parent(n)] = heap->heap[n];
		heap->heap[n] = t;
		n = heap_parent(n);
	}
}

/**
 * Creates a new edge struct with given endpoint vertices and weight
 */
struct Edge makeEdge(int v1, int v2, double weight) {
	struct Edge edge;
	edge.v1 = v1;
	edge.v2 = v2;
	edge.weight = weight;
	return edge;
}

//https://stackoverflow.com/questions/4892680/sorting-a-vector-of-structs
//Used for sorting edges. Sort edges according to weight in ascending order.
bool compareEdges(const Edge &e1, const Edge e2) {
	return e1.weight < e2.weight;
}

/**
 * Finds a minimum spanning tree in a given graph using Prim's algorithm.
 */
struct MinimumSpanningTree prim(int numPoints, vector< vector<struct Edge> > adjLists) {
	//Create the heap and insert the first vertex with a value of zero
	int s = 0;
	struct Heap heap = makeHeap();
	struct HeapElem e1 = makeHeapElem(s, 0);
	heap_insert(e1, &heap);

	double *dist = new double[numPoints];
	double *prev = new double[numPoints];
	bool *inTree = new bool[numPoints]; //simulates the set S of vertices in the tree
	for (int i1 = 0; i1 < numPoints; i1++) {
		dist[i1] = 1e10; //not exactly infinity, but this will be big enough
		prev[i1] = -1;
		inTree[i1] = false;
	}
	dist[s] = 0;
	inTree[s] = true;
	while (heap.heap.size() != 0) {
		struct HeapElem elem = heap_removeMin(&heap);
		int v = elem.vertex;
		inTree[v] = true;

		//Get v's adjacency list
		vector<struct Edge> adjList = adjLists[v];
		for (int i1 = 0; i1 < adjList.size(); i1++) {
			struct Edge edge = adjList[i1];

			//Get the other end of the edge
			int w = edge.v1;
			if (edge.v1 == v) {
				w = edge.v2;
			}

			//Only consider the edges to vertices not in the tree
			if (inTree[w]) {
				continue;
			}

			if (dist[w] > edge.weight) {
				dist[w] = edge.weight;
				prev[w] = v;
				struct HeapElem he = makeHeapElem(w, dist[w]);
				heap_insert(he, &heap);
			}
		}
	}

	//Free up the inTree and prev arrays
	delete[] inTree;
	delete[] prev; //TODO we don't need prev, so we can remove it later

	//Calculate the tree weight and the max edge weight
	double treeWeight = 0;
	double maxEdgeWeight = 0;
	for (int i1 = 0; i1 < numPoints; i1++) {
		treeWeight += dist[i1];
		if (dist[i1] > maxEdgeWeight) {
			maxEdgeWeight = dist[i1];
		}
	}

	//Free up the dist array
	delete[] dist;

	//Build the MST struct to return
	struct MinimumSpanningTree mst;
	mst.treeWeight = treeWeight;
	mst.maxEdgeWeight = maxEdgeWeight;

	return mst;
}

/**
 * Finds a minimum spanning tree in a given graph using Kruskal's algorithm.
 */
struct MinimumSpanningTree kruskal(int numPoints, vector< vector<struct Edge> > adjLists, vector<struct Edge> edges) {
	//Sort the edges by weight, least to greatest
	sort(edges.begin(), edges.end(), compareEdges);

	//Disjoint set data structure
	struct DisjointSet disjointSet = makeDisjointSet(numPoints);

	for (int i1 = 0; i1 < numPoints; i1++) {
		ds_makeSet(i1, disjointSet);
	}

	vector<struct Edge> tree;
	double treeWeight;

	for (int i1 = 0; i1 < edges.size(); i1++) {
		int u = edges[i1].v1;
		int v = edges[i1].v2;
		if (ds_find(u, disjointSet) != ds_find(v, disjointSet)) {
			tree.push_back(edges[i1]);
			treeWeight += edges[i1].weight;
			ds_union(u, v, disjointSet);
		}
	}

	//Clear up memory from disjoint set
	freeDisjointSet(disjointSet);

	//Build the MST struct to return
	struct MinimumSpanningTree mst;
	mst.treeWeight = treeWeight;
	mst.maxEdgeWeight = tree[tree.size() - 1].weight;

	return mst;
}

int main(int argc, char *argv[]) {
	//Make sure the correct arguments are given
	if (argc != 5) {
		cout << "Please run as the following: ./randmst 0 numpoints numtrials dimension" << endl;
		return 0;
	}

	//Get command line arguments
	int debug = stoi(argv[1]);
	int numPoints = stoi(argv[2]);
	int numTrials = stoi(argv[3]);
	int dimension = stoi(argv[4]);

	if (debug == DEBUG) {
		numPoints = 9;
	}

	if (debug == DEBUG) {
		//Create adjacency list and list of all edges with weights
		vector< vector<struct Edge> > adjLists(numPoints);
		vector<struct Edge> edges;

		//Testing code, test the graph from class
		adjLists[0].push_back(makeEdge(0, 1, 1));
		adjLists[0].push_back(makeEdge(0, 3, 3));
		adjLists[1].push_back(makeEdge(1, 0, 1));
		adjLists[1].push_back(makeEdge(1, 4, 5));
		adjLists[1].push_back(makeEdge(1, 2, 5));
		adjLists[2].push_back(makeEdge(2, 1, 5));
		adjLists[2].push_back(makeEdge(2, 5, 2));
		adjLists[3].push_back(makeEdge(3, 0, 3));
		adjLists[3].push_back(makeEdge(3, 4, 4));
		adjLists[3].push_back(makeEdge(3, 6, 2));
		adjLists[4].push_back(makeEdge(4, 1, 5));
		adjLists[4].push_back(makeEdge(4, 3, 4));
		adjLists[4].push_back(makeEdge(4, 5, 1));
		adjLists[4].push_back(makeEdge(4, 7, 5));
		adjLists[5].push_back(makeEdge(5, 2, 2));
		adjLists[5].push_back(makeEdge(5, 4, 1));
		adjLists[5].push_back(makeEdge(5, 8, 7));
		adjLists[6].push_back(makeEdge(6, 3, 2));
		adjLists[6].push_back(makeEdge(6, 7, 3));
		adjLists[7].push_back(makeEdge(7, 6, 3));
		adjLists[7].push_back(makeEdge(7, 4, 5));
		adjLists[7].push_back(makeEdge(7, 8, 6));
		adjLists[8].push_back(makeEdge(8, 7, 6));
		adjLists[8].push_back(makeEdge(8, 5, 7));

		struct Edge e1;
		e1.v1 = 0;
		e1.v2 = 1;
		e1.weight = 1;
		struct Edge e2;
		e2.v1 = 1;
		e2.v2 = 2;
		e2.weight = 5;
		struct Edge e3;
		e3.v1 = 0;
		e3.v2 = 3;
		e3.weight = 3;
		struct Edge e4;
		e4.v1 = 1;
		e4.v2 = 4;
		e4.weight = 5;
		struct Edge e5;
		e5.v1 = 2;
		e5.v2 = 5;
		e5.weight = 2;
		struct Edge e6;
		e6.v1 = 3;
		e6.v2 = 4;
		e6.weight = 4;
		struct Edge e7;
		e7.v1 = 4;
		e7.v2 = 5;
		e7.weight = 1;
		struct Edge e8;
		e8.v1 = 3;
		e8.v2 = 6;
		e8.weight = 2;
		struct Edge e9;
		e9.v1 = 4;
		e9.v2 = 7;
		e9.weight = 5;
		struct Edge e10;
		e10.v1 = 5;
		e10.v2 = 8;
		e10.weight = 7;
		struct Edge e11;
		e11.v1 = 6;
		e11.v2 = 7;
		e11.weight = 3;
		struct Edge e12;
		e12.v1 = 7;
		e12.v2 = 8;
		e12.weight = 6;

		edges.push_back(e1);
		edges.push_back(e2);
		edges.push_back(e3);
		edges.push_back(e4);
		edges.push_back(e5);
		edges.push_back(e6);
		edges.push_back(e7);
		edges.push_back(e8);
		edges.push_back(e9);
		edges.push_back(e10);
		edges.push_back(e11);
		edges.push_back(e12);

		//Run Prim and Kruskal
		struct MinimumSpanningTree mstPrim = prim(numPoints, adjLists);
		cout << "Tree weight from Prim: " << mstPrim.treeWeight << endl;
		cout << "Max edge weight from Prim: " << mstPrim.maxEdgeWeight << endl;

		struct MinimumSpanningTree mstKruskal = kruskal(numPoints, adjLists, edges);
		cout << "Tree weight from Kruskal: " << mstKruskal.treeWeight << endl;
		cout << "Max edge weight from Kruskal: " << mstKruskal.maxEdgeWeight << endl;
	} else {
		//Create random number generator and seed with current time
		mt19937 generator(time(0));
		uniform_real_distribution<double> distribution(0.0, 1.0);

		double averageTreeWeight = 0;
		double overallMaxEdgeWeight = 0;
		double averageMaxEdgeWeight = 0;
		for (int i1 = 0; i1 < numTrials; i1++) {
			//Create adjacency list and list of all edges with weights
			vector< vector<struct Edge> > adjLists(numPoints); //[numPoints];
			//vector<struct Edge> edges; //UNCOMMENT FOR KRUSKAL, commented out to save memory

			if (dimension == 0) {
				double edgeThreshold = (double) (32.0 / numPoints);
				for (int i1 = 0; i1 < numPoints; i1++) {
					for (int i2 = 0; i2 < i1; i2++) { //Loop over each pair of edges
						double weight = distribution(generator);

						//Generate the edge
						struct Edge edge = makeEdge(i1, i2, weight);

						if (weight < edgeThreshold) {
							//edges.push_back(edge); //UNCOMMENT FOR KRUSKAL
						
							//Insert each vertex into the other's adj. list
							adjLists[i1].push_back(edge);
							adjLists[i2].push_back(edge);
						}
					}
				}
			} else if (dimension == 2 || dimension == 3 || dimension == 4) {
				double edgeThreshold;
				if (dimension == 2) {
					edgeThreshold = (double) (4.0 / pow((double) numPoints, 1.0 / 2.0));
				} else if (dimension == 3) {
					edgeThreshold = (double) (2.5 / pow((double) numPoints, 1.0 / 3.0));
				} else if (dimension == 4) {
					edgeThreshold = (double) (2.1 / pow((double) numPoints, 1.0 / 4.0));
				}

				double **points = new double*[numPoints];
				for (int i1 = 0; i1 < numPoints; i1++) {
					points[i1] = new double[dimension];
					for (int i2 = 0; i2 < dimension; i2++) {
						points[i1][i2] = distribution(generator);
					}
				}

				//Loop over each pair of points and calculate the Euclidean distance
				for (int i1 = 0; i1 < numPoints; i1++) {
					for (int i2 = 0; i2 < i1; i2++) {
						//Calculat the edge weight as the euclidean distance btwn points
						double weight = 0;
						for (int i3 = 0; i3 < dimension; i3++) {
							weight += (points[i1][i3] - points[i2][i3]) * (points[i1][i3] - points[i2][i3]);
						}
						weight = sqrt(weight);

						if (weight < edgeThreshold) {
							//Generate the edge
							struct Edge edge = makeEdge(i1, i2, weight);
							//edges.push_back(edge); //UNCOMMENT FOR KRUSKAL

							//Insert each edge into the other's adj. list
							adjLists[i1].push_back(edge);
							adjLists[i2].push_back(edge);
						}
					}
				}

				//Free up the points array, we don't need it anymore
				for (int i1 = 0; i1 < numPoints; i1++) {
					delete[] points[i1];
				}
				delete[] points;
			}

			struct MinimumSpanningTree mst = prim(numPoints, adjLists);
			//struct MinimumSpanningTree mst = kruskal(numPoints, adjLists, edges); //UNCOMMENT FOR KRUSKAL
			averageTreeWeight += mst.treeWeight;
			averageMaxEdgeWeight += mst.maxEdgeWeight;
			if (mst.maxEdgeWeight > overallMaxEdgeWeight) {
				overallMaxEdgeWeight = mst.maxEdgeWeight;
			}
		}

		//Calculate averages and display results
		averageTreeWeight /= numTrials;
		averageMaxEdgeWeight /= numTrials;

		cout << "Average tree weight: " << averageTreeWeight << endl;
		cout << "Overall max edge weight: " << overallMaxEdgeWeight << endl;
		cout << "Average max edge weight: " << averageMaxEdgeWeight << endl;
	}

	return 0;
}
