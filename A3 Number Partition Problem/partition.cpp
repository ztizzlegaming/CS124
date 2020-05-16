#include <iostream>
#include <vector>
#include <random>
#include <chrono>

#define PRODUCTION 0
#define SIMULATION 1
#define DEBUG 2

#define KARMARKAR_KARP 0
#define REPEATED_RANDOM 1
#define HILL_CLIMBING 2
#define SIMULATED_ANNEALING 3
#define PREPARTITIONED_REPEATED_RANDOM 11
#define PREPARTITIONED_HILL_CLIMBING 12
#define PREPARTITIONED_SIMULATED_ANNEALING 13

using namespace std;

//Variables for random number generation
mt19937 generator(time(0));
uniform_real_distribution<double> distribution(0.0, 1.0);

//Holds one element stored in a heap, which is just a value for this problem
struct HeapElem {
	long value;
};

//Holds the elements of a heap
struct Heap {
	vector<struct HeapElem> heap;
};

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
struct HeapElem makeHeapElem(long value) {
	struct HeapElem elem;
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
	if (l < heap->heap.size() && heap->heap[l].value != -1 && heap->heap[l].value > heap->heap[n].value) {
		smallest = l;
	} else {
		smallest = n;
	}

	if (r < heap->heap.size() && heap->heap[r].value != -1 && heap->heap[r].value > heap->heap[smallest].value) {
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
struct HeapElem heap_removeMax(struct Heap *heap) {
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
	while (n != 0 && heap->heap[heap_parent(n)].value < heap->heap[n].value) {
		struct HeapElem t = heap->heap[heap_parent(n)];
		heap->heap[heap_parent(n)] = heap->heap[n];
		heap->heap[n] = t;
		n = heap_parent(n);
	}
}

/**
 * Calculate the residue of a list of numbers using the Karmarkar-Karp algorithm
 */
long karmarkar_karp(long* nums, int n) {
	//Add everything to the heap
	struct Heap heap = makeHeap();
	for (int i1 = 0; i1 < n; i1++) {
		struct HeapElem elem = makeHeapElem(nums[i1]);
		heap_insert(elem, &heap);
	}

	//Run the Karmarkar-Karp algorithm
	for (int i1 = 0; i1 < n - 1; i1++) {
		//Take the two largest elements in the heap
		struct HeapElem e1 = heap_removeMax(&heap);
		struct HeapElem e2 = heap_removeMax(&heap);

		//Calculate the difference and put it back on the heap
		long diff = e1.value - e2.value;
		struct HeapElem newElem = makeHeapElem(diff);
		heap_insert(newElem, &heap);
	}

	//Get the final value, this is the residue
	struct HeapElem finalElem = heap_removeMax(&heap);
	return finalElem.value;
}

/**
 * Transfer every element of b into a. Implemented for arrays of longs.
 * Sets a[i] = b[i] for i = 1, ..., n
 */
void transfer(long* a, long* b, int n) {
	for (int i1 = 0; i1 < n; i1++) {
		a[i1] = b[i1];
	}
}

/**
 * Generates a random uniform from zero to one
 */
double random_uniform() {
	//return (double) rand() / RAND_MAX;
	return distribution(generator);
}

/**
 * Generate a random S array, i.e. a random array of -1 and 1.
 */
long* generate_random_s(int n) {
	long* s = new long[n];
	for (int i1 = 0; i1 < n; i1++) {
		if (random_uniform() < 0.5) {
			s[i1] = -1;
		} else {
			s[i1] = 1;
		}
	}
	return s;
}

/**
 * Generate a random neighbor of a given S array.
 */
long* random_neighbor_s(long* s, int n) {
	long* new_s = new long[n];
	transfer(new_s, s, n);

	long i = random_uniform() * (double) n;
	long j = random_uniform() * (double) n;
	while (i == j) {
		j = random_uniform() * (double) n;
	}

	new_s[i] = -new_s[i];
	if (random_uniform() < 0.5) {
		new_s[j] = -new_s[j];
	}

	return new_s;
}

/**
 * Calculate the residue for the given list of numbers and a given S array
 */
long s_residue(long* nums, long* s, int n) {
	long residue = 0;
	for (int i1 = 0; i1 < n; i1++) {
		residue += s[i1] * nums[i1];
	}
	residue = abs(residue);
	return residue;
}

/**
 * Generate a solution S using the repeated random method.
 * Return the residue of the solution found
 */
long repeated_random(long* nums, int n, int max_iter) {
	//Declare variables
	long *s, *s_prime;
	int i1;
	long s_resid, s_prime_resid, residue;
	
	//Start with a random solution S
	s = generate_random_s(n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate a new random solution S'
		s_prime = generate_random_s(n);

		//Calculate residue(S) and residue(S')
		s_resid = s_residue(nums, s, n);
		s_prime_resid = s_residue(nums, s_prime, n);

		//If the 
		if (s_prime_resid < s_resid) {
			transfer(s, s_prime, n);
		}

		delete[] s_prime;
	}

	residue = s_residue(nums, s, n);

	delete[] s;

	return residue;
}

/**
 * Generate a solution S using the hill climbing method.
 * Return the residue of the solution found
 */
long hill_climbing(long* nums, int n, int max_iter) {
	//Declare variables
	long *s, *s_prime;
	int i1;
	long s_resid, s_prime_resid, residue;

	//Start with a random solution S
	s = generate_random_s(n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate a random neighbor S'
		s_prime = random_neighbor_s(s, n);

		//Calculate residue(S) and residue(S')
		s_resid = s_residue(nums, s, n);
		s_prime_resid = s_residue(nums, s_prime, n);

		//If residue(S') < residue(S), S = S'
		if (s_prime_resid < s_resid) {
			transfer(s, s_prime, n);
		}

		delete[] s_prime;

		//End loop early if we have found an optimal solution
		if (s_prime_resid == 0) {
			break;
		}
	}

	//Calculate the residue
	residue = s_residue(nums, s, n);

	delete[] s;

	return residue;
}

/**
 * Suggested T function for use with simulated annealing
 */
double T(long iter) {
	return 10000000000.0 * pow(0.8, floor(iter / 300));
}

/**
 * Generate a solution S using the simulated annealing method.
 * Return the residue of the solution found
 */
long simulated_annealing(long* nums, int n, int max_iter) {
	//Declare variables
	long *s, *s_prime, *s_dprime;
	int i1;
	long s_resid, s_prime_resid, s_dprime_resid, residue;
	double p;

	//Start with a random solution S
	s = generate_random_s(n);

	//Set S'' = S
	s_dprime = new long[n];
	transfer(s_dprime, s, n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate S', a random neighbor of S
		s_prime = random_neighbor_s(s, n);
		
		//Calculate residue(S) and residue(S')
		s_resid = s_residue(nums, s, n);
		s_prime_resid = s_residue(nums, s_prime, n);

		p = exp(-(s_prime_resid - s_resid) / T(i1));

		//If residue(S') < residue(S) then S = S'
		if (s_prime_resid < s_resid) {
			transfer(s, s_prime, n);
		} else if (random_uniform() < p) { //S = S' with prob p
			transfer(s, s_prime, n);
		}

		delete[] s_prime;

		s_dprime_resid = s_residue(nums, s_dprime, n);

		//If residue(S) < residue(S'') then S'' = S
		if (s_resid < s_dprime_resid) {
			transfer(s_dprime, s, n);
		}

		//End the loop early if we have found an optimal solution
		if (s_dprime_resid == 0) {
			break;
		}
	}

	delete[] s;

	//Calculate the final residue
	residue = s_residue(nums, s_dprime, n);

	delete[] s_dprime;

	return residue;
}

/**
 * Generate a random partutation, i.e. an array P = {p_1, ... p_n} where each
 * p_i is in {1, 2, ..., n}.
 */
long* random_part(int n) {
	long* part = new long[n];
	for (int i1 = 0; i1 < n; i1++) {
		part[i1] = random_uniform() * n + 1;
	}
	return part;
}

/**
 * Generate a random neighbor of a given partutation
 */
long* random_neighbor_part(long* part, int n) {
	long* new_part = new long[n];
	transfer(new_part, part, n);

	long i = random_uniform() * n + 1;
	long j = random_uniform() * n + 1;

	while (part[i - 1] == j) {
		j = random_uniform() * n + 1;
	}

	new_part[i - 1] = j;

	return new_part;
}

/**
 * Converts an array of numbers and a given P array back to the
 * standard representation.
 */
long* convert_part_to_standard(long* nums, long* part, int n) {
	long* a_prime = new long[n];
	for (int i1 = 0; i1 < n; i1++) {
		a_prime[i1] = 0;
	}
	for (int j = 0; j < n; j++) {
		a_prime[part[j] - 1] = a_prime[part[j] - 1] + nums[j];
	}
	return a_prime;
}

/**
 * Generate a solution P using the partutation representation and the
 * repeated random method.
 * Return the residue of the solution found.
 */
long repeated_random_part(long* nums, int n, int max_iter) {
	//Declare variables
	long *p, *p_prime, *p_std, *p_prime_std;
	int i1, p_resid, p_prime_resid, residue; 

	//Start with a random solution
	p = random_part(n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate a new random solution
		p_prime = random_part(n);

		//Convert P and P' back to standard representation
		p_std = convert_part_to_standard(nums, p, n);
		p_prime_std = convert_part_to_standard(nums, p_prime, n);

		//Calculate residue(P) and residue(P') using Karmarkar-Karp
		p_resid = karmarkar_karp(p_std, n);
		p_prime_resid = karmarkar_karp(p_prime_std, n);

		//Free up memory
		delete[] p_std;
		delete[] p_prime_std;

		//If the new residue is lower, replace the solution with the new
		if (p_prime_resid < p_resid) {
			transfer(p, p_prime, n);
		}

		//Free up new partutation
		delete[] p_prime;

		//If we reach zero, we have an optimal solution so end the loop
		if (p_prime_resid == 0) {
			break;
		}
	}

	//Convert P back to standard representation
	p_std = convert_part_to_standard(nums, p, n);

	//Calculate the final residue, free memory, and return
	residue = karmarkar_karp(p_std, n);

	delete[] p_std;
	delete[] p;

	return residue;
}

/**
 * Generate a solution P using the partutation representation and the
 * hill climbing method.
 * Return the residue of the solution found.
 */
long hill_climbing_part(long* nums, int n, int max_iter) {
	//Declare variables
	long *p, *p_prime, *p_std, *p_prime_std;
	int i1, p_resid, p_prime_resid, residue;

	//Start with a random solution
	p = random_part(n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate a random neighbor
		p_prime = random_neighbor_part(p, n);

		//Convert P and P' back to standard representation
		p_std = convert_part_to_standard(nums, p, n);
		p_prime_std = convert_part_to_standard(nums, p_prime, n);

		//Calculate residue(P) and residue(P') using Karmarkar-Karp
		p_resid = karmarkar_karp(p_std, n);
		p_prime_resid = karmarkar_karp(p_prime_std, n);

		delete[] p_std;
		delete[] p_prime_std;

		//If the new residue is lower, replace the solution with the new
		if (p_prime_resid < p_resid) {
			transfer(p, p_prime, n);
		}

		delete[] p_prime;

		//If we reach zero, we have an optimal solution so end the loop
		if (p_prime_resid == 0) {
			break;
		}
	}

	//Convert P to standard representation
	p_std = convert_part_to_standard(nums, p, n);

	//Calculate the final residue, free memory, and return
	residue = karmarkar_karp(p_std, n);

	delete[] p_std;
	delete[] p;

	return residue;
}

/**
 * Generate a solution P using the partutation representation and the
 * simulated annealing method.
 * Return the residue of the solution found.
 */
long simulated_annealing_part(long* nums, int n, int max_iter) {
	long *p, *p_prime, *p_dprime, *p_std, *p_prime_std, *p_dprime_std;
	int i1, p_resid, p_prime_resid, p_dprime_resid, residue;
	double pr;

	//Start with a random solution
	p = random_part(n);

	//Set S'' = S
	p_dprime = new long[n];
	transfer(p_dprime, p, n);

	//Loop for maximum number of iterations
	for (i1 = 0; i1 < max_iter; i1++) {
		//Generate S', a random neighbor of S
		p_prime = random_neighbor_part(p, n);

		//Convert P and P' to standard representation
		p_std = convert_part_to_standard(nums, p, n);
		p_prime_std = convert_part_to_standard(nums, p_prime, n);

		//Calculate residue(P) and residue(P')
		p_resid = karmarkar_karp(p_std, n);
		p_prime_resid = karmarkar_karp(p_prime_std, n);

		//Free up memory
		delete[] p_std;
		delete[] p_prime_std;

		pr = exp(-(p_prime_resid - p_resid) / T(i1));

		//If residue(P') < residue(P) then P = P'
		if (p_prime_resid < p_resid) {
			transfer(p, p_prime, n);
		} else if (random_uniform() < pr) {
			transfer(p, p_prime, n);
		}

		delete[] p_prime;

		//Convert P'' to standard representation
		p_dprime_std = convert_part_to_standard(nums, p_dprime, n);
		p_dprime_resid = karmarkar_karp(p_dprime_std, n);
		delete[] p_dprime_std;

		//If residue(P) < residue(P'') then P'' = P
		if (p_resid < p_dprime_resid) {
			transfer(p_dprime, p, n);
		}

		if (p_dprime_resid == 0) {
			break;
		}
	}

	delete[] p;

	//Convert P'' to standard representation and calculate residue
	p_dprime_std = convert_part_to_standard(nums, p_dprime, n);
	residue = karmarkar_karp(p_dprime_std, n);
	delete[] p_dprime_std;

	delete[] p_dprime;

	return residue;
}

long now_nano() {
	return std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now().time_since_epoch()).count();
}

int main(int argc, char* argv[]) {
	srand(time(0));

	if (argc != 4) {
		cout << "Error: incorrect number of command line arguments" << endl;
		return 0;
	}

	int n = 100;
	long max_iter = 25000;

	long debug = stoi(argv[1]);
	long algorithm = stoi(argv[2]);
	char *input_file = argv[3];

	if (debug == PRODUCTION) {
		//Initialize the input file
		FILE* file;
		file = fopen(input_file, "r");
		if (file == NULL) {
			printf("Error opening input file.\n");
			return 0;
		}

		//Loop and get each number in input file, store in nums array
		long s;
		long* nums = new long[n];
		for (int i1 = 0; i1 < n; i1++) {
			s = fscanf(file, "%ld", &nums[i1]);
			if (s != 1) {
				cout << "Error: cannot read file." << endl;
				return 0;
			}
		}

		//Find the residue from the specified method
		long residue;
		if (algorithm == KARMARKAR_KARP) {
			residue = karmarkar_karp(nums, n);
		} else if (algorithm == REPEATED_RANDOM) {
			residue = repeated_random(nums, n, max_iter);
		} else if (algorithm == HILL_CLIMBING) {
			residue = hill_climbing(nums, n, max_iter);
		} else if (algorithm == SIMULATED_ANNEALING) {
			residue = simulated_annealing(nums, n, max_iter);
		} else if (algorithm == PREPARTITIONED_REPEATED_RANDOM) {
			residue = repeated_random_part(nums, n, max_iter);
		} else if (algorithm == PREPARTITIONED_HILL_CLIMBING) {
			residue = hill_climbing_part(nums, n, max_iter);
		} else if (algorithm == PREPARTITIONED_SIMULATED_ANNEALING) {
			residue = simulated_annealing_part(nums, n, max_iter);
		} else {
			cout << "Error: Algorithm not recognized." << endl;
		}

		cout << residue << endl;

		delete[] nums;
	} else if (debug == SIMULATION) {
		//Variables for use in simulation
		long max = 1000000000000;
		int trials = 100;
		long *nums = new long[n];

		//Initialize random number generator
		uniform_int_distribution<long> intdist(1, max);

		//Variables to hold residual and times
		double kk_avg = 0;
		double kk_time = 0;
		double rr_avg = 0;
		double rr_time = 0;
		double hc_avg = 0;
		double hc_time = 0;
		double sa_avg = 0;
		double sa_time = 0;
		double rr_part_avg = 0;
		double rr_part_time = 0;
		double hc_part_avg = 0;
		double hc_part_time = 0;
		double sa_part_avg = 0;
		double sa_part_time = 0;

		//Variables for before and after timings of each algorithm
		long before, after;

		//Loop through the trials
		for (int trial = 0; trial < trials; trial++) {
			//Generate the random array
			for (int i1 = 0; i1 < n; i1++) {
				nums[i1] = intdist(generator);// random_uniform() * max;
			}

			//Run each algorithm and time it

			before = now_nano();
			kk_avg += karmarkar_karp(nums, n);
			after = now_nano();
			kk_time += after - before;

			before = now_nano();
			rr_avg += repeated_random(nums, n, max_iter);
			after = now_nano();
			rr_time += after - before;

			before = now_nano();
			hc_avg += hill_climbing(nums, n, max_iter);
			after = now_nano();
			hc_time += after - before;

			before = now_nano();
			sa_avg += simulated_annealing(nums, n, max_iter);
			after = now_nano();
			sa_time += after - before;

			before = now_nano();
			rr_part_avg += repeated_random_part(nums, n, max_iter);
			after = now_nano();
			rr_part_time += after - before;

			before = now_nano();
			hc_part_avg += hill_climbing_part(nums, n, max_iter);
			after = now_nano();
			hc_part_time += after - before;

			before = now_nano();
			sa_part_avg += simulated_annealing_part(nums, n, max_iter);
			after = now_nano();
			sa_part_time += after - before;

			cout << trial << endl;
		}
		delete[] nums;

		//Calculate the statistics

		kk_avg /= trials;
		rr_avg /= trials;
		hc_avg /= trials;
		sa_avg /= trials;
		rr_part_avg /= trials;
		hc_part_avg /= trials;
		sa_part_avg /= trials;

		kk_time = kk_time / trials / 1000000000;
		rr_time = rr_time / trials / 1000000000;
		hc_time = hc_time / trials / 1000000000;
		sa_time = sa_time / trials / 1000000000;
		rr_part_time = rr_part_time / trials / 1000000000;
		hc_part_time = hc_part_time / trials / 1000000000;
		sa_part_time = sa_part_time / trials / 1000000000;

		cout << "KK: " << kk_avg << ", " << kk_time << endl;
		cout << "RR: " << rr_avg << ", " << rr_time << endl;
		cout << "HC: " << hc_avg << ", " << hc_time << endl;
		cout << "SA: " << sa_avg << ", " << sa_time << endl;
		cout << "RR Part: " << rr_part_avg << ", " << rr_part_time << endl;
		cout << "HC Part: " << hc_part_avg << ", " << hc_part_time << endl;
		cout << "SA Part: " << sa_part_avg << ", " << sa_part_time << endl;
	} else if (debug == DEBUG) {
		
	}

	return 0;
}
