//CS124 Programming Assignemnt 2
//Qiang Fei and Jordan Turley
//
//Compile: gcc strassen.c -o strassen -O3 -std=c99
//Run:     ./strassen [debug code] [dimension] [input file]
//
//To run regularly, set [debug code] = 0. For example:
//	./strassen 0 4 input.txt
//To run with randomly generated matrices A and B, set [debug code] = 1 and [input file] = anything. For example:
//	./strassen 1 4 a
//To run simulation to find n0, set [debug code] = 2 and [input file] = anything. May need to recompile with -DUSE_CLOCK. For example:
//	./strassen 2 512 a
//To find triangles in random graph, set [debug code] = 3 and [input file] = anything. For example:
//	./strassen 3 1024 a

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

//Define debug codes
#define DEBUG 1
#define SIMULATION 2
#define GRAPH 3

//Define optimal crossover point
#define CROSSOVER 128

//Setup timing stuff
//Taken from timing code from CS205, timing.c and timing.h:
//https://harvard-iacs.github.io/2020-CS205/homeworks/HWB/
//Uncomment to use timing in simulation, may need to compile with -DUSE_CLOCK
/*
#if defined(CLOCK_HIGHRES)
#  define CLOCK CLOCK_HIGHRES
#elif defined(CLOCK_REALTIME)
#  define CLOCK CLOCK_REALTIME
#elif defined(USE_CLOCK)
#else
#  error No suitable clock found.  Check docs for clock_gettime.
#endif
*/

/**
 * Initializes a square 2D array of size n
 */
int** init_matrix(int n) {
	int** matrix = (int **) malloc(sizeof(int *) * n);
	for (int i1 = 0; i1 < n; i1++) {
		matrix[i1] = (int *) malloc(sizeof(int) * n);
	}
	return matrix;
}

/**
 * Frees up the memory of a given square 2D matrix of size n
 */
void free_matrix(int** matrix, int n) {
	for (int i1 = 0; i1 < n; i1++) {
		free(matrix[i1]);
	}
	free(matrix);
}

/**
 * Multiplies two square matrices, a and b, of size n using the naive
 * row-column method. Returns the product as a square matrix of size n.
 */
int** multiply_naive(int** a, int** b, int n) {
	int** result = init_matrix(n);

	int i1, i2, k;
	for (i1 = 0; i1 < n; i1++) {
		for (i2 = 0; i2 < n; i2++) {
			result[i1][i2] = 0;
			for (k = 0; k < n; k++) {
				result[i1][i2] += a[i1][k] * b[k][i2];
			}
		}
	}

	return result;
}

/**
 * Add (a + b) or subtract (a - b) two matrices of given size n. Adding
 * or subtracting is determined by bSign.
 */
int** add_matrix(int** a, int** b, int n, int bSign) {
	int** result = init_matrix(n);
	int i1, i2;
	for (i1 = 0; i1 < n; i1++) {
		for (i2 = 0; i2 < n; i2++) {
			result[i1][i2] = a[i1][i2] + b[i1][i2] * bSign;
		}
	}
	return result;
}

/**
 * Take a given matrix a of given size n and create a new matrix of size n / 2
 * consisting of the rows and columns of a given by the starting and ending indices.
 */
int** block_matrix(int** a, int n, int startRow, int endRow, int startCol, int endCol) {
	//Create the result matrix
	int n_2 = n / 2;
	int** result = init_matrix(n_2);

	//Loop over the correct part of a and store it in result
	int i1, i2, r = 0, c = 0;
	for (i1 = startRow; i1 < endRow; i1++) {
		c = 0;
		for (i2 = startCol; i2 < endCol; i2++) {
			result[r][c] = a[i1][i2];
			c++;
		}
		r++;
	}

	return result;
}

/**
 * Print out a matrix of a given size n to the console
 */
void print_matrix(int** matrix, int n) {
	int i1, i2;
	for (i1 = 0; i1 < n; i1++) {
		for (i2 = 0; i2 < n; i2++) {
			printf("%d ", matrix[i1][i2]);
		}
		printf("\n");
	}
}

/**
 * Multiplies two matrices a_ and b_ of size n_, using Strassen's algorithm
 * if the matrices are larger than the crossover point, or using the naive
 * algorithm if the matrices are smaller than the crossover point.
 */
int** multiply_matrix(int** a_, int** b_, int n_, int crossover) {
	//Once we are smaller than the crossover, multiply the naive way
	if (n_ <= crossover) {
		int** c = multiply_naive(a_, b_, n_);
		return c;
	}

	//Create matrices a, b, and size n
	//If the matrices can be evenly divided into blocks, a = a_, b = b_, and n = n_
	//If they cannot, we pad the matrix with a row and column of zeros and set a, b, and n
	//to the new matrix and new size
	int** a;
	int** b;
	int n;

	//If the matrix size is odd, we must pad the edges with zeros
	if (n_ % 2 == 0) {
		//n is divisible by two, 
		a = a_;
		b = b_;
		n = n_;
	} else {
		//Add a row and column of zeros to a and b
		int** anew = init_matrix(n_ + 1);
		int** bnew = init_matrix(n_ + 1);
		
		//Store the value of a_ and b_ in anew and bnew
		int i1, i2;
		for (i1 = 0; i1 < n_; i1++) {
			for (i2 = 0; i2 < n_; i2++) {
				anew[i1][i2] = a_[i1][i2];
				bnew[i1][i2] = b_[i1][i2];
			}
		}

		//Set the new row and column to zeros
		for (i1 = 0; i1 < n_ + 1; i1++) {
			anew[n_][i1] = 0;
			anew[i1][n_] = 0;
			bnew[n_][i1] = 0;
			bnew[i1][n_] = 0;
		}

		//Set a and b to the new padded matrices
		a = anew;
		b = bnew;

		n = n_ + 1;
	}

	int n_2 = n / 2;

	//Block matrices into 4 equal parts each
	int** a_11 = block_matrix(a, n, 0, n_2, 0, n_2);
	int** a_12 = block_matrix(a, n, 0, n_2, n_2, n);
	int** a_21 = block_matrix(a, n, n_2, n, 0, n_2);
	int** a_22 = block_matrix(a, n, n_2, n, n_2, n);

	int** b_11 = block_matrix(b, n, 0, n_2, 0, n_2);
	int** b_12 = block_matrix(b, n, 0, n_2, n_2, n);
	int** b_21 = block_matrix(b, n, n_2, n, 0, n_2);
	int** b_22 = block_matrix(b, n, n_2, n, n_2, n);

	//Make the additions and subtractions we need for the M matrices
	int** a_11_plus_a_22 = add_matrix(a_11, a_22, n_2, 1);
	int** b_11_plus_b_22 = add_matrix(b_11, b_22, n_2, 1);
	int** a_21_plus_a_22 = add_matrix(a_21, a_22, n_2, 1);
	int** b_12_minus_b_22 = add_matrix(b_12, b_22, n_2, -1);
	int** b_21_minus_b_11 = add_matrix(b_21, b_11, n_2, -1);
	int** a_11_plus_a_12 = add_matrix(a_11, a_12, n_2, 1);
	int** a_21_minus_a_11 = add_matrix(a_21, a_11, n_2, -1);
	int** b_11_plus_b_12 = add_matrix(b_11, b_12, n_2, 1);
	int** a_12_minus_a_22 = add_matrix(a_12, a_22, n_2, -1);
	int** b_21_plus_b_22 = add_matrix(b_21, b_22, n_2, 1);

	//Create the m_i matrices that will be used to calculate the c_ij matrices
	int** m_1 = multiply_matrix(a_11_plus_a_22, b_11_plus_b_22, n_2, crossover);
	int** m_2 = multiply_matrix(a_21_plus_a_22, b_11, n_2, crossover);
	int** m_3 = multiply_matrix(a_11, b_12_minus_b_22, n_2, crossover);
	int** m_4 = multiply_matrix(a_22, b_21_minus_b_11, n_2, crossover);
	int** m_5 = multiply_matrix(a_11_plus_a_12, b_22, n_2, crossover);
	int** m_6 = multiply_matrix(a_21_minus_a_11, b_11_plus_b_12, n_2, crossover);
	int** m_7 = multiply_matrix(a_12_minus_a_22, b_21_plus_b_22, n_2, crossover);

	//Free up the variables we don't need
	free_matrix(a_11, n_2);
	free_matrix(a_12, n_2);
	free_matrix(a_21, n_2);
	free_matrix(a_22, n_2);
	free_matrix(b_11, n_2);
	free_matrix(b_12, n_2);
	free_matrix(b_21, n_2);
	free_matrix(b_22, n_2);
	free_matrix(a_11_plus_a_22, n_2);
	free_matrix(b_11_plus_b_22, n_2);
	free_matrix(a_21_plus_a_22, n_2);
	free_matrix(b_12_minus_b_22, n_2);
	free_matrix(b_21_minus_b_11, n_2);
	free_matrix(a_11_plus_a_12, n_2);
	free_matrix(a_21_minus_a_11, n_2);
	free_matrix(b_11_plus_b_12, n_2);
	free_matrix(a_12_minus_a_22, n_2);
	free_matrix(b_21_plus_b_22, n_2);

	//Temp variables for forming C_ij matrices
	int** m_1_plus_m_4 = add_matrix(m_1, m_4, n_2, 1);
	int** m_1_plus_m_4_minus_m_5 = add_matrix(m_1_plus_m_4, m_5, n_2, -1);
	int** m_1_minus_m_2 = add_matrix(m_1, m_2, n_2, -1);
	int** m_1_minus_m_2_plus_m_3 = add_matrix(m_1_minus_m_2, m_3, n_2, 1);

	//Calculate the c_ij matrices that will form the result matrix c
	int** c_11 = add_matrix(m_1_plus_m_4_minus_m_5, m_7, n_2, 1);
	int** c_12 = add_matrix(m_3, m_5, n_2, 1);
	int** c_21 = add_matrix(m_2, m_4, n_2, 1);
	int** c_22 = add_matrix(m_1_minus_m_2_plus_m_3, m_6, n_2, 1);

	//Now that we have the c_ij matrices we don't need all these, free them
	free_matrix(m_1_plus_m_4, n_2);
	free_matrix(m_1_plus_m_4_minus_m_5, n_2);
	free_matrix(m_1_minus_m_2, n_2);
	free_matrix(m_1_minus_m_2_plus_m_3, n_2);
	free_matrix(m_1, n_2);
	free_matrix(m_2, n_2);
	free_matrix(m_3, n_2);
	free_matrix(m_4, n_2);
	free_matrix(m_5, n_2);
	free_matrix(m_6, n_2);
	free_matrix(m_7, n_2);

	int** c = init_matrix(n_);
	//Combine the c_ij matrices into the full c matrix
	int i1, i2;
	if (n == n_) {
		//Set c to each of the blocks c_ij
		for (i1 = 0; i1 < n_2; i1++) {
			for (i2 = 0; i2 < n_2; i2++) {
				c[i1][i2] = c_11[i1][i2];
				c[i1][i2 + n_2] = c_12[i1][i2];
				c[i1 + n_2][i2] = c_21[i1][i2];
				c[i1 + n_2][i2 + n_2] = c_22[i1][i2];
			}
		}
	} else {
		//We have to be more careful here to not include the zeros
		for (i1 = 0; i1 < n_2; i1++) {
			for (i2 = 0; i2 < n_2; i2++) {
				c[i1][i2] = c_11[i1][i2];
			}
		}
		for (i1 = 0; i1 < n_2; i1++) {
			for (i2 = n_2; i2 < n_; i2++) {
				c[i1][i2] = c_12[i1][i2 - n_2];
			}
		}
		for (i1 = n_2; i1 < n_; i1++) {
			for (i2 = 0; i2 < n_2; i2++) {
				c[i1][i2] = c_21[i1 - n_2][i2];
			}
		}
		for (i1 = n_2; i1 < n_; i1++) {
			for (i2 = n_2; i2 < n_; i2++) {
				c[i1][i2] = c_22[i1 - n_2][i2 - n_2];
			}
		}

		//Free up our new a and b that were padded with zeros
		free_matrix(a, n);
		free_matrix(b, n);
	}

	//Free up the c_ij matrices
	free_matrix(c_11, n_2);
	free_matrix(c_12, n_2);
	free_matrix(c_21, n_2);
	free_matrix(c_22, n_2);

	return c;
}

int main(int argc, char* argv[]) {
	//Ensure we have the correct command line arguments
	if (argc != 4) {
		printf("Please run as follows: ./strassen 0 [dimension] [inputfile]\n");
		printf("For debugging run: ./strassen 1 [dimension] [anything]");
		printf("Debugging will multiply two arbitrary n-dimensional matrices.");
		printf("The value of [inputfile] does not matter for debugging.");
		return 0;
	}

	//Get the debug code, dimension, and input file for the matrices
	int debug = atoi(argv[1]);
	int dimension = atoi(argv[2]);
	char* inputfile = argv[3];

	//Initialize the a and b matrices
	int** a = init_matrix(dimension);
	int** b = init_matrix(dimension);

	//Fill in a, b.
	//For debug, fill in randomly with 0 or 1
	//For graph, generate graph adjacency matrix
	//Otherwise, read from input file
	int i1, i2;
	if (debug == DEBUG || debug == SIMULATION) {
		srand(time(NULL));
		for (i1 = 0; i1 < dimension; i1++) {
			for (i2 = 0; i2 < dimension; i2++) {
				a[i1][i2] = rand() % 2;
				b[i1][i2] = rand() % 2;
			}
		}

		//print_matrix(a, dimension);
		//printf("\n");
		//print_matrix(b, dimension);
		//printf("\n");
	} else if (debug == 0) {
		//Read a and b in from the given input file
		FILE* file;
		file = fopen(inputfile, "r");
		if (file == NULL) {
			printf("Error opening input file.\n");
			return 0;
		}

		int s; //File read status code
		//Read the first half of the file in to a
		for (i1 = 0; i1 < dimension; i1++) {
			for (i2 = 0; i2 < dimension; i2++) {
				//Read the next int from the input file
				s = fscanf(file, "%d", &a[i1][i2]);
				if (s != 1) {
					printf("Error reading file.\n");
					return 0;
				}
			}
		}

		//Read the second half of the file in to b
		for (i1 = 0; i1 < dimension; i1++) {
			for (i2 = 0; i2 < dimension; i2++) {
				s = fscanf(file, "%d", &b[i1][i2]);
				if (s != 1) {
					printf("Error reading file.\n");
				}
			}
		}
	}

	if (debug == SIMULATION) { //Run the simulation with different values of n0
		/* Uncomment to run simulation. May need to compile with -DUSE_CLOCK
		struct timespec start, stop;
		int maxN0 = 512;
		int trials;
		for (int n0 = 2; n0 <= maxN0; n0 += 2) {
			long double avgTime = 0;
			for (trials = 0; trials < 5; trials++) {
				clock_gettime(CLOCK, &start);
				int** c = multiply_matrix(a, b, dimension, n0);
				clock_gettime(CLOCK, &stop);

				long double elapsed;
				elapsed = stop.tv_nsec - (double) start.tv_nsec;
				elapsed *= 1.0E-9L;
				elapsed += stop.tv_sec - (double) start.tv_sec;

				avgTime += elapsed;
				printf("n0 %d Took %Lf\n", n0, elapsed);
				free_matrix(c, dimension);
			}
			avgTime /= 5;
			printf("n0 %d avg time %Lf\n", n0, avgTime);
		}
		*/
	} else if (debug == GRAPH) { //Cube the adjacency matrix to find triangles in the graph
		srand(time(NULL));
		double p;
		double r;
		int trials = 500, trial;
		int triangles;
		for (p = 0.01; p <= 0.05; p += 0.01) {
			printf("p = %.2f\n", p);
			for (trial = 0; trial < trials; trial++) {
				//Generate adjacency matrix
				for (i1 = 0; i1 < dimension; i1++) {
					for (i2 = 0; i2 < i1; i2++) {
						r = (double) rand() / (double) RAND_MAX;
						if (r <= p) {
							a[i1][i2] = 1;
							a[i2][i1] = 1;
						} else {
							a[i1][i2] = 0;
							a[i2][i1] = 0;
						}
					}
				}

				int** a2 = multiply_matrix(a, a, dimension, CROSSOVER);
				int** a3 = multiply_matrix(a2, a, dimension, CROSSOVER);
				free_matrix(a2, dimension);

				triangles = 0;
				for (i1 = 0; i1 < dimension; i1++) {
					triangles += a3[i1][i1];
				}

				triangles /= 6;

				free_matrix(a3, dimension);
				printf("Triangles: %d\n", triangles);
			}
		}
	} else { //Compute a * b and print out the diagonal
		int** c = multiply_matrix(a, b, dimension, CROSSOVER);
		
		//print_matrix(c, dimension);
		//printf("\n");

		for (i1 = 0; i1 < dimension; i1++) {
			printf("%d\n", c[i1][i1]);
		}
		

		free_matrix(c, dimension);
	}

	free_matrix(a, dimension);
	free_matrix(b, dimension);

	return 0;
}
