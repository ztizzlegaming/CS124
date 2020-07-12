# Programming Assignment 2 - Strassen Algorithm
## Harvard University, CS 124
### Qiang Fei, Jordan Turley

We compared the Strassen algorithm with the naive matrix multiplication algorithm for varying matrix sizes.

### Compile

    gcc strassen.c -o strassen

### Run

Debugging mode, multiply two random matrices, where dimensions is the size of matrices and anything is literally anything (0, asdf, whatever):

    ./strassen 1 [dimension] [anything]

Regular run, read matrices of a given size from a given input file:

    ./strassen 0 [dimension] [input file]
