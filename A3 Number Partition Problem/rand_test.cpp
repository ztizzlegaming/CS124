#include <iostream>
#include <random>
//#include <time.h>

using namespace std;

mt19937 generator(time(0));
uniform_real_distribution<double> distribution(0.0, 1.0);

void x() {
	cout << distribution(generator) << endl;
	cout << distribution(generator) << endl;
	cout << distribution(generator) << endl;
	cout << distribution(generator) << endl;
	cout << distribution(generator) << endl;
}


int main() {
	x();
	cout << endl;
	x();
	return 0;
}