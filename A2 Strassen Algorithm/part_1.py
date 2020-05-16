import numpy as np
import matplotlib.pyplot as plt

def T(n, n0):
	if n > n0:
		return 7 * T(n / 2, n0) + 18 * (n / 2) ** 2
	else:
		return n * n * (2 * n - 1)

n = 1057
n0Vals = tuple(range(1, n + 1))
tVals = []
for n0 in n0Vals:
	tVals.append(T(n, n0))

minIdx = np.argmin(tVals)
print(n0Vals[minIdx])

plt.plot(n0Vals, tVals)
plt.xlabel('$n_0$')
plt.ylabel('T(n), Number of Arithmetic Operations')
plt.title('Arithmetic Operations vs. $n_0$ for n = ' + str(n))
plt.show()
