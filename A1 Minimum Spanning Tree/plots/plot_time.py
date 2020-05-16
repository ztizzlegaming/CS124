import numpy as np
import matplotlib.pyplot as plt

def process(l):
	res = []
	for x in l:
		p = x.split('m')
		s = int(p[0]) * 60
		s += float(p[1][:-1])
		res.append(s)
	return res

n = list(range(7, 19))
g1 = ['0m0.004s', '0m0.006s', '0m0.007s', '0m0.015s', '0m0.038s', '0m0.104s', '0m0.338s', '0m1.279s', '0m4.835s', '0m19.522s', '1m13.332s', '4m41.436s']
g2 = ['0m0.005s', '0m0.006s', '0m0.008s', '0m0.013s', '0m0.022s', '0m0.059s', '0m0.146s', '0m0.439s', '0m1.539s', '0m5.789s', '0m20.366s', '1m18.514s']
g3 = ['0m0.004s', '0m0.005s', '0m0.007s', '0m0.013s', '0m0.026s', '0m0.064s', '0m0.164s', '0m0.527s', '0m1.881s', '0m6.751s', '0m27.188s', '1m45.231s']
g4 = ['0m0.004s', '0m0.005s', '0m0.007s', '0m0.013s', '0m0.028s', '0m0.071s', '0m0.190s', '0m0.620s', '0m2.145s', '0m7.675s', '0m30.936s', '1m59.095s']

g1 = process(g1)
g2 = process(g2)
g3 = process(g3)
g4 = process(g4)

plt.plot(n, g1)
plt.plot(n, g2)
plt.plot(n, g3)
plt.plot(n, g4)
plt.xlabel('Number of vertices ($2^x$)')
plt.ylabel('Runtime (seconds)')
plt.title('Running Time to find MST for One Trial vs. Number of Vertices')
plt.legend(['Uniform', '2-Dimensional', '3-Dimensional', '4-Dimensional'])
plt.show()
