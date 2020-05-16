import matplotlib.pyplot as plt

f = open('graph_out_500.txt')
d = {}
per = 500
for i1 in range(5):
	l = f.readline().strip()
	p = l.split(' = ')[1]
	triangles = []
	for i2 in range(per):
		l = f.readline().strip()
		tri = int(l.split(': ')[1])
		triangles.append(tri)
	d[p] = triangles

choose = 178433024
for p in (0.01, 0.02, 0.03, 0.04, 0.05):
	triangles = d[str(p)]
	cumAvg = []
	avg = 0
	for i1 in range(per):
		avg += triangles[i1]
		cumAvg.append(avg / (i1 + 1))

	tr = choose * (p ** 3)

	print(p, cumAvg[-1], tr, cumAvg[-1] - tr)

	x = range(per)
	plt.close()
	plt.plot(x, cumAvg)
	plt.axhline(choose * (p ** 3), c = 'r')
	plt.xlabel('Number of Trials Averaged')
	plt.ylabel('Number of Triangles')
	plt.title('Number of Triangles vs. Number of Trials Averaged for $p = ' + str(p) + '$')
	plt.show()
