import matplotlib.pyplot as plt

n = 512
f = open('timing_' + str(n) + '_full.txt')
lines = f.readlines()
f.close()

times = []
for i1 in range(0, len(lines), 6):
	l = lines[i1 + 5].strip()
	time = l.split()[-1]
	times.append(float(time))

n0s = tuple(range(2, len(times) * 2 + 1, 2))

plt.plot(n0s, times)
plt.axvline(128, c = 'r', alpha = 0.2)
plt.xlabel('$n_0$')
plt.ylabel('Average Runtime, 5 Trials')
plt.yscale('log')
plt.title('Average Runtime vs. $n_0$ for Matrix Size $n = ' + str(n) + '$')
plt.show()
