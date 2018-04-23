import math
import matplotlib.pyplot as plt



quick = [0.0344403, 0.097614, 0.118179, 0.247941, 0.220659, 0.298444, 0.424821, 0.44849, 0.580806, 1.07838]
radix = [0.1417, 0.125746, 0.142445, 0.153742, 0.241969, 0.548323, 0.417442, 0.44637, 0.50140, 0.66169]


plt.xlabel("Size of input array)")
plt.ylabel("Running time in seconds")

t = linspace(1024, 524288, 10)

plt.plot(t, quick, 'r', label='Quick-Sort') # plotting t, a separately
plt.plot(t, radix, 'b', label='Radix-Sort') # plotting t, b separately

plt.legend(loc="best")

plt.savefig("result.png")
plt.show()

