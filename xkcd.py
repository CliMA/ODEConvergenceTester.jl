import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 1, 100)
mu = 0.7
sigma = 0.1
y = np.exp(-0.5*((x - mu)/sigma)**2)

plt.xkcd(scale=5, length=400)
plt.xticks([])
plt.yticks([])
plt.ylabel('Some variable')
plt.xlabel('Time')
plt.title('A poor choice of `tspan`')
plt.text(0, 0.3, 'Convergence study window')
plt.plot([0.2, 0.1], [0.25, 0.15], 'red')
plt.plot([0.0, 0.2], [0.1, 0.1], 'black') # horizontal line
plt.plot([0.0, 0.0], [0.05, 0.15], 'black') # vertical bar
plt.plot([0.2, 0.2], [0.05, 0.15], 'black') # vertical bar
plt.text(0, 0.8, 'Interesting phenomena')
plt.plot(x, y)
plt.plot([0.2, 0.5], [0.7, 0.5], 'red')
plt.savefig('xkcd_plot.png', dpi=300)
