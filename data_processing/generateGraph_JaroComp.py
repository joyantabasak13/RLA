import numpy as np
import matplotlib.pyplot as plt


x_full = np.linspace(0.001, 0.02, 20)
x_enn = np.linspace(0.05, 1.95, 20)
x_knn = np.linspace(2, 40, 20)

y_full = np.random.rand(len(x_full))
y_enn = np.random.rand(len(x_enn))
y_knn = np.random.rand(len(x_knn))


fig, ax_full = plt.subplots()

full = ax_full.plot(x_full, y_full, color = 'red', label = 'full')
ax_full.spines['bottom'].set_color('red')
ax_full.tick_params(axis='x', colors='red')
ax_full.xaxis.label.set_color('red')


ax_enn = ax_full.twiny()
enn = ax_enn.plot(x_enn, y_enn, color = 'blue', label = 'enn')
ax_enn.xaxis.set_ticks_position('bottom')
ax_enn.xaxis.set_label_position('bottom')
ax_enn.spines['bottom'].set_position(('axes', -0.15))
ax_enn.spines['bottom'].set_color('blue')
ax_enn.tick_params(axis='x', colors='blue')
ax_enn.xaxis.label.set_color('blue')


ax_knn = ax_full.twiny()
knn = ax_knn.plot(x_knn, y_knn, color = 'green', label = 'knn')
ax_knn.xaxis.set_ticks_position('bottom')
ax_knn.xaxis.set_label_position('bottom')
ax_knn.spines['bottom'].set_position(('axes', -0.3))
ax_knn.spines['bottom'].set_color('green')
ax_knn.tick_params(axis='x', colors='green')
ax_knn.xaxis.label.set_color('green')


lines = full + enn + knn
labels = [l.get_label() for l in lines]
ax_full.legend(lines, labels)

plt.tight_layout()

plt.show()