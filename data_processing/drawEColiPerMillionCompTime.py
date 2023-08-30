# Importing the required modules
import matplotlib.pyplot as plt
import numpy as np

# Creating a new figure
plt.figure(dpi=300)

# Creating an array from the list
x= np.array([100,200,300,400,500,600,700,800,900,1000])
y1= np.array([1.143, 2.23, 3.15, 4.33, 5.298, 6.41234, 7.43769, 8.45738, 9.47989, 10.36])
y2= np.array([7.207, 19.93, 37.89, 61.766, 90.336, 123.698, 162.625, 206.105, 255.068, 309.23])

# Plotting Numpy array
plt.yticks(np.arange(10, 360, 50))
plt.plot(x, y1, label ='Linear Jaro Algorithm')
plt.plot(x, y2, '-.', label ='Quadratic Jaro Algorithm')

# plt.yscale("ln")
# Adding details to the plot
# plt.title('Plot NumPy array')
plt.xlabel('String Lengths')
plt.ylabel('Elapsed Time (seconds)')
plt.legend()

# Displaying the plot
# plt.show()

plt.savefig('Figure_1.png')