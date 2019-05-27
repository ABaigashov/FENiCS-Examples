from scipy.interpolate import lagrange
import numpy as np

x = np.array([0, 1, 2])
y = x**3
poly = lagrange(x, y)
