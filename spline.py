from scipy.interpolate import splev, splrep, InterpolatedUnivariateSpline, interp1d
import sympy as sym
import matplotlib.pyplot as plt
import numpy as np


x_points = [0, 0.3, 0.4, 0.1, 0.6, 0.7]
y_points = [1., 2., 5., 10., 99., 88]
x = np.linspace(0, 10, num=11, endpoint=True)
y = np.cos(-x**2/9.0)
f = interp1d(x, y)
f2 = interp1d(x, y, kind='cubic')

xnew = np.linspace(0, 10, num=41, endpoint=True)

print(f(xnew))

plt.plot(x, y, 'o', xnew, f(xnew), '-', xnew, f2(xnew), '--')
plt.legend(['data', 'linear', 'cubic'], loc='best')
plt.show()


# x = sym.symbols('x[0]')
# f = - sym.diff(u*sym.diff(u, x), x)
#u_code = sym.printing.ccode(u)



# import sympy as sym
#
# # Use SymPy to compute f from the manufactured solution u
# x = sym.symbols('x[0]')
# u = 1 + x**2
# f = - sym.diff(u*sym.diff(u, x), x)
# print("f = ", f)
#
# f = sym.simplify(f)
# print("f = ", f)
#
# f_code = sym.printing.ccode(f)
# print("f = ", f_code)
