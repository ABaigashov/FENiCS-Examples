from scipy.interpolate import splev, splrep
import sympy as sym

def f(x):
    x_points = [ 0, 1, 2, 3, 4, 5]
    u_points = [12,14,22,39,58,77]

    tmp = splrep(x_points, y_points)

    return splev(x, tmp)



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
