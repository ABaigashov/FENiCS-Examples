

import fenics as fen

def q(u):
    return 1 + 0*u**2

import sympy as sym

x, y = sym.symbols('x[0], x[1]')

u = 2 + x + 0*2*y

f = - sym.diff(q(u)*sym.diff(u, x), x) - sym.diff(q(u)*sym.diff(u, y), y)

f = sym.simplify(f)

f=x+y

u_code = sym.printing.ccode(u)

f_code = sym.printing.ccode(f)

print('u =', u_code)
print('f =', f_code)



mesh = fen.UnitSquareMesh(8, 8)
V = fen.FunctionSpace(mesh, 'P', 1)


def boundary(x, on_boundary):
    return on_boundary



u = fen.Function(V)
v = fen.TestFunction(V)

u_D = fen.Expression(u_code, degree=1)

f = fen.Expression(f_code, degree=1)

bc =fen.DirichletBC(V, u_D, boundary)


F = q(u)*fen.dot(fen.grad(u), fen.grad(v))*fen.dx+u*u*v*fen.dx - f*v*fen.dx
# Plot solution

fen.solve(F == 0, u, bc)

fen.plot(u)
# Save solution to file

fen.File("results/U.pvd") << u
# Hold plot

#задание координат точки
p = fen.Point(0.8,0.2)

#вывод значения решения в точке
print(u(p.x(), p.y()))


#plot(mesh)
