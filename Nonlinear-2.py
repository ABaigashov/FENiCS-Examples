

import fenics as fen

from mshr import *

import sympy as sym

x, y = sym.symbols('x[0], x[1]')

R=50

R0=10

rho0=0.1

domain=Circle(fen.Point(0,0),R)

mesh=generate_mesh(domain, 32)

V = fen.FunctionSpace(mesh, 'P', 1)

def boundary(x, on_boundary):
	return on_boundary

bc =fen.DirichletBC(V, fen.Constant(0), boundary)

f = fen.Expression('x[0]*x[0]+x[1]*x[1]<=100 ? 1 - (x[0]*x[0] + x[1]*x[1])/100 : 0', degree=3)

u = fen.Function(V)
v = fen.TestFunction(V)


F = fen.dot(fen.grad(u), fen.grad(v))*fen.dx-f*v*fen.dx+
# Plot solution

fen.solve(F == 0, u, bc)

fen.plot(u)
# Save solution to file

fen.File("results/U.pvd") << u
# Hold plot

#задание координат точки
p = fen.Point(0.0,0.0)

#вывод значения решения в точке
print(u(p.x(), p.y()))


#plot(mesh)
