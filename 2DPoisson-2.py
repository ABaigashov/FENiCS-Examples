
from fenics import *
from mshr import *
import numpy as np

R=50

R0=10

rho0=0.1

tol=1E-14

#mesh=RectangleMesh(Point(0,0), Point(1, np.pi*2), 50, 50, 'right/left')
domain=Rectangle(Point(0.0001,0), Point(1, np.pi))

mesh=generate_mesh(domain, 64)

u_R=Constant(0)

g=Constant(0.0)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0], 1, tol)

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)

bc_R = DirichletBC(V, u_R, boundary_R)

#bc_L=DirichletBC(V, u_L, boundary_L)

#bcs = [bc_L, bc_R]

bcs=[bc_R]

# Define variational problem
u = Function(V)
v = TestFunction(V)

f = Expression('x[0]*x[0]<=0.25 ? 10: 0', degree=2)

f1=Expression('1/(x[0])', degree=2)

f2=Expression('x[0]*x[0]*x[0]',degree=2)

f3=Expression('x[0]',degree=2)

F=f2*f1*f1*u.dx(1)*v.dx(1)*dx+f2*u.dx(0)*v.dx(0)*dx-f2*f*v*dx+g*v*ds

solve(F == 0, u, bcs)

#задание координат точки
p = Point(0.001,0)

#вывод значения решения в точке
print(u(p.x(), p.y()))


# Save to file and plot
File("results/Solution-2.pvd") << u
