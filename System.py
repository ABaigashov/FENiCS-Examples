

from fenics import *
from mshr import *
import numpy as np

R=50

R0=10

rho0=0.1

tol=1E-14

#mesh=RectangleMesh(Point(0,0), Point(1, np.pi*2), 50, 50, 'right/left')
domain=Rectangle(Point(0.001,0), Point(1, np.pi))

mesh=generate_mesh(domain, 64)

u_R=Expression(('0','0'),degree=2)

g=Constant(0.0)

P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)

u_1, u_2 = split(u)

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0], 1, tol)

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)

bc_R = DirichletBC(V, u_R, boundary_R)

bcs=[bc_R]

f = Expression('x[0]*x[0]<=1 ? 1 : 0', degree=2)

f1=Expression('1/(x[0])', degree=2)

f2=Expression('x[0]',degree=2)

f3=Expression('x[0]*x[0]*sin(x[1])',degree=2)

F=f2*u_1.dx(0)*v_1.dx(0)*dx+f2*f1*f1*u_1.dx(1)*v_1.dx(1)*dx \
+f3*u_2.dx(0)*v_2.dx(0)*dx+f3*f1*f1*u_2.dx(1)*v_2.dx(1)*dx \
-f2*f*v_1*dx-f3*f*v_2*dx

solve(F == 0, u, bcs)


# Hold plot

#задание координат точки
p = Point(0.001,0.5)

#вывод значения решения в точке
print(u(p.x(), p.y()))

vtkfile_u_1 = File('results/u_1.pvd')
vtkfile_u_2 = File('results/u_2.pvd')

# Save solution to file

vtkfile_u_1 << u_1
vtkfile_u_2 << u_2
#plot(mesh)
