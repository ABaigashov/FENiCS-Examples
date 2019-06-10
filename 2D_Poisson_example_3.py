
from fenics import *
from mshr import *
import numpy as np

# Пример решения линейного уравнения Пуассона,
# записанный без использования встроенных операторов интегрирования и
# дифференцирования вариационной задачив прямоугольной области с
# аналитически заданной граничной функцией и произвольно сгенерированной сеткой.
# Вывод решения осуществляется в двух направлениях vtkfile и Ploting pylab

R=1
R0=0.1
rho0=0.1

domain=Circle(Point(0,0),R)
mesh=generate_mesh(domain, 32)

plot(mesh)

V = FunctionSpace(mesh, 'P', 1)

# Define boundary condition

def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Constant(0), boundary)

# Define variational problem
u = TrialFunction(V)
v = TestFunction(V)

f = Expression('x[0]*x[0]+x[1]*x[1]<=0.01 ? 1 : 0', degree=2)

b=u.dx(1)*v.dx(1)*dx+u.dx(0)*v.dx(0)*dx

L=f*v*dx

# Compute solution
u = Function(V)
solve(b == L, u, bc)

#задание координат точки
p = Point(0.5,0)

#вывод значения решения в точке
print(u(p.x(), p.y()))

# Save to file and plot
File("results/Solution-1.pvd") << u
