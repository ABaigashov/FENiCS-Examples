

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym

R0=0.0001

R1=10

r0=0.5

rho0=10

tol=1E-14

#mesh=RectangleMesh(Point(0,0), Point(1, np.pi*2), 50, 50, 'right/left')
domain=Rectangle(Point(R0,0), Point(R1, np.pi))

mesh=generate_mesh(domain, 128)

u_R=Expression(('0','1','0'),degree=2)

g=Constant(0.0)

P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2, v_3 = TestFunctions(V)

u = Function(V)

u_1, u_2, u_3 = split(u)

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0], R1, tol)

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)

bc_R = DirichletBC(V, u_R, boundary_R)

bcs=[bc_R]

f = Expression('x[0]*x[0]<=0.25 ? -10*(0.25-x[0]*x[0]) : 0', degree=2)

f0=Expression('1/(x[0])', degree=2)

f1=Expression('x[0]*x[0]*sin(x[1])',degree=2)

f2=Expression('x[0]*x[0]*x[0]*sin(x[1])',degree=2)

f3=Expression('x[0]',degree=2)

F=f1*u_1.dx(0)*v_1.dx(0)*dx+f1*f0*f0*u_1.dx(1)*v_1.dx(1)*dx \
-f1*v_1*u_3.dx(0)*u_1.dx(0)*dx \
+f2*u_2.dx(0)*v_2.dx(0)*dx+f2*f0*f0*u_2.dx(1)*v_2.dx(1)*dx \
+f3*u_3.dx(0)*v_3.dx(0)*dx+f3*f0*f0*u_3.dx(1)*v_3.dx(1)*dx \
-f3*v_3*u_1.dx(0)*u_1.dx(0)*dx \
-f1*f*u_2*u_2*exp(-2*u_1)*v_1*dx-f2*u_2*u_2*u_2*exp(-2*u_1)*f*v_2*dx \
-f3*f*u_2*u_2*exp(-2*u_1)*v_3*dx

#F=f2*u_1.dx(0)*v_1.dx(0)*dx+f2*f1*f1*u_1.dx(1)*v_1.dx(1)*dx \
#-f2*v_1*u_2.dx(0)*u_1.dx(0)*dx \
#+f3*u_2.dx(0)*v_2.dx(0)*dx+f3*f1*f1*u_2.dx(1)*v_2.dx(1)*dx \
#-f3*v_2*u_2.dx(0)*u_2.dx(0)*dx \
#-f2*f*exp(2*u_2-2*u_1)*v_1*dx-f3*exp(2*u_2-2*u_1)*f*v_2*dx

solve(F == 0, u, bcs)


# Hold plot

#задание координат точки
theta0=0

p=Point(R0,theta0)
F1=u(p.x(), p.y())[0]
F2=u(p.x(), p.y())[1]
F3=u(p.x(), p.y())[2]

print(F1,F2,F3)

data_x=[]

data_y1=[]

data_y2=[]

data_y3=[]


for i in range (100):
	x1=R0+0.001*i*(R1-R0)
	p=Point(x1,theta0)
	F1=u(p.x(), p.y())[0]
	F2=u(p.x(), p.y())[1]
	F3=u(p.x(), p.y())[2]
	data_x.append(x1)
	data_y1.append(F1)
	data_y2.append(F2)
	data_y3.append(F3)


fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y1)
plt.savefig('Profile1.png')
plt.show()

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y2)
plt.savefig('Profile2.png')
plt.show()

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y1)
plt.savefig('Profile1.png')
plt.show()

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y3)
plt.savefig('Profile3.png')
plt.show()


_u_1, _u_2,_u_3 = u.split()

vtkfile_u_1 = File('results/u_1.pvd')
vtkfile_u_2 = File('results/u_2.pvd')
vtkfile = File('results/u.pvd')


# Save solution to file

vtkfile_u_1 << _u_1
vtkfile_u_2 << _u_2

#plot(mesh)
