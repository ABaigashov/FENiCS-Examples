

from fenics import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt
import sympy as sym
import math
from NumSrcF import numer_src_function

R0=0.0001

R1=30


tol=1E-14

#mesh=RectangleMesh(Point(0,0), Point(1, np.pi*2), 50, 50, 'right/left')
domain=Rectangle(Point(R0,0), Point(R1, np.pi/2))

mesh=generate_mesh(domain, 64)

SRC = FunctionSpace(mesh, 'CG', 1)  # so that dofs are only in mesh vertices

n = SRC.dim()

print(n)

d = mesh.geometry().dim()

print(d)

F_dof_coordinates = SRC.tabulate_dof_coordinates().reshape(n,d)

F_dof_coordinates.resize((n, d))
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

g0=0.50

R0=5.0

x=sym.Symbol('x')
y=sym.Symbol('y')

g=g0*(1-x*x/(R0*R0))

G = open('data/H.txt', 'w')
for j in range(n):
    X=F_dof_x[j]
    Y=F_dof_y[j]
    enthalpy0=g.subs(x,X)
    enthalpy0=enthalpy0.subs(y,Y)
    a=[]
    a.append(X)
    a.append(Y)
    a.append(enthalpy0)
    G.write(str(a[0]))
    G.write(' ')
    G.write(str(a[1]))
    G.write(' ')
    G.write(str(a[2]))
    G.write(' ')
    G.write('\n')

G.close()

Name1='data/EOS.txt'
Name2='data/H.txt'

result=numer_src_function(Name1, Name2,n)

G = open('data/Src.txt', 'w')
for j in range(n):
	a=[]
	X=F_dof_x[j]
	density=result[0][j]
	pressure=result[1][j]
	Y=F_dof_y[j]

	a.append(X)
	a.append(Y)
	a.append(density)
	a.append(pressure)
	G.write(str(a[0]))
	G.write(' ')
	G.write(str(a[1]))
	G.write(' ')
	G.write(str(a[2]))
	G.write(' ')
	G.write(str(a[3]))
	G.write('\n')

G.close()

u_R=Expression(('0','0'),degree=2)

g=Constant(0.0)

P1 = FiniteElement('P', triangle, 1)
element = MixedElement([P1, P1])
V = FunctionSpace(mesh, element)

v_1, v_2 = TestFunctions(V)

u = Function(V)

u_1, u_2 = split(u)

def boundary_R(x, on_boundary):
	return on_boundary and near(x[0], R1, tol)

def boundary_L(x, on_boundary):
	return on_boundary and near(x[0], 0, tol)

bc_R = DirichletBC(V, u_R, boundary_R)

bcs=[bc_R]

'''
f = Expression('x[0]*x[0]<=0.25 ? 10*(0.25-x[0]*x[0]) : 0', degree=2)

f00=Expression('x[0]*x[0]<=0.25 ? 2*(0.25-x[0]*x[0]) : 0', degree=2)
'''

f = Function(SRC)

f00 = Function(SRC)

f.vector()[:] = result[0]
f00.vector()[:] = result[1]

f0=Expression('1/(x[0])', degree=2)

f1=Expression('x[0]*x[0]',degree=2)

f2=Expression('x[0]*x[0]*x[0]',degree=2)

f3=Expression('x[0]',degree=2)

F=-f1*u_1.dx(0)*v_1.dx(0)*dx-f1*f0*f0*u_1.dx(1)*v_1.dx(1)*dx \
+f1*v_1*u_1.dx(0)*u_2.dx(0)*dx \
-f2*u_2.dx(0)*v_2.dx(0)*dx-f2*f0*f0*u_2.dx(1)*v_2.dx(1)*dx \
+f2*v_2*u_2.dx(0)*u_2.dx(0)*dx \
-f1*4*(np.pi)*(f+3*f00)*exp(2*u_2-2*u_1)*v_1*dx \
-f2*exp(2*u_2-2*u_1)*4*4*(np.pi)*f00*v_2*dx

#F=f2*u_1.dx(0)*v_1.dx(0)*dx+f2*f1*f1*u_1.dx(1)*v_1.dx(1)*dx \
#-f2*v_1*u_2.dx(0)*u_1.dx(0)*dx \
#+f3*u_2.dx(0)*v_2.dx(0)*dx+f3*f1*f1*u_2.dx(1)*v_2.dx(1)*dx \
#-f3*v_2*u_2.dx(0)*u_2.dx(0)*dx \
#-f2*f*exp(2*u_2-2*u_1)*v_1*dx-f3*exp(2*u_2-2*u_1)*f*v_2*dx

solve(F == 0, u, bcs)


# Hold plot

#задание координат точки
theta0=0

data_x=[]

data_y1=[]

data_y2=[]

for i in range (1000):
	x1=0.0001+0.001*i*(R1)
	p=Point(x1,theta0)
	F1=u(p.x(), p.y())[0]
	F2=u(p.x(), p.y())[1]
	data_x.append(x1)
	data_y1.append(F1)
	data_y2.append(F2)



fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y1)
plt.savefig('results/Profile1.png')
plt.show()

fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.plot(data_x,data_y2)
plt.savefig('results/Profile2.png')
plt.show()


for k in range(3):

    p0=Point(0.0001,0)
    nu0=u(p0.x(), p0.y())[0]

    print(nu0)

    G = open('data/H1.txt', 'w')

    for i in range(n):
        a=[]
        p=Point(F_dof_x[i], F_dof_y[i])
        nu=u(p.x(), p.y())[0]
        H=g0+nu0-nu
        X=F_dof_x[i]
        Y=F_dof_y[i]
        a.append(H)
        G.write(str(X))
        G.write(' ')
        G.write(str(Y))
        G.write(' ')
        G.write(str(a[0]))
        G.write('\n')

    G.close()


    v_1, v_2 = TestFunctions(V)

    u = Function(V)

    u_1, u_2 = split(u)

    Name1='data/EOS.txt'
    Name2='data/H1.txt'

    result=numer_src_function(Name1, Name2,n)

    f = Function(SRC)

    f00 = Function(SRC)

    f.vector()[:] = result[0]
    f00.vector()[:] = result[1]

    F=-f1*u_1.dx(0)*v_1.dx(0)*dx-f1*f0*f0*u_1.dx(1)*v_1.dx(1)*dx \
    +f1*v_1*u_1.dx(0)*u_2.dx(0)*dx \
    -f2*u_2.dx(0)*v_2.dx(0)*dx-f2*f0*f0*u_2.dx(1)*v_2.dx(1)*dx \
    +f2*v_2*u_2.dx(0)*u_2.dx(0)*dx \
    -f1*4*(np.pi)*(f+3*f00)*exp(2*u_2-2*u_1)*v_1*dx \
    -f2*exp(2*u_2-2*u_1)*4*4*(np.pi)*f00*v_2*dx

    solve(F == 0, u, bcs)

    theta0=0

    data_x=[]

    data_y1=[]

    data_y2=[]

    for i in range (1000):
    	x1=0.0001+0.001*i*(R1)
    	p=Point(x1,theta0)
    	F1=u(p.x(), p.y())[0]
    	F2=u(p.x(), p.y())[1]
    	data_x.append(x1)
    	data_y1.append(F1)
    	data_y2.append(F2)


    fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
    plt.plot(data_x,data_y1)
    plt.savefig('results/Profile1-1.png')
    plt.show()

    fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
    plt.plot(data_x,data_y2)
    plt.savefig('results/Profile2-1.png')
    plt.show()










'''
_u_1, _u_2,_u_3 = u.split()

vtkfile_u_1 = File('results/u_1.pvd')
vtkfile_u_2 = File('results/u_2.pvd')
vtkfile = File('results/u.pvd')


# Save solution to file

vtkfile_u_1 << _u_1
vtkfile_u_2 << _u_2
'''
#plot(mesh)
