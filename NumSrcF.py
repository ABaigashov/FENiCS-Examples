from fenics import *
import numpy as np
from mshr import *
import sympy as smp
from enthalpy0 import generate_src_function
import matplotlib.pyplot as plt


# solve -laplace(u) = f in [0, 1]^2 with homog. Dirichlet boundary conditions
# and funny 'numerical' rhs

#mesh = RectangleMesh(-1,-1,1,1,10,10,'right')

def numer_src_function(Name1,Name2,n):
    # source term calculation

    f_values = np.zeros(n)

    g_values=np.zeros(n)

    G = open(Name1, 'r')

    lines = G.read().split()

    G.close()

    F1=open(Name2,'r')

    lines1 = F1.read().split()

    F1.close()

    x = []
    y= []
    enth = []

    for i in range(0, len(lines1)-2, 3):
        enth.append(float(lines1[i+2]))
        x.append(float(lines1[i]))
        y.append(float(lines1[i+1]))

    for i in range(n):
        enthalpy0=enth[i]
        result=generate_src_function(lines,enthalpy0,x[i],y[i])
        f_values[i]=result[2]
        g_values[i]=result[3]

    return f_values, g_values

'''
domain=Rectangle(Point(0,0), Point(10, 10))

mesh=generate_mesh(domain, 16)

F = FunctionSpace(mesh, 'CG', 1)  # so that dofs are only in mesh vertices

n = F.dim()

print(n)

d = mesh.geometry().dim()

print(d)

f = Function(F)

f1 = Function(F)

    # dofmap approach
F_dof_coordinates = F.tabulate_dof_coordinates().reshape(n,d)

F_dof_coordinates.resize((n, d))
F_dof_x = F_dof_coordinates[:, 0]
F_dof_y = F_dof_coordinates[:, 1]

g0=0.45

R0=5

x=smp.Symbol('x')
y=smp.Symbol('y')

g=g0*(1-x*x/(R0*R0))



N=10

X0=1
Y0=0

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

result=numer_src_function(Name1, Name2)

print(result[0][104])


f.vector()[:] = result[0]
f1.vector()[:] = result[1]




fig=plt.figure(figsize=(8,6), facecolor='pink', frameon=True)
plt.scatter(F_dof_x,result[0])
plt.scatter(F_dof_x,result[1])
plt.savefig('results/Src.png')
plt.show()

'''
# # vertex to dofmap approach
# vertex_values = np.zeros(mesh.num_vertices())
# for vertex in vertices(mesh):
#   x = vertex.x(0)
#   y = vertex.x(1)
#   vertex_values[vertex.index()] = abs(x + y)
#
# f.vector()[:] = vertex_values[dof_to_vertex_map(F)]
# plot(f, title='Source term')

# variational problem definition

'''
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

a = inner(grad(u), grad(v))*dx
L = f*v*dx
bc = DirichletBC(V, Constant(0.), DomainBoundary())

u = Function(V)
solve(a == L, u, bc)

plot(u, title='Solution')
'''
#interactive()
