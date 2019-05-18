import fenics
import mshr

R=10

domain =  mshr.Rectangle(fenics.Point(-1., -1.), fenics.Point(1., 1.)) #\
         # mshr.Rectangle(fenics.Point(2., 1.25), fenics.Point(3., 1.75)) #\
         # - mshr.Circle(fenics.Point(1, 4), .25) \
         # + mshr.Circle(fenics.Point(7, 7), 1.25)

mesh = mshr.generate_mesh(domain, 16)

V = fenics.FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = fenics.Expression('x[0]*x[0] - x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
	return on_boundary

bc = fenics.DirichletBC(V, u_D, boundary)

# Define variational problem
u = fenics.TrialFunction(V)
v = fenics.TestFunction(V)
f = fenics.Constant(0.0)
a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
L = f*v*fenics.dx

# Compute solution
u = fenics.Function(V)
s = fenics.solve(a == L, u, bc)

print(type(s))
print(type(u))

#--------------------------Ploting pylab------------------------
#get array componets and triangulation :
v = u.compute_vertex_values(mesh)
x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

import pylab as plb

ax = plb.axes()
cm = plb.get_cmap('viridis')

ax.tricontourf(x, y, t, v, 10, cmap = cm)
ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

plb.ylim(-1, 1)
plb.xlim(-1, 1)

# Output in the file
print("gen_mesh.png")
plb.savefig("results/gen_mesh.%s" % "png", bbox_inches="tight")
