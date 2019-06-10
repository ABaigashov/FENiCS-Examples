import fenics as fen

# Пример решения уравнения Пуассона в прямоугольной области с
# аналитически заданной граничной функцией.
# Вывод решения осуществляется в двух направлениях vtkfile и Ploting pylab

# Create mesh and define function space
mesh = fen.UnitSquareMesh(8, 8)
V = fen.FunctionSpace(mesh, 'P', 1)

# Define boundary condition
u_D = fen.Expression('1 + x[0]*x[0] + x[1]*x[1]', degree=2)

def boundary(x, on_boundary):
	return on_boundary

bc = fen.DirichletBC(V, u_D, boundary)

# Define variational problem
u = fen.TrialFunction(V)
v = fen.TestFunction(V)
f = fen.Constant(-6.0)
a = fen.dot(fen.grad(u), fen.grad(v))*fen.dx
L = f*v*fen.dx

# Compute solution
u = fen.Function(V)
fen.solve(a == L, u, bc)


# #--------------------------Ploting PVD------------------------
# # Plot solution and mesh
# plot(u)
# plot(mesh)
#
# # Save solution to file in VTK format
# vtkfile = File('Lenear_Poisson/Lenear_Poisson.pvd')
# vtkfile << u
#
# # Compute maximum error at vertices
# vertex_values_u_D = u_D.compute_vertex_values(mesh)
# vertex_values_u = u.compute_vertex_values(mesh)



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

# Output in the file
# print("Lenear_Poisson.pdf")
# plb.savefig("results/lenear_Poisson.%s" % "pdf", bbox_inches="tight")
print("Lenear_Poisson.png")
plb.savefig("results/lenear_Poisson.%s" % "png", bbox_inches="tight")
