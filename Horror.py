
from fenics import *
from mshr import *

R=10

b=1

a=0.001



domain = Circle(Point(0, 0,0), R)

cylinder = Circle(Point(0, 0,0), b) - Circle(Point(0, 0,0), a)



mesh = generate_mesh(domain, 32)

#plot(mesh)

x = mesh.coordinates()[:,0]
y = mesh.coordinates()[:,1]
t = mesh.cells()

import pylab as plb

ax = plb.axes()
cm = plb.get_cmap('viridis')

ax.triplot(x, y, t, '-', color='k', lw=0.2, alpha=0.4)

# Output in the file
# print("Lenear_Poisson.pdf")
# plb.savefig("results/lenear_Poisson.%s" % "pdf", bbox_inches="tight")
print("Lenear_Poisson.png")
plb.savefig("results/Horror_new.%s" % "png", bbox_inches="tight")
