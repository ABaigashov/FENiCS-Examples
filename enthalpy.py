import numpy as np
import sympy as smp
import math

def generate_src_function(lines,f,X0,Y0):
    x=smp.Symbol('x')
    y=smp.Symbol('y')
    density = []
    pressure = []
    enthalpy = []


    for i in range(0, len(lines)-2, 3):
        enthalpy.append(float(lines[i]))
        density.append(float(lines[i+1]))
        pressure.append(float(lines[i+2]))

    def energy(X0,Y0,enthalpy0):
        for i in range(len(pressure)):
            if enthalpy0<0:
                result1=0
                result2=0
                break
            else:
                for j in range(len(pressure)):
                    if enthalpy0<enthalpy[j]:
                        K1=(pressure[j]-pressure[j-1])/(math.exp(enthalpy[j])-math.exp(enthalpy[j-1]))
                        K2=(density[j]-density[j-1])/(math.exp(enthalpy[j])-math.exp(enthalpy[j-1]))
                        result1=pressure[j-1]+(math.exp(enthalpy0)-math.exp(enthalpy[j-1]))*K1
                        result2=density[j-1]+(math.exp(enthalpy0)-math.exp(enthalpy[j-1]))*K2
                        break
        return X0,Y0,result2,result1

    enthalpy0=f.subs(x,X0)
    enthalpy0=enthalpy0.subs(y,Y0)
    result=energy(X0,Y0,enthalpy0)
    return result

'''
Name = 'data/EOS.txt'
x=smp.Symbol('x')
y=smp.Symbol('y')
f0=0.45

R0=5

f=f0*(1-x*x/(R0*R0))

N=10

X0=1
Y0=0

F = open(Name, 'r')

lines = F.read().split()

F.close()

bum=generate_src_function(lines,f,X0,Y0)

print(bum)

X1=2
Y1=1
bum=generate_src_function(lines,f,X1,Y1)
print(bum)
'''

'''
G = open('data/Table.txt', 'w')
for i in range(N):
    X=0.001+i*R0/N
    RES=generate_src_function(lines,f,X,Y0)
    a=[]
    a.append(RES[0])
    a.append(RES[1])
    a.append(RES[2])
    a.append(RES[3])
    G.write(str(a[0]))
    G.write(' ')
    G.write(str(a[1]))
    G.write(' ')
    G.write(str(a[2]))
    G.write(' ')
    G.write(str(a[3]))
    G.write('\n')

G.close()
'''
