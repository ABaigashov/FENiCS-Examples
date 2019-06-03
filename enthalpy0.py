import numpy as np
import sympy as smp
import math

def generate_src_function(lines,enthalpy0,X0,Y0):

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
    result=energy(X0,Y0,enthalpy0)
    return result



'''
F = open('data/EOS.txt', 'r')

lines = F.read().split()

F.close()

x = []
y= []
enth = []

F1=open('data/H1.txt','r')

lines1 = F1.read().split()

F1.close()

for i in range(0, len(lines1)-2, 3):
    enth.append(float(lines1[i+2]))
    x.append(float(lines1[i]))
    y.append(float(lines1[i+1]))



G = open('data/Table1.txt', 'w')
for i in range(len(x)):
    X=x[i]
    Y=y[i]
    enthalpy0=enth[i]
    RES=generate_src_function(lines,enthalpy0,X,Y)
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
