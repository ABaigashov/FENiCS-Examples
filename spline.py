import scipy as sp
import matplotlib.pyplot as plt
import sympy as sym
def mnkGP(x,y):
              d=2 # степень полинома
              fp, residuals, rank, sv, rcond = sp.polyfit(x, y, d, full=True) # Модель
              f = sp.poly1d(fp) # аппроксимирующая функция
              # print('Коэффициент -- a %s  '%round(fp[0],4))
              # print('Коэффициент-- b %s  '%round(fp[1],4))
              # print('Коэффициент -- c %s  '%round(fp[2],4))
              y1=[fp[0]*x[i]**2+fp[1]*x[i]+fp[2] for i in range(0,len(x))] # значения функции a*x**2+b*x+c
              so=round(sum([abs(y[i]-y1[i]) for i in range(0,len(x))])/(len(x)*sum(y))*100,4) # средняя ошибка
              # print('Average quadratic deviation '+str(so))
              fx = sp.linspace(x[0], x[-1] + 1, len(x)) # можно установить вместо len(x) большее число для интерполяции
              # plt.plot(x, y, 'o', label='Original data', markersize=10)
              # plt.plot(fx, f(fx), linewidth=2)
              # plt.grid(True)
              # plt.show()

x=[10, 14, 18, 22, 26, 30, 34, 38, 42, 46, 50, 54, 58, 62, 66, 70, 74, 78, 82, 86]
y=[0.1, 0.0714, 0.0556, 0.0455, 0.0385, 0.0333, 0.0294, 0.0263, 0.0238, 0.0217,
   0.02, 0.0185, 0.0172, 0.0161, 0.0152, 0.0143, 0.0135, 0.0128, 0.0122,
   0.0116] # данные для проверки по функции y=1/x

#f_code = sym.printing.ccode(mnkGP(x,y))
print(mnkGP(x,y))
