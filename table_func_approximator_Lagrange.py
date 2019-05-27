from scipy.interpolate import lagrange
from numpy.polynomial.polynomial import Polynomial
import sympy as sym

def approximator(f):
    """ Функция возвращает, аппроксимированную полиномами Лагранжа,
        таблично-заданную функцию y = f(x), в терминах C++.

        На вход функции подается открытый для чтения файл формата .txt
        (f = open('table_func.txt', 'r'))
    """
    # Чтение и разбиение строк списка
    lines = f.read().split()
    f.close()

    # Создание списков для абсциссы и ординаты
    abscissa = []
    ordinate = []
    for i in range(0, len(lines)-1, 2):
        abscissa.append(float(lines[i]))
        ordinate.append(float(lines[i+1]))

    # Представление функции полиномом и определение коэффициентов полинома
    poly = lagrange(abscissa, ordinate)
    poly_coef = Polynomial(poly).coef

    # Определение степени полинома
    N = len(poly_coef)

    # Составление полинома
    x = sym.symbols('x[0]')
    func = 0
    for i in range(0, N, 1):
        func = func + poly_coef[N-i-1]*x**i

    # Принтинг функции в язык С++
    f_code = sym.printing.ccode(func)

    return f_code

#
# print(approximator(f))
