import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

"""Verifica della riproducibilità della radiazione trasmessa dal tubo."""

PATH = 'C:/Users/Lorenzo/Desktop/Lab/RX'
FILENAME = 'riproducibilita.txt'

PATH = os.path.join(PATH, FILENAME)
x1, x2, x3, x4, x5, x6, x7 = np.loadtxt(PATH, unpack = True)

V = np.array([40, 50, 60, 70, 80, 90, 100]) #kVp
N = 6
m1 = x1.sum()/len(x1)
m2 = x2.sum()/len(x2)
m3 = x3.sum()/len(x3)
m4 = x4.sum()/len(x4)
m5 = x5.sum()/len(x5)
m6 = x6.sum()/len(x6)
m7 = x7.sum()/len(x7)

s1 = np.sqrt(1/len(x1)*1/(len(x1)-1)*((x1 - m1)**2).sum())
s2 = np.sqrt(1/len(x2)*1/(len(x2)-1)*((x2 - m2)**2).sum())
s3 = np.sqrt(1/len(x3)*1/(len(x3)-1)*((x3 - m3)**2).sum())
s4 = np.sqrt(1/len(x4)*1/(len(x4)-1)*((x4 - m4)**2).sum())
s5 = np.sqrt(1/len(x5)*1/(len(x5)-1)*((x5 - m5)**2).sum())
s6 = np.sqrt(1/len(x6)*1/(len(x6)-1)*((x6 - m6)**2).sum())
s7 = np.sqrt(1/len(x7)*1/(len(x7)-1)*((x7 - m7)**2).sum())

err1 = np.sqrt((0.03*m1)**2 + s1**2)
err2 = np.sqrt((0.03*m2)**2 + s2**2)
err3 = np.sqrt((0.03*m3)**2 + s3**2)
err4 = np.sqrt((0.03*m4)**2 + s4**2)
err5 = np.sqrt((0.03*m5)**2 + s5**2)
err6 = np.sqrt((0.03*m6)**2 + s6**2)
err7 = np.sqrt((0.03*m7)**2 + s7**2)

m = [m1, m2, m3, m4, m5, m6, m7]
err = [err1, err2, err3, err4, err5, err6, err7]

if __name__ == '__main__':
    def retta(x, a, b):
        return a*x + b

    init_values = [1., 0.]

    pars, covm = curve_fit(retta, V, m, init_values, err)
    a0, b0 = pars
    da, db = np.sqrt(covm.diagonal())

    chisq = (((m - retta(V, *pars))/err)**2).sum()
    ndof = len(V) - 2

    print(f'a = {a0:.2f} +- {da:.2f} V/muGy\n')
    print(f'b = {b0:.2f} +- {db:.2f} muGy\n')
    print(f'chisq/ndof = {chisq:.3f}/{ndof}\n')

    xx = np.linspace(V[0], V[-1], 100)
    plt.subplot(2, 1, 1)
    plt.errorbar(V, m, err, marker = 'o', linestyle = '')
    plt.title('Riproducibilità')
    plt.ylabel('Dose [$\mu Gy$]')
    plt.plot(xx, retta(xx, *pars), color = 'r')

    plt.minorticks_on()
    plt.subplot(2, 1, 2)
    plt.plot(V, (m - retta(V, *pars))/err, marker = 'o', linestyle = '')
    plt.plot(V, np.zeros(len(V)), linestyle = '--')
    plt.xlabel('V [kVp]')
    plt.ylabel('Residui')
    plt.minorticks_on()

    plt.show()






