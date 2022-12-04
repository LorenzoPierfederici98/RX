import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import riprod as rip

"""Verifica dell'andamento lineare della dose in aria al variare di corrente anodica,
tempo di ritardo per ogni kVp."""

PATH = 'C:/Users/Lorenzo/Desktop/Lab/RX'
FILENAME = '60kVp.txt'

PATH = os.path.join(PATH, FILENAME)

I = np.array([10., 20., 32., 40., 50., 63., 80., 100.]) #mA
t = np.array([10., 20., 32., 40., 50., 80., 100.])

D_I = np.loadtxt(PATH, max_rows = 8) #dose a tempo fissato a 100ms
D_t = np.loadtxt(PATH, unpack = True, skiprows = 9)  #dose a corrente fissata a 100mA
N = 6 #numero di misure per kVp in riprod.py

#SWITCH 0 tempo fissato e considero D_I e err_I
#SWITCH 1 corrente fissata e considero D_t e err_t
SWITCH = 1

FILENAME = FILENAME.replace('.txt', '')
if FILENAME == '40kVp':
    s = rip.s1
elif FILENAME == '60kVp':
    s = rip.s3
elif FILENAME == '80kVp':
    s = rip.s5
    SWITCH = 0 #80 e 100 kVp non hanno misure a tempi diversi
elif FILENAME == '100kVp':
    s = rip.s7
    SWITCH = 0

#3% + fluttuazioni della singola misura
err_I = np.sqrt((0.03*D_I)**2 + N*s**2)
err_t = np.sqrt((0.03*D_t)**2 + N*s**2)

def retta(x, a, b):
    """Funzione di fit lineare."""
    return a*x + b

init_values = [1., 0.]

plt.subplot(2,1,1)
if SWITCH == 0:
    plt.title('V = ' + FILENAME + '  ' + 't = 100 ms')
    x = I
    y = D_I
    dy = err_I
    xx = np.linspace(I[0], I[-1], 100)
    UNIT = 'mA'
elif SWITCH == 1:
    plt.title('V = ' + FILENAME + '  ' + 'I = 100 mA')
    if FILENAME == '60kVp':
        t = np.insert(t, 5, 63.)
    x = t
    y = D_t
    dy = err_t
    xx = np.linspace(t[0], t[-1], 100)
    UNIT = 'ms'

pars, covm = curve_fit(retta, x, y, init_values, dy)
chisq = (((y - retta(x, *pars))/dy)**2).sum()
ndof = len(x) - 2

a0, b0 = pars
da, db = np.sqrt(covm.diagonal())
print(f'a = {a0:.3f} +- {da:.3f} muGy/{UNIT}\n')
print(f'b = {b0:.3f} +- {db:.3f} muGy\n')
print(f'chisq/ndof = {chisq:.3f}/{ndof}\n')

plt.ylabel('Dose [$mu Gy$]')
plt.errorbar(x, y, dy, marker = 'o', linestyle = '')
plt.plot(xx, retta(xx, *pars), color = 'r')
plt.minorticks_on()

plt.subplot(2,1,2)
plt.plot(x, (y - retta(x, *pars))/dy, marker = 'o', linestyle = '--')
plt.plot(x, np.zeros(len(x)))
plt.ylabel('Residui')
if SWITCH == 0:
    plt.xlabel('I [mA]')
elif SWITCH == 1:
    plt.xlabel('t [ms]')
plt.minorticks_on()

plt.show()









