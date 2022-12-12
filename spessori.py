import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import riprod as rip

# spessori in cm e relative incertezze
spessori = np.array([i*0.120 for i in range(7)])
dspessori = np.array([i*0.005 for i in range(7)])


print(dspessori)
# dosex corrispondenti a x spessori
# dosex è un array di dosi, una per ogni V = [40, 60, 70, 80, 90, 100] kVp, I = 100 mA e t = 100 ms
# e ripetuto per più spessori
dose1 = np.loadtxt('spessori.txt', max_rows=6, unpack=True)
dose2 = np.loadtxt('spessori.txt', skiprows=7, max_rows=6, unpack=True)
dose3 = np.loadtxt('spessori.txt', skiprows=14, max_rows=6, unpack=True)
dose4 = np.loadtxt('spessori.txt', skiprows=21, max_rows=6, unpack=True)
dose5 = np.loadtxt('spessori.txt', skiprows=28, max_rows=6, unpack=True)
dose6 = np.loadtxt('spessori.txt', skiprows=35, max_rows=6, unpack=True)
dose_vuoto = np.loadtxt('spessori.txt', skiprows=42, unpack=True)

# dose_x è la dose a kVp fissato, da 40 a 100 kVp, al variare degli spessori
dose_1 = np.array([dose_vuoto[0], dose1[0], dose2[0],
                  dose3[0], dose4[0], dose5[0], dose6[0]])
dose_2 = np.array([dose_vuoto[1], dose1[1], dose2[1],
                  dose3[1], dose4[1], dose5[1], dose6[1]])
dose_3 = np.array([dose_vuoto[2], dose1[2], dose2[2],
                  dose3[2], dose4[2], dose5[2], dose6[2]])
dose_4 = np.array([dose_vuoto[3], dose1[3], dose2[3],
                  dose3[3], dose4[3], dose5[3], dose6[3]])
dose_5 = np.array([dose_vuoto[4], dose1[4], dose2[4],
                  dose3[4], dose4[4], dose5[4], dose6[4]])
dose_6 = np.array([dose_vuoto[5], dose1[5], dose2[5],
                  dose3[5], dose4[5], dose5[5], dose6[5]])

# incertezze sulle dosi, vedi lin.py
s_40 = rip.s1
s_60 = rip.s3
s_70 = rip.s4
s_80 = rip.s5
s_90 = rip.s6
s_100 = rip.s7

N = 6
err_1 = np.sqrt((0.03*dose_1)**2 + N*s_40**2)
err_2 = np.sqrt((0.03*dose_2)**2 + N*s_60**2)
err_3 = np.sqrt((0.03*dose_3)**2 + N*s_70**2)
err_4 = np.sqrt((0.03*dose_4)**2 + N*s_80**2)
err_5 = np.sqrt((0.03*dose_5)**2 + N*s_90**2)
err_6 = np.sqrt((0.03*dose_6)**2 + N*s_100**2)

r_Al = 2.70  # densità alluminio g/cm^3
r_Pb = 11.34  # densità piombo g/cm^3

mu_Al = 0.073215  # coeff. att. alluminio cm^2/g
mu_Pb = 0.10675  # coeff. att. piombo cm^2/g

t = r_Al*spessori  # t g/cm^2
dt = r_Al*dspessori


def esponenziale(x, a, b):
    return a*np.exp(-b*x)


init_values1 = [dose_1[0], 0.5685]
init_values2 = [dose_2[0], 0.2778]
init_values3 = [dose_3[0], 0.2398]
init_values4 = [dose_4[0], 0.2018]
init_values5 = [dose_5[0], 0.1861]
init_values6 = [dose_6[0], 0.1704]

pars1, covm1 = curve_fit(esponenziale, t, dose_1, init_values1, err_1)
pars2, covm2 = curve_fit(esponenziale, t, dose_2, init_values2, err_2)
pars3, covm3 = curve_fit(esponenziale, t, dose_3, init_values3, err_3)
pars4, covm4 = curve_fit(esponenziale, t, dose_4, init_values4, err_4)
pars5, covm5 = curve_fit(esponenziale, t, dose_5, init_values5, err_5)
pars6, covm6 = curve_fit(esponenziale, t, dose_6, init_values6, err_6)

a1, b1 = pars1
da1, db1 = np.sqrt(covm1.diagonal())
a2, b2 = pars2
da2, db2 = np.sqrt(covm2.diagonal())
a3, b3 = pars3
da3, db3 = np.sqrt(covm3.diagonal())
a4, b4 = pars4
da4, db4 = np.sqrt(covm4.diagonal())
a5, b5 = pars5
da5, db5 = np.sqrt(covm5.diagonal())
a6, b6 = pars6
da6, db6 = np.sqrt(covm6.diagonal())

dy1 = np.sqrt((-pars1[1]*esponenziale(t, a1, b1)*dt)**2 + err_1**2)
dy2 = np.sqrt((-pars2[1]*esponenziale(t, a2, b2)*dt)**2 + err_2**2)
dy3 = np.sqrt((-pars3[1]*esponenziale(t, a3, b3)*dt)**2 + err_3**2)
dy4 = np.sqrt((-pars4[1]*esponenziale(t, a4, b4)*dt)**2 + err_4**2)
dy5 = np.sqrt((-pars5[1]*esponenziale(t, a5, b5)*dt)**2 + err_5**2)
dy6 = np.sqrt((-pars6[1]*esponenziale(t, a6, b6)*dt)**2 + err_6**2)

pars1, covm1 = curve_fit(esponenziale, t, dose_1, init_values1, dy1)
pars2, covm2 = curve_fit(esponenziale, t, dose_2, init_values2, dy2)
pars3, covm3 = curve_fit(esponenziale, t, dose_3, init_values3, dy3)
pars4, covm4 = curve_fit(esponenziale, t, dose_4, init_values4, dy4)
pars5, covm5 = curve_fit(esponenziale, t, dose_5, init_values5, dy5)
pars6, covm6 = curve_fit(esponenziale, t, dose_6, init_values6, dy6)

chisq1 = (((dose_1 - esponenziale(t, *pars1))/dy1)**2).sum()
chisq2 = (((dose_2 - esponenziale(t, *pars2))/dy2)**2).sum()
chisq3 = (((dose_3 - esponenziale(t, *pars3))/dy3)**2).sum()
chisq4 = (((dose_4 - esponenziale(t, *pars4))/dy4)**2).sum()
chisq5 = (((dose_5 - esponenziale(t, *pars5))/dy5)**2).sum()
chisq6 = (((dose_6 - esponenziale(t, *pars6))/dy6)**2).sum()
ndof = len(dose_1) - 2


print(f'a = {a1:.3f} +- {da1:.3f} muGy\t b = {b1:.3f} +- {db1:.3f} cm^2/g\t chisq/ndof = {chisq1:.2f}/{ndof}\t V = 40 kVp\n')
print(f'a = {a2:.3f} +- {da2:.3f} muGy\t b = {b2:.3f} +- {db2:.3f} cm^2/g\t chisq/ndof = {chisq2:.2f}/{ndof}\t V = 60 kVp\n')
print(f'a = {a3:.3f} +- {da3:.3f} muGy\t b = {b3:.3f} +- {db3:.3f} cm^2/g\t chisq/ndof = {chisq3:.2f}/{ndof}\t V = 70 kVp\n')
print(f'a = {a4:.3f} +- {da4:.3f} muGy\t b = {b4:.3f} +- {db4:.3f} cm^2/g\t chisq/ndof = {chisq4:.2f}/{ndof}\t V = 80 kVp\n')
print(f'a = {a5:.3f} +- {da5:.3f} muGy\t b = {b5:.3f} +- {db5:.3f} cm^2/g\t chisq/ndof = {chisq5:.2f}/{ndof}\t V = 90 kVp\n')
print(f'a = {a6:.3f} +- {da6:.3f} muGy\t b = {b6:.3f} +- {db6:.3f} cm^2/g\t chisq/ndof = {chisq6:.2f}/{ndof}\t V = 100 kVp\n')


HVL1 = np.log(2)/b1
dHVL1 = np.log(2)/b1**2 * db1
HVL2 = np.log(2)/b2
dHVL2 = np.log(2)/b2**2 * db2
HVL3 = np.log(2)/b3
dHVL3 = np.log(2)/b3**2 * db3
HVL4 = np.log(2)/b4
dHVL4 = np.log(2)/b4**2 * db4
HVL5 = np.log(2)/b5
dHVL5 = np.log(2)/b5**2 * db5
HVL6 = np.log(2)/b6
dHVL6 = np.log(2)/b6**2 * db6

print(f'HVL = {10*HVL1/r_Al:.2f} +- {10*dHVL1/r_Al:.2f} mm\t HVL = {10*HVL2/r_Al:.2f} +- {10*dHVL2/r_Al:.2f} mm\t HVL = {10*HVL3/r_Al:.2f} +- {10*dHVL3/r_Al:.2f} mm\n')
print(f'HVL = {10*HVL4/r_Al:.2f} +- {10*dHVL4/r_Al:.2f} mm\t HVL = {10*HVL5/r_Al:.2f} +- {10*dHVL5/r_Al:.2f} mm\t HVL = {10*HVL6/r_Al:.2f} +- {10*dHVL6/r_Al:.2f} mm\n')
xx = np.linspace(t[0], t[-1], 100)


fig, axs = plt.subplots(2, 3)
fig.suptitle('Fit degli spessori massici vs dosi misurate al variare dei kVp')
axs[0, 0].set_title('V = 40 kVp')
axs[0, 0].errorbar(t, dose_1, err_1, dt, marker='o', linestyle='')
axs[0, 0].plot(xx, esponenziale(xx, *pars1), color='r')
axs[0, 0].text(
    0.3, 0.7, f'$D_0 = ({a1:.1f} \pm {da1:.1f}) \mu Gy$\n $\mu = ({b1:.2f} \pm {db1:.2f}) cm^2/g$\n $\chi/ndof = {chisq1:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[0, 0].transAxes)

axs[0, 1].set_title('V = 60 kVp')
axs[0, 1].errorbar(t, dose_2, err_2, dt, marker='o', linestyle='')
axs[0, 1].plot(xx, esponenziale(xx, *pars2), color='r')
axs[0, 1].text(
    0.4, 0.7, f'$D_0 = ({a2:.1f} \pm {da2:.1f}) \mu Gy$\n $\mu = ({b2:.2f} \pm {db2:.2f}) cm^2/g$\n $\chi/ndof = {chisq2:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[0, 1].transAxes)


axs[0, 2].set_title('V = 70 kVp')
axs[0, 2].errorbar(t, dose_3, err_3, dt, marker='o', linestyle='')
axs[0, 2].plot(xx, esponenziale(xx, *pars3), color='r')
axs[0, 2].text(
    0.4, 0.7, f'$D_0 = ({a3:.1f} \pm {da3:.1f}) \mu Gy$\n $\mu = ({b3:.2f} \pm {db3:.2f}) cm^2/g$\n $\chi/ndof = {chisq3:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[0, 2].transAxes)

axs[1, 0].set_title('V = 80 kVp')
axs[1, 0].errorbar(t, dose_4, err_4, dt, marker='o', linestyle='')
axs[1, 0].plot(xx, esponenziale(xx, *pars4), color='r')
axs[1, 0].text(
    0.4, 0.7, f'$D_0 = ({a4:.1f} \pm {da4:.1f}) \mu Gy$\n $\mu = ({b4:.2f} \pm {db4:.2f}) cm^2/g$\n $\chi/ndof = {chisq4:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[1, 0].transAxes)

axs[1, 1].set_title('V = 90 kVp')
axs[1, 1].errorbar(t, dose_5, err_5, dt, marker='o', linestyle='')
axs[1, 1].plot(xx, esponenziale(xx, *pars5), color='r')
axs[1, 1].text(
    0.4, 0.7, f'$D_0 = ({a5:.1f} \pm {da5:.1f}) \mu Gy$\n $\mu = ({b5:.2f} \pm {db5:.2f}) cm^2/g$\n $\chi/ndof = {chisq5:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[1, 1].transAxes)

axs[1, 2].set_title('V = 100 kVp')
axs[1, 2].errorbar(t, dose_6, err_6, dt, marker='o', linestyle='')
axs[1, 2].plot(xx, esponenziale(xx, *pars6), color='r')
axs[1, 2].text(
    0.4, 0.7, f'$D_0 = ({a6:.1f} \pm {da6:.1f}) \mu Gy$\n $\mu = ({b6:.2f} \pm {db6:.2f}) cm^2/g$\n $\chi/ndof = {chisq6:.1f}/{ndof}$', bbox=dict(facecolor='white', alpha=0.2), fontsize='medium', transform=axs[1, 2].transAxes)

for ax in axs.flat:
    ax.set(xlabel='t [$g/cm^2$]', ylabel='Dose [$\mu$Gy]')

plt.show()
