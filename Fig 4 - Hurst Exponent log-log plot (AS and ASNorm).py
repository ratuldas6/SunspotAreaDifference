import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def findRSquare(y,y_new):
    n = len(y)
    sumsqreg = 0
    totalsumsq = 0
    mean_y = 0
    for i in range(n):
        pasq = (y[i]-y_new[i])**2
        sumsqreg += pasq
        mean_y += y[i]
    mean_y = mean_y/n
    for j in range(n):
        masq = (y[i]-mean_y)**2
        totalsumsq += masq
    rsq = 1 - (sumsqreg/totalsumsq)
    return rsq

data = pd.read_csv('logrs_vs_logk.csv')

logk_total = data['log_s']
logrs_total = data['log_rs']
logrs_total_diff = data['log_rs_diff']

#get r/s values for k = 11 to 1699
logk = []
logrs = []
logrsd = []
for i in range(1699):
    logk.append(logk_total[i])
    logrs.append(logrs_total[i])
    logrsd.append(logrs_total_diff[i])
#build list for linear regressions for ASNorm
logk1 = []
logrs1 = []
logk2 = []
logrs2 = []
for i in range(130):
    logk1.append(logk_total[i])
    logrs1.append(logrs_total[i])
for i in range(223,631):
    logk2.append(logk_total[i])
    logrs2.append(logrs_total[i])
#build list for linear regressions for AS
logk1d = []
logrs1d = []
logk2d = []
logrs2d = []
for i in range(111):
    logk1d.append(logk_total[i])
    logrs1d.append(logrs_total_diff[i])
for i in range(263,651):
    logk2d.append(logk_total[i])
    logrs2d.append(logrs_total_diff[i])

#regression
coef1 = np.polyfit(logk1, logrs1, 1)
poly1d_fn1 = np.poly1d(coef1)
coef2 = np.polyfit(logk2, logrs2, 1)
poly1d_fn2 = np.poly1d(coef2)
coef3 = np.polyfit(logk1d, logrs1d, 1)
poly1d_fn3 = np.poly1d(coef3)
coef4 = np.polyfit(logk2d, logrs2d, 1)
poly1d_fn4 = np.poly1d(coef4)

H1 = 10000*coef1[0]
H1 = H1 - H1%1
H1 = H1/10000
RSquared1 = findRSquare(logrs1,poly1d_fn1(logk1))
H_error1 = 10000*(1-RSquared1)*abs(H1)
H_error1 = H_error1 - H_error1%1
H_error1 = H_error1/10000
c1 = coef1[1]

H2 = 10000*coef2[0]
H2 = H2 - H2%1
H2 = H2/10000
RSquared2 = findRSquare(logrs2,poly1d_fn2(logk2))
H_error2 = 10000*(1-RSquared2)*abs(H2)
H_error2 = H_error2 - H_error2%1
H_error2 = H_error2/10000
c2 = coef2[1]

H3 = 10000*coef3[0]
H3 = H3 - H3%1
H3 = H3/10000
RSquared3 = findRSquare(logrs1d,poly1d_fn3(logk1d))
H_error3 = 10000*(1-RSquared3)*abs(H3)
H_error3 = H_error3 - H_error3%1
H_error3 = H_error3/10000
c3 = coef3[1]

H4 = 10000*coef4[0]
H4 = H4 - H4%1
H4 = H4/10000
RSquared4 = findRSquare(logrs2d,poly1d_fn4(logk2d))
H_error4 = 10000*(1-RSquared4)*abs(H4)
H_error4 = H_error4 - H_error4%1
H_error4 = H_error4/10000
c4 = coef4[1]

#plot
plt.subplot(211)
plt.plot(logk, logrsd, marker='.', linewidth=1, markersize=1, c=(0,0,1), label=r'${AS}$')
plt.text(3.5, 1, r'$(a)$', size=20)
#regression lines
plt.plot([0.75,4], [H3*0.75+c3,H3*4+c3],'--k', linewidth=1, color='crimson', label=r'$\tau = 11-110$'+' months')
plt.plot([0.75,4], [H4*0.75+c4,H4*4+c4],'--k', linewidth=1, color='darkred', label=r'$\tau = 263-650$'+' months')
plt.text(1.2, 0.78, '$H = $' + str(H3)+ u'$\u00B1$' + str(H_error3),size=12, rotation=18)
plt.text(2.54, 1.40, '$H = $' + str(H4)+ u'$\u00B1$' + str(H_error4),size=12, rotation=19)
#vertical temporal markers
plt.plot([np.log10(135),np.log10(135)],[0,3],marker='.', linewidth=0.5, markersize=2, c=(0,0,0))
plt.plot([np.log10(271),np.log10(271)],[0,3],marker='.', linewidth=0.5, markersize=2, c=(0,0,0))
plt.text(2.15, 0.04, r'$\tau = 135$'+' months', size=10, rotation=90)
plt.text(2.45, 0.04, r'$\tau = 271$'+' months', size=10, rotation=90)
#limit window
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.xlim(0.74, 4)
plt.ylim(0, 3)
#labels and legend
plt.ylabel(r'${log(R/S)}$',fontsize=11)
leg = plt.legend()
plt.legend(loc='upper left',fontsize=14,shadow=True)

#2nd plot
plt.subplot(212)
plt.plot(logk, logrs, marker='.', linewidth=1, markersize=1, c=(1,0,0), label=r'${AS_{Norm}}$')
plt.text(3.5, 1, r'$(b)$', size=20)
#regression lines
plt.plot([0.75,4], [H1*0.75+c1,H1*4+c1],'--k', linewidth=1, color='purple', label=r'$\tau = 11-130$'+' months')
plt.plot([0.75,4], [H2*0.75+c2,H2*4+c2],'--k', linewidth=1, color='forestgreen', label=r'$\tau = 223-630$'+' months')
plt.text(1.1, 0.7, '$H = $' + str(H1) + u'$\u00B1$' + str(H_error1),size=12, rotation=19)
plt.text(2.5, 1.4, '$H = $' + str(H2)+ u'$\u00B1$' + str(H_error2),size=12, rotation=20)
#vertical temporal markers
plt.plot([np.log10(135),np.log10(135)],[0,3],marker='.', linewidth=0.5, markersize=2, c=(0,0,0))
plt.plot([np.log10(271),np.log10(271)],[0,3],marker='.', linewidth=0.5, markersize=2, c=(0,0,0))
plt.text(2.15, 0.04, r'$\tau = 135$'+' months', size=10, rotation=90)
plt.text(2.45, 0.04, r'$\tau = 271$'+' months', size=10, rotation=90)
#limit window
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.xlim(0.75, 4)
plt.ylim(0, 3)
#labels and legend
plt.xlabel(r'${log(\tau)}$',fontsize=13)
plt.ylabel(r'${log(R/S)}$',fontsize=13)
leg = plt.legend()
plt.legend(loc='upper left',fontsize=14,shadow=True)
plt.tight_layout()

plt.show()
