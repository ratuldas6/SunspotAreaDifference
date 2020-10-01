import pandas as pd
import numpy as np
from matplotlib import pyplot as plt

def findCurveLength(series,m,k):
    add_terms = (len(series)-m)/k
    add_terms = int(add_terms - add_terms%1)
    N = add_terms
    norm = (len(series)-1)/(add_terms*k)
    series.insert(0,0)
    L = 0
    for i in range(0,N+1):
        L += (abs(series[m+i*k]-series[m+(i-1)*k])*norm)/k
    del series[0]
    return L

def avgCL(series,k):
    aCL = 0
    for i in range(1,k+1):
        aCL+=findCurveLength(series,i,k)
    aCL = aCL/k
    return aCL

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

data = pd.read_csv('area_m_1874-2016.csv')

sunspot_areadiff = data['Difference']
sunspot_totalarea = data['Total']
AS = []
ASD = []
for i in range(len(sunspot_totalarea)):
    ASD.append(sunspot_areadiff[i])
    if sunspot_totalarea[i] == 0:
        AS.append(0)
    else:
        AS.append(sunspot_areadiff[i]/sunspot_totalarea[i])

AS_avgCLk = []
ASD_avgCLk = []
AS_k = []

for i in range(2,55):
    AS_avgCLk.append(np.log(avgCL(AS,i)))
    ASD_avgCLk.append(np.log(avgCL(ASD,i)))
    AS_k.append(np.log(i))
    
#linear_regression
coef = np.polyfit(AS_k, AS_avgCLk, 1)
poly1d_fn = np.poly1d(coef)
c = coef[1]

D = 10000*coef[0]
D = D - D%1
D = D/10000
RSquared = findRSquare(AS_avgCLk,poly1d_fn(AS_k))
D_error = 10000*(1-RSquared)*abs(D)
D_error = D_error - D_error%1
D_error = D_error/10000

coef1 = np.polyfit(AS_k, ASD_avgCLk, 1)
poly1d_fn1 = np.poly1d(coef1)
c1 = coef1[1]

D1 = 10000*coef1[0]
D1 = D1 - D1%1
D1 = D1/10000
RSquared1 = findRSquare(ASD_avgCLk,poly1d_fn1(AS_k))
D1_error = 10000*(1-RSquared1)*abs(D1)
D1_error = D1_error - D1_error%1
D1_error = D1_error/10000

#plots
plt.plot(AS_k, ASD_avgCLk, marker='.', linestyle='None', markersize=4, c=(0,0,1), label=r'${AS}$')
plt.plot(AS_k, AS_avgCLk, marker='.', linestyle='None', markersize=4, c=(1,0,0), label=r'${AS_{Norm}}$')
plt.plot([0,5], [c,D*5+c],'--k', linewidth=0.7, c=(1,0,0))
plt.plot([0,5], [c1,D1*5+c1],'--k', linewidth=0.7, c=(0,0,1))

plt.text(2.5, 0, '$D = $' + str(abs(D)) + u'$\u00B1$' + str(D_error),size=13,rotation=-23)
plt.text(2.5, 6.7, '$D = $' + str(abs(D1)) + u'$\u00B1$' + str(D1_error),size=13,rotation=-23)

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.xlim(0, 4.5)
plt.ylim(-2, 15)
plt.xlabel(r'${log(\tau)}$',fontsize=13)
plt.ylabel(r'${log\langle{L(\tau)}\rangle}$',fontsize=13)

#legend
leg = plt.legend()
plt.legend(loc='upper right',fontsize=16,shadow=True)
plt.tight_layout()

plt.show()

