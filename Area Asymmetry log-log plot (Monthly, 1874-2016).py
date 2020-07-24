import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import rc as rcee

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

plt.style.use('seaborn')

data = pd.read_csv('area_m_1874-2016.csv')

sunspot_areadiff = data['Difference']
sunspot_totalarea = data['Total']
AS = []
for i in range(len(sunspot_totalarea)):
    if sunspot_totalarea[i] == 0:
        AS.append(0)
    else:
        AS.append(sunspot_areadiff[i]/sunspot_totalarea[i])

AS_avgCLk = []
AS_k = []

for i in range(2,33):
    AS_avgCLk.append(np.log(avgCL(AS,i))/np.log(2))
    AS_k.append(np.log(i)/np.log(2))
    
#linear_regression
coef = np.polyfit(AS_k, AS_avgCLk, 1)
poly1d_fn = np.poly1d(coef)

D = 10000*coef[0]
D = D - D%1
D = D/10000
RSquared = findRSquare(AS_avgCLk,poly1d_fn(AS_k))
D_error = 10000*(1-RSquared)*abs(D)
D_error = D_error - D_error%1
D_error = D_error/10000

plt.plot(AS_k, AS_avgCLk, marker='.', linestyle='None', c=(0,0,0))
#plt.plot(AS_k, AS_avgCLk, 'yo', AS_k, poly1d_fn(AS_k),'--k')
plt.plot(AS_k, poly1d_fn(AS_k),'--k', linewidth=0.5, c=(0,0.7,1))
plt.text(3, 6, 'D = ' + str(abs(D)) + u'\u00B1' + str(D_error))

xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.xlim(0, 5.5)
plt.ylim(-1, 9)
plt.title('Hemispherical sunspot asymmetry log-log plot (1874-2016)')
plt.xlabel('log2(k)')
plt.ylabel('log2(<L(k)>)')

plt.show()
