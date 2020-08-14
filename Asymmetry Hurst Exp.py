import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

def hurst(ts):
    ts = list(ts)
    N = len(ts)
    if N < 20:
        raise ValueError("Time series is too short! Input series ought to have at least 20 samples!")

    max_k = int(np.floor(N/2))
    R_S_dict = []
    for k in range(10,max_k+1):
        R,S = 0,0
        # split ts into subsets
        subset_list = [ts[i:i+k] for i in range(0,N,k)]
        if np.mod(N,k)>0:
            subset_list.pop()
            #tail = subset_list.pop()
            #subset_list[-1].extend(tail)
        # calc mean of every subset
        mean_list=[np.mean(x) for x in subset_list]
        for i in range(len(subset_list)):
            cumsum_list = pd.Series(subset_list[i]-mean_list[i]).cumsum()
            R += max(cumsum_list)-min(cumsum_list)
            S += np.std(subset_list[i])
        R_S_dict.append({"R":R/len(subset_list),"S":S/len(subset_list),"n":k})
    
    log_R_S = []
    log_n = []
    
    for i in range(301):
        R_S = (R_S_dict[i]["R"]+np.spacing(1)) / (R_S_dict[i]["S"]+np.spacing(1))
        log_R_S.append(np.log(R_S))
        log_n.append(np.log(R_S_dict[i]["n"]))

    Hurst_exponent = np.polyfit(log_n,log_R_S,1)[0]
    return Hurst_exponent, log_R_S, log_n

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
for i in range(len(sunspot_totalarea)):
    if sunspot_totalarea[i] == 0:
        AS.append(0)
    else:
        AS.append(sunspot_areadiff[i]/sunspot_totalarea[i])

Hexp, AS_H, AS_k = hurst(AS)

#linear_regression
coef = np.polyfit(AS_k, AS_H, 1)
poly1d_fn = np.poly1d(coef)

H = 10000*coef[0]
H = H - H%1
H = H/10000
RSquared = findRSquare(AS_H,poly1d_fn(AS_k))
H_error = 10000*(1-RSquared)*abs(H)
H_error = H_error - H_error%1
H_error = H_error/10000

#plot
plt.plot(AS_k, AS_H, marker='.', linewidth=0.5, markersize=2, c=(0,0,0))
plt.plot(AS_k, poly1d_fn(AS_k),'--k', linewidth=0.5, c=(0,0.7,1))

plt.plot([np.log(135),np.log(135)],[1,4],marker='.', linewidth=0.5, markersize=2, c=(0,0,1))
plt.plot([np.log(271),np.log(271)],[1,4],marker='.', linewidth=0.5, markersize=2, c=(1,0,0))
plt.text(3, max(AS_H)/1.2, '$H = $' + str(abs(H)) + u'$\u00B1$' + str(H_error))
plt.text(np.log(135), 2, r'$\tau = 135 months$', size=7.5)
plt.text(np.log(271), 2, r'$\tau = 271 months$', size=7.5)
         
xmin, xmax = plt.xlim()
ymin, ymax = plt.ylim()
plt.xlim(2, 6.3)
plt.ylim(1, 4)
plt.xlabel(r'${log(\tau)}$')
plt.ylabel(r'${log(R/S)}$')

plt.show()
