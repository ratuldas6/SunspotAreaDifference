import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates

data = pd.read_csv('area_m_1874-2016.csv')

data['Date'] = pd.to_datetime(data['Date'])
data.sort_values('Date', inplace=True)

sunspot_date = data['Date']
month_nos = [x for x in range(1,len(sunspot_date)+1)] 
sunspot_areadiff = data['Difference']
sunspot_totalarea = data['Total']
AS = []
for i in range(len(sunspot_totalarea)):
    if sunspot_totalarea[i] == 0:
        AS.append(0)
    else:
        AS.append(sunspot_areadiff[i]/sunspot_totalarea[i])

#plot1
plt.subplot(211)
plt.plot(month_nos, sunspot_areadiff, linestyle='-',linewidth=0.5,marker='.',ms=3, c=(0,0,1))
plt.ylabel('Asymmetry (N-S)',fontsize=13)
plt.tight_layout()
#plot2
plt.subplot(212)
plt.plot(month_nos, AS, linestyle='-',linewidth=0.5,marker='.',ms=3, c=(1,0,0))
'''linestyle='solid', marker='None'''
'''markersize = 3, c=(0,0,0),'''
scale_factor = 1.5
ymin, ymax = plt.ylim()
plt.ylim(ymin * scale_factor, ymax * scale_factor)
plt.xlabel('Time (months)',fontsize=13)
plt.ylabel('Normalised Asymmetry',fontsize=13)

plt.tight_layout()
plt.show()
