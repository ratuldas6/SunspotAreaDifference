import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates

plt.style.use('seaborn')

data = pd.read_csv('area_m_1874-2016.csv')

data['Date'] = pd.to_datetime(data['Date'])
data.sort_values('Date', inplace=True)

sunspot_date = data['Date']
sunspot_areadiff = data['Difference']
sunspot_totalarea = data['Total']
AS = []
for i in range(len(sunspot_totalarea)):
    if sunspot_totalarea[i] == 0:
        AS.append(0)
    else:
        AS.append(sunspot_areadiff[i]/sunspot_totalarea[i])

plt.plot_date(sunspot_date, AS, markersize = 3, c=(0,0,0))

plt.gcf().autofmt_xdate()

scale_factor = 3
ymin, ymax = plt.ylim()
plt.ylim(ymin * scale_factor, ymax * scale_factor)
plt.title('Monthly Hemispherical Sunspot Area Difference by Total (1874-2016)')
plt.xlabel('Time')
plt.ylabel('AS')
plt.grid()

plt.show()
