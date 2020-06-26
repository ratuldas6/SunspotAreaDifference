import pandas as pd
import numpy as np
from datetime import datetime, timedelta
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates

plt.style.use('seaborn')

data = pd.read_csv('area_1874-2016.csv')

data['Date'] = pd.to_datetime(data['Date'])
data.sort_values('Date', inplace=True)

sunspot_date = data['Date']
sunspot_areadiff = data['Difference']

plt.plot_date(sunspot_date, sunspot_areadiff, markersize = 1, c=(0,0,0))

plt.gcf().autofmt_xdate()

plt.title('Daily Hemispherical Sunspot Area Difference (1874-2016)')
plt.xlabel('Time')
plt.ylabel('Difference in Area')
plt.grid()

plt.show()
