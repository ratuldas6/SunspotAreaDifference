import pandas as pd
from matplotlib import pyplot as plt
from matplotlib import dates as mpl_dates

def makeBinsN(startpt,sample_range,N):
    bins = []
    size = sample_range/N
    for i in range(N+1):
        x = startpt + i*size
        bins.append(x)
    return bins

def unifyBins(bins):
    classmarks = []
    for i in range(len(bins)-1):
        cm = (bins[i]+bins[i+1])/2
        classmarks.append(cm)
    return classmarks

def binFreq(AS,bins):
    freqs = [0 for x in range(len(bins)-1)]
    for i in range(len(bins)-1):
        for s in AS:
            if s >= bins[i] and s <= bins[i+1]:
                freqs[i]+=1
    return freqs

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

bins = makeBinsN(-1,2,25)
binuni = unifyBins(bins)
binf = binFreq(AS,bins)
lowerline = [0 for x in range(-2,2)]

fig=plt.figure()
ax=fig.add_subplot()
ax.set_title('Line Plot for AS (Monthly, 1874-2016)')
ax.set_xlabel('AS')
ax.set_ylabel('Frequency')
ax.plot(binuni,binf)

plt.show()
