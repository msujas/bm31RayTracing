import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams.update({'font.size':14})

file = 'resultsXASdf_PtCol.dat'
df = pd.read_csv(file,sep = '\t',index_col=0)
hr = False
if 'harmonicRatio' in df.columns:
    hr = True

#df['harmonicN'] = df['harmonic'].map({True:1,False:0})
df2 = df.loc[(df['harmonic'] == False)]
df2Pt = df2.loc[df2['coating1'] == 'Pt']
df2Rh = df2.loc[df2['coating1'] == 'Rh']
df2Si = df2.loc[df2['coating1'] == 'Si']

print(df)
energiesPt = df2Pt['energy(eV)']
energiesSi = df2Si['energy(eV)']
energiesRh = df2Rh['energy(eV)']
intPlotPt = df2Pt['finalPhotons/s']
intPlotRh = df2Rh['finalPhotons/s']
intPlotSi = df2Si['finalPhotons/s']
if hr:
    hrPlotPt = df2Pt['harmonicRatio']
    #hrPlotRh = df2Rh['harmonicRatio']
    #hrPlotSi = df2Si['harmonicRatio']
harmonics = df2['harmonic']
coatings = df2['coating1']
coatingDct = {'Pt':1,'Rh':2,'Si':3}
coatingDctR = {1:'Pt',2:'Rh',3:'Si'}
coatingMap = [coatingDct[coat] for coat in coatings]
fig,ax = plt.subplots(2,1,dpi = 150, figsize = (6.4,9.6))
ax[0].plot(energiesPt,intPlotPt, 'o-', label = 'Pt', color='blue') #plot of (final photons/s)/(fwhmH * fwhmV)
#ax[0].plot(energiesRh,intPlotRh, 'o-', label = 'Rh', color='red') #plot of (final photons/s)/(fwhmH * fwhmV)
#ax[0].plot(energiesSi,intPlotSi, 'o-', label = 'Si', color='green') #plot of (final photons/s)/(fwhmH * fwhmV)
if hr:
    ax[1].plot(energiesPt,hrPlotPt, 'o-', label = 'Pt', color='blue') 
    #ax[1].plot(energiesRh,hrPlotRh, 'o-', label = 'Rh', color='red')
    #ax[1].plot(energiesSi,hrPlotSi, 'o-', label = 'Si', color='green')
ax[0].set_xlabel('energy (eV)')
ax[0].set_ylabel('intensity (photons/s)')
ax[1].set_xlabel('energy (eV)')
ax[1].set_ylabel('harmonic ratio')
ax[1].set_yscale('log')

ax[0].legend()
fig.tight_layout()
#plt.savefig('intensity_energy_mirror.png')
plt.show()
