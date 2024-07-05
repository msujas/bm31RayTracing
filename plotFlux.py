import numpy as np
import matplotlib.pyplot as plt
file = 'spectrum5000_150000.dat'
q = 1.60217663e-19
e,flux,power = np.loadtxt(file,unpack = True)
ekev = e/1000
fig,ax = plt.subplots(dpi=150)
ax.plot(ekev,flux, label = 'flux')
ax2=ax.twinx()
ax2.plot(ekev,flux*e*q,color = 'red',label = 'power density')
ax.set_ylabel(r'flux (photons/(s 0.1% bw)')
ax2.set_ylabel(r'power density (w/(0.1% bw)')
ax.set_xlabel('energy (keV)')
fig.legend(loc= 'upper right', bbox_to_anchor = (0.9,0.9))
ax.set_ylim(0)
ax2.set_ylim(0)
ax.set_xlim(5,150)
plt.tight_layout()
plt.savefig('fluxPlot.png')
plt.show()