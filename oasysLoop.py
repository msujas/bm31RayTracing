import bm31_oasys
import Shadow
import numpy as np
from bm31_oasys import fluxEnergy, fluxDensity, initialPhotons, finalPhotons
import matplotlib.pyplot as plt
import os
import pandas as pd

energies = np.arange(5000,6001,1000)
harmonics = [False]*len(energies)

harmonicEnergies = energies*3
harmonics += [True]*len(harmonicEnergies)
energies = np.append(energies,harmonicEnergies)


colMirrorAngles = [3.002]*len(energies)
torroidalMirrorAngles = [3.002]*len(energies) #mrad from surface
secondCrystalRots = [0]*len(energies)
coating1s = ['Pt']*len(energies)
coating2s = coating1s
mirror1types = ['spherical']*len(energies)
mirror2types = ['torroidal']*len(energies)
results = {}
eResults = {}
beams = {}
createdRays = {}
nrays = 1000000
eRange = 50
plot = True

resultsFile = 'resultsXAS.dat'

for n, (energy, colAngle, torroidalMirrorAngle, secondCrystalRot,harmonic, c1, c2, m1, m2) in enumerate(zip(energies, colMirrorAngles,
                                                                                            torroidalMirrorAngles,secondCrystalRots,harmonics,
                                                                                            coating1s, coating2s, mirror1types, mirror2types)):
    results[n],eResults[n], beams[n], createdRays[n] = bm31_oasys.run(energy=energy, colMirrorRad=colAngle, eRange=eRange,
                                                      torrAnglemRad=torroidalMirrorAngle, secondCrystalRot=secondCrystalRot,  
                                                      nrays = nrays, writeBeam=True, autoStart=True, harmonic = harmonic, coating1=c1, 
                                                      coating2=c2)

intPlot = []
if os.path.exists(resultsFile):
    os.remove(resultsFile)

df = pd.DataFrame(columns = ['energy(eV)', 'harmonic', 'collimatingMirrorAngle(mrad)','torroidalMirrorAngle(mrad)','secondCrystalRotation(°)',
                             'coating1','coating2',r'sourceFluxDensity(photons/(s.0.1%bw))', 'sourcePhotons/s','finalPhotons/s',
                             'totalCreatedRays','intensity','created/accepted','fwhm_h(mm)','fwhm_v(mm)','energyFWHM(eV)'])
for n in results:
    intensity = results[n]['intensity']
    
    intRatio = intensity/nrays
    energy = energies[n]
    harmonic = harmonics[n]
    fwhmH = results[n]['fwhm_h']*10
    fwhmV = results[n]['fwhm_v']*10
    fwhmE = eResults[n]['fwhm']
    colAngle = colMirrorAngles[n]
    torrAngle = torroidalMirrorAngles[n]
    secCrystRot = secondCrystalRots[n]
    energyIndex = np.abs(fluxEnergy-energy).argmin()
    fluxInitial = fluxDensity[energyIndex]
    coat1 = coating1s[n]
    coat2 = coating2s[n]
    cr = createdRays[n]
    print()
    fluxEnd = intRatio*fluxInitial #this is approximating equal flux density in the energy range
    NphotonsI = initialPhotons(fluxInitial,eRange,energy) #approximate

    NphotonsF = finalPhotons(NphotonsI, createdRays[n], intensity)  #approximate
    photonDensityF = NphotonsF/(fwhmH*fwhmV)
    intPlot.append(photonDensityF)
    if fwhmE == None:
        fwhmEstring = f"energy fwhm: NA\n"
        fwhmE = 'NA'
    else:
        fwhmEstring = f"energy fwhm: {fwhmE:.6f} eV\n"
    string = (f"energy: {energy} eV\n"
    f"harmoic: {harmonic}\n"
    f"collimating mirror angle: {colAngle} mrad\n"
    f"torroidal mirror angle: {torrAngle} mrad\n"
    f"second crystal rotation: {secCrystRot} °\n"
    f"coating1: {coat1}\n"
    f"coating2: {coat2}\n"
    f"source flux density: {fluxInitial:.6e} photons/(s 0.1%bw)\n"
    f"source total photons/s: {NphotonsI:.6e}\n"
    f"total created rays: {cr}\n"
    f"intensity: {intensity:.1e}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
    f"created/accepted: {cr/nrays}\n"
    f"final photons/s: {NphotonsF:.6e}\n"
    f"fwhm_h: {fwhmH:.6f} mm\n"
    f"fwhm_v: {fwhmV:.6f} mm\n"
    f"{fwhmEstring}")
    print(string)
    df.loc[n] = [energy,harmonic,colAngle,torrAngle,secCrystRot,coat1,coat2,fluxInitial,NphotonsI,
                 NphotonsF,cr,intensity,cr/nrays,fwhmH,fwhmV,fwhmE]

    f = open(resultsFile,'a')
    f.write(string)
    f.close()
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title="Real space")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)

df.to_csv('resultsXASdf.dat',sep='\t')
hcolourdct = {True:'blue',False:'red'}
hcolours = [hcolourdct[x] for x in harmonics]
plt.scatter(energies,intPlot,c = harmonics,cmap = 'bwr') #plot of (final photons/s)/(fwhmH * fwhmV)
plt.xlabel('energy (eV)')
plt.ylabel('photon density (photons/(s mm$^2$)')
plt.legend(title = 'harmonic')
plt.show()
