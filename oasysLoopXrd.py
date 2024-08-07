import bm31_oasys_xrd
import Shadow
import numpy as np
from bm31_oasys import fluxEnergy, fluxDensity, initialPhotons, finalPhotons
import os
import pandas as pd

energies =  [47500]*3 +[47500*3] #np.linspace(47000,51000,5)
meridionalFs = [2000,5000,100000, 100000]
harmonics = [False,False,False, True]
nrays = 500000
focalEnergy = 47000
results=  {}
eResults = {}
beams = {}
createdRays = {}
sr = bm31_oasys_xrd.saggitalRadius(focalEnergy)
eRange = 200

plot = False

resultsFile = 'resultsXRD.dat'

for n,(e,m,h) in enumerate(zip(energies, meridionalFs,harmonics)):
    results[n], eResults[n], beams[n], createdRays[n] = bm31_oasys_xrd.run(energy=e, nrays= nrays, focalEnergy=focalEnergy, eRange=eRange,
                                                           meridionalDist = m,  autoStart=True, imageDist=bm31_oasys_xrd.f2 + 20, harmonic = h)

if os.path.exists(resultsFile):
    os.remove(resultsFile)
df = pd.DataFrame(columns = ['energy(eV)', 'harmonic', 'focalDist(cm)','merFocalDist(cm)', r'sourceFluxDensity(photons/(s.0.1%bw))', 
                             'sourcePhotons/s','finalPhotons/s','finalPhotonDensity(p/(s_mm2))','totalCreatedRays','intensity','created/accepted','fwhm_h(mm)','fwhm_v(mm)',
                             'energyFWHM(eV)'])
for n in results:
    e = energies[n]
    m = meridionalFs[n]
    cr = createdRays[n]
    intensity = results[n]['intensity']
    fwhmH = results[n]['fwhm_h']*10 #mm
    fwhmV = results[n]['fwhm_v']*10 #mm
    fwhmE = eResults[n]["fwhm"]
    energyIndex = np.abs(fluxEnergy-e).argmin()
    fluxInitial = fluxDensity[energyIndex]
    NphotonsI = initialPhotons(fluxInitial,eRange,e) #approximate
    NphotonsF = finalPhotons(NphotonsI, cr, intensity) #approximate
    photonDensityF = NphotonsF/(fwhmH*fwhmV)
    print()
    f2 = bm31_oasys_xrd.srTof2(e,harmonics[n],sr)
    string = (f"{e} eV\n"
    f"meridional fdist: {m} cm\n"
    f"focal distance: {f2:.1f} cm\n"
    f"harmonic: {harmonics[n]}\n"
    f"source flux density: {fluxInitial:.6e} photons/(s 0.1%bw)\n"
    f"source total photons/s: {NphotonsI:.6e}\n"
    f"total created rays: {cr}\n"
    f"created/accepted: {cr/nrays}\n"
    f"intensity: {intensity:.1e}\n"
    f"final photons/s: {NphotonsF:.6e}\n"
    f"final photon denstiy (p/(s mm\u00b2)): {photonDensityF:.6e}\n"
    f"fwhm_h: {fwhmH} mm\n"
    f"fwhm_v: {fwhmV} mm\n"
    f"energy fwhm: {fwhmE} eV\n"
    f"intensity/fwhm_v: {intensity/fwhmV:.1f}\n"
    f"intensity/(fwhm_v*fwhm_h): {intensity/(fwhmV*fwhmH):.1f}\n\n")
    print(string)
    f = open(resultsFile,'a')
    f.write(string)
    f.close()
    df.loc[n] = [e,harmonics[n],f2,m,fluxInitial,NphotonsI,NphotonsF,photonDensityF,cr,intensity,cr/nrays,fwhmH,fwhmV,fwhmE]
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title=f"xz e={e}, MerFoc={m}")
        if harmonics[n]:
            nbins = 10001
        else:
            nbins = 201
        Shadow.ShadowTools.histo1(beams[n],11,nbins = nbins, nolost=  1, ref = 23)#, title=f"energy e={energies[n]}, MerFoc={meridionalFs[n]}")


df.to_csv('resultsXRDdf.dat',sep='\t')