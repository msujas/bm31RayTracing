import bm31_oasys
import Shadow
import numpy as np
from bm31_oasys import fluxEnergy, fluxDensity, initialPhotons, finalPhotons
import os

energies = [9000,9000, 9000,27000,27000, 27000]
colMirrorAngles = [3.002,3.002,2,3.002, 3.002,2]
torroidalMirrorAngles = [3.002,3.002,3.002,3.002, 3.002, 3.002] #mrad from surface
secondCrystalRots = [0,0.001,0,0,0.001,0]
results = {}
eResults = {}
beams = {}
createdRays = {}
nrays = 1000000
eRange = 50
plot = True
harmonics = [False, False, False, True, True, True]

resultsFile = 'resultsXAS.dat'

for n, (energy, colAngle, torroidalMirrorAngle, secondCrystalRot,harmonic) in enumerate(zip(energies, colMirrorAngles,
                                                                                            torroidalMirrorAngles,secondCrystalRots,harmonics)):
    results[n],eResults[n], beams[n], createdRays[n] = bm31_oasys.run(energy=energy, colMirrorRad=colAngle, eRange=eRange,
                                                      torrAnglemRad=torroidalMirrorAngle, secondCrystalRot=secondCrystalRot,  
                                                      nrays = nrays, writeBeam=True, autoStart=True, harmonic = harmonic)
if os.path.exists(resultsFile):
    os.remove(resultsFile)

for n in results:
    intensity = results[n]['intensity']
    intRatio = intensity/nrays
    energy = energies[n]
    energyIndex = np.abs(fluxEnergy-energy).argmin()
    fluxInitial = fluxDensity[energyIndex]
    print()
    fluxEnd = intRatio*fluxInitial #this is approximating equal flux density in the energy range
    NphotonsI = initialPhotons(fluxInitial,eRange,energy) #approximate

    NphotonsF = finalPhotons(NphotonsI, createdRays[n], intensity)  #approximate
    string = (f"energy: {energy} eV\n"
    f"collimating mirror angle: {colMirrorAngles[n]} mrad\n"
    f"torroidal mirror angle: {torroidalMirrorAngles[n]} mrad\n"
    f"second crystal rotation: {secondCrystalRots[n]} °\n"
    f"harmoic: {harmonics[n]}\n"
    f"source flux density: {fluxInitial:.6e} photons/(s 0.1%bw)\n"
    f"source total photons/s: {NphotonsI:.6e}\n"
    f"intensity: {results[n]['intensity']:.1f}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
    f"final photons/s {NphotonsF:.6e}\n"
    f"fwhm_h: {results[n]['fwhm_h']*10:.6f} mm\n"
    f"fwhm_v: {results[n]['fwhm_v']*10:.6f} mm\n"
    f"energy fwhm: {eResults[n]['fwhm']:.6f} eV\n\n")
    print(string)
    f = open(resultsFile,'a')
    f.write(string)
    f.close()
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title="Real space")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)




