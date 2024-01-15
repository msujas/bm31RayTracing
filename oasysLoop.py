import bm31_oasys
import Shadow
import numpy as np

energies = [9000,9000,15000]
colMirrorAngles = [3.002,3.002,3.002]
torroidalMirrorAngles = [3.002,3.002,3.002] #mrad from surface
secondCrystalRots = [0,0.0005,0]
monoEnergies = [9000,9000,15000]
results = {}
eResults = {}
beams = {}
nrays = 1000000
eRange = 100
plot = True
harmonic = False

fluxFile = r'C:\Users\kenneth1a\Documents\mirrors/spectrum5000_150000.dat'
fluxArray = np.loadtxt(fluxFile,comments='#', unpack=True)
fluxEnergy = fluxArray[0]
fluxDensity = fluxArray[1]
powerDensity = fluxArray[2]

for n, (energy, mEnergy, colAngle, torroidalMirrorAngle, secondCrystalRot) in enumerate(zip(energies, monoEnergies, colMirrorAngles,
                                                                                            torroidalMirrorAngles,secondCrystalRots)):
    results[n],eResults[n], beams[n] = bm31_oasys.run(energy=energy, monoEnergy=mEnergy, colMirrorRad=colAngle, eRange=eRange,
                                                      torrAnglemRad=torroidalMirrorAngle, secondCrystalRot=secondCrystalRot,  
                                                      nrays = nrays, writeBeam=True, autoStart=True, harmonic = harmonic)


for n in results:
    intensity = results[n]['intensity']
    intRatio = intensity/nrays
    energy = energies[n]
    energyIndex = np.abs(fluxEnergy-energy).argmin()
    fluxInitial = fluxDensity[energyIndex]
    print(f'source flux density: {fluxInitial:.6e}')
    fluxEnd = intRatio*fluxInitial #this is approximating equal flux density in the energy range
    NphotonsI = fluxInitial*eRange/(energy/1000) #approximate
    print(f'source total photons: {NphotonsI:.6e}')
    NphotonsF = NphotonsI*intRatio #approximate
    string = (f"intensity: {results[n]['intensity']:.1f}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
    f"final flux: {fluxEnd:.6e}\n"
    f"final photons {NphotonsF:.6e}\n"
    f"fwhm_h: {results[n]['fwhm_h']*10:.6f} mm\n"
    f"fwhm_v: {results[n]['fwhm_v']*10:.6f} mm\n"
    f"energy fwhm: {eResults[n]['fwhm']:.6f} eV\n")
    print(string)
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title="Real space")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)
    print()



