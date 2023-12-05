import bm31_oasys
import Shadow

energies = [9000,9000,15000]
colMirrorAngles = [3.002,3.002,3.002]
torroidalMirrorAngles = [3.002,3.002,3.002] #mrad from surface
secondCrystalRots = [0,0.0005,0]
monoEnergies = [9000,9000,15000]
results = {}
eResults = {}
beams = {}
nrays = 1000000

plot = True

for n, (energy, mEnergy, colAngle, torroidalMirrorAngle, secondCrystalRot) in enumerate(zip(energies, monoEnergies, colMirrorAngles,
                                                                                            torroidalMirrorAngles,secondCrystalRots)):
    results[n],eResults[n], beams[n] = bm31_oasys.run(energy=energy, monoEnergy=mEnergy, colMirrorRad=colAngle, 
                                                      torrAnglemRad=torroidalMirrorAngle, secondCrystalRot=secondCrystalRot,  
                                                      nrays = nrays, writeBeam=True, autoStart=True)


for n in results:
    string = (f"intensity: {results[n]['intensity']:.1f}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
    f"fwhm_h: {results[n]['fwhm_h']*10:.6f} mm\n"
    f"fwhm_v: {results[n]['fwhm_v']*10:.6f} mm\n"
    f"energy fwhm: {eResults[n]['fwhm']:.6f} eV\n")
    print(string)
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title="Real space")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)



