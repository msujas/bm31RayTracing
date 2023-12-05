import bm31_oasys_xrd
import Shadow
import numpy as np

energies = np.linspace(47000,51000,5)
meridionalFs = [1000000]*5
nrays = 500000
monoEnergies = energies
focalEnergy = 47000
results=  {}
eResults = {}
beams = {}
sr = bm31_oasys_xrd.saggitalRadius(focalEnergy)


plot = True

for n,(e,m,me) in enumerate(zip(energies, meridionalFs,monoEnergies)):
    results[n], eResults[n], beams[n] = bm31_oasys_xrd.run(energy=e, monoEnergy=me,nrays= nrays, focalEnergy=focalEnergy,
                                                           meridionalDist = m,  autoStart=True)

for n,(e,m,me) in enumerate(zip(energies, meridionalFs,monoEnergies)):

    intensity = results[n]['intensity']
    fwhmH = results[n]['fwhm_h']*10 #mm
    fwhmV = results[n]['fwhm_v']*10 #mm
    fwhmE = eResults[n]["fwhm"]
    f2 = bm31_oasys_xrd.srTof2(me,sr)
    string = (f"{e} eV\n"
    f"meridional fdist: {m} cm\n"
    f"mono energy: {me} eV\n"
    f"focal distance: {f2:.1f} cm\n"
    f"intensity: {intensity:.1f}\n"
    f"fwhm_h: {fwhmH} mm\n"
    f"fwhm_v: {fwhmV} mm\n"
    f"energy fwhm: {fwhmE} eV\n"
    f"intensity/fwhm_v: {intensity/fwhmV:.1f}\n"
    f"intensity/(fwhm_v*fwhm_h) {intensity/(fwhmV*fwhmH):.1f}\n")
    print(string)
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title=f"xz e={e}, MerFoc={m}")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)#, title=f"energy e={energies[n]}, MerFoc={meridionalFs[n]}")