import bm31_oasys_xrd
import Shadow

energies = [49000,49000, 49000,49000]
meridionalFs = [1811.6,3000,10000,1000000]
nrays = 500000

results=  {}
eResults = {}
beams = {}


f1 = 3032.8
f2 = 1811.6
plot = True

for n,(e,m) in enumerate(zip(energies, meridionalFs)):
    results[n], eResults[n], beams[n] = bm31_oasys_xrd.run(energy=e, monoEnergy=e,nrays= nrays, meridionalDist = m, autoStart=True)

for n,(e,m) in enumerate(zip(energies, meridionalFs)):
    print(f'{e} eV')
    print(f'meridional fdist: {m}')
    intensity = results[n]['intensity']
    fwhmH = results[n]['fwhm_h']*10 #mm
    fwhmV = results[n]['fwhm_v']*10 #mm
    fwhmE = eResults[n]["fwhm"]
    print(f'intensity: {intensity:.1f}')
    print('fwhm_h:',fwhmH, 'mm')
    print('fwhm_v:',fwhmV, 'mm')
    print(f'energy fwhm: {fwhmE} eV')
    print(f'intensity/fwhm_v: {intensity/fwhmV:.1f}')
    print(f'intensity/(fwhm_v*fwhm_h) {intensity/(fwhmV*fwhmH):.1f}')
    if plot:
        Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title=f"xz e={e}, MerFoc={m}")
        Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)#, title=f"energy e={energies[n]}, MerFoc={meridionalFs[n]}")