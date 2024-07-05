import bm31_oasys
import Shadow
import numpy as np
from bm31_oasys import fluxEnergy, fluxDensity, initialPhotons, finalPhotons
import matplotlib.pyplot as plt
import os
import pandas as pd
import time

plot = False
def main(plot = False):
    energiesBase = np.arange(5000,30001,500)
    harmonicsBase = [False]*len(energiesBase)
    harmonicEnergies = energiesBase*3
    harmonicsBase += [True]*len(harmonicEnergies)
    energiesBase = np.append(energiesBase,harmonicEnergies)
    energies = energiesBase
    harmonics = harmonicsBase
    energiesRh = np.arange(5000,30001,500)
    harmonicsRh = [False]*len(energiesRh)
    harmonicEnergiesRh = energiesRh*3
    harmonicsRh += [True]*len(harmonicEnergiesRh)
    energiesRh = np.append(energiesRh,harmonicEnergiesRh)
    energies = energiesBase
    harmonics = harmonicsBase
    #energies = np.append(energies,energiesRh)
    #harmonics = np.append(harmonics,harmonicsRh)
    
    m1Angles = [3]*len(energiesBase)
    m2Angles = [3]*len(energiesBase) #mrad from surface
    #m1Angles += [3]*len(energiesRh)
    #m2Angles += [3]*len(energiesRh)

    
    coating1s = ['Rh']*len(energiesBase)
    #coating1s += ['Rh']*len(energiesRh)
    #coating1s += ['Si']*len(energiesBase)
    coating2s = coating1s

    mirror1types = ['spherical']*len(energies)
    mirror2types = ['spherical']*len(energies)
    #mirror2types += ['spherical']*len(energiesBase)

    torrMajors = [None]*len(energies) #meridional radius. 1667159 - 2mrad, 1067159.1907 - 3mrad, None automatically calculates based on angle
    #torrMajors += [1067159]*len(energiesBase)
    torrMinors = [None]*len(energies)
    #torrMinors += [6.4533]*len(energiesBase)  #saggital radius. 4.24533 - 2mrad, 6.4533 - 3mrad
    secondCrystalRots = [0]*len(energies)

    beams = {}
    nrays = 1000000
    eRange = 50
    slitx = 0.4
    slitz = 0.1
    runMirrors = True

    resultsFile = 'resultsXAS.dat'
    dfFileName = 'resultsXASdf.dat'
    intPlot = []
    if os.path.exists(resultsFile):
        os.remove(resultsFile)
    df = pd.DataFrame(columns = ['energy(eV)', 'harmonic','mirrors', 'collimatingMirrorAngle(mrad)','torroidalMirrorAngle(mrad)','secondCrystalRotation(°)',
                                'coating1','coating2',r'sourceFluxDensity(photons/(s.0.1%bw))', 'sourcePhotons/s','finalPhotons/s', 'finalPhotonDensity(p/(s_mm2))',
                                'totalCreatedRays','intensity','created/accepted','fwhm_h(mm)','fwhm_v(mm)','energyFWHM(eV)'])
    if not len(energies)== len(m1Angles)==  len(m2Angles) == len(secondCrystalRots) ==len(harmonics) ==  len(coating1s) ==  \
    len(coating2s) == len(mirror1types) == len(mirror2types) == len(torrMinors) == len(torrMajors):
        sleeptime = 5
        print(f'input lists are not all equal in length, please check, starting in {sleeptime} s')
        time.sleep(sleeptime)
    print('program:')
    print(energies)
    print(harmonics)
    print(coating1s)
    time.sleep(10)
    for n, (energy, m1angle, m2angle, secondCrystalRot,harmonic, c1, c2, m1, m2,tMin,tMaj) in enumerate(zip(energies, m1Angles,
                                                                                                m2Angles,secondCrystalRots,harmonics,
                                                                                                coating1s, coating2s, mirror1types, mirror2types,
                                                                                                torrMinors,torrMajors)):
        print(f'{n} of {len(energies)}')
        results,eResults, beams[n], createdRays = bm31_oasys.run(energy=energy, colMirrorRad=m1angle, eRange=eRange,
                                                        torrAnglemRad=m2angle, secondCrystalRot=secondCrystalRot,  
                                                        nrays = nrays, writeBeam=True, autoStart=True, harmonic = harmonic, coating1=c1, 
                                                        coating2=c2, torrMajor=tMaj, torrMinor=tMin, slitx=slitx, slitz=slitz,
                                                        mirror1type=m1, mirror2type=m2, runMirrors=runMirrors)
        
        intensity = results['intensity']
        intRatio = intensity/nrays
        try:
            fwhmH = results['fwhm_h']*10
            fwhmHString = f"fwhm_h: {fwhmH:.6f} mm\n"
        except:
            fwhmH = -1
            fwhmHString = 'fwhm_h: NA\n'
        try:
            fwhmV = results['fwhm_v']*10
            
        except:
            fwhmV = -1
            fwhmVString = 'fwhm_v: NA\n'
        fwhmVString = f"fwhm_v: {fwhmV:.6f} mm\n"
        fwhmE = eResults['fwhm']

        energyIndex = np.abs(fluxEnergy-energy).argmin()
        fluxInitial = fluxDensity[energyIndex]
        print()
        fluxEnd = intRatio*fluxInitial #this is approximating equal flux density in the energy range
        NphotonsI = initialPhotons(fluxInitial,eRange,energy) #approximate

        NphotonsF = finalPhotons(NphotonsI, createdRays, intensity)  #approximate
        if fwhmV > 0 and fwhmH > 0:
            photonDensityF = NphotonsF/(fwhmH*fwhmV)
        else:
            photonDensityF = -1
        intPlot.append(photonDensityF)
        if fwhmE == None:
            fwhmEstring = "energy fwhm: NA\n"
            fwhmE = 'NA'
        else:
            fwhmEstring = f"energy fwhm: {fwhmE:.6f} eV\n"


        string = (f"energy: {energy} eV\n"
        f"harmoic: {harmonic}\n"
        f"run mirros: {runMirrors}\n"
        f"collimating mirror angle: {m1angle} mrad\n"
        f"torroidal mirror angle: {m2angle} mrad\n"
        f"second crystal rotation: {secondCrystalRot} °\n"
        f"coating1: {c1}\n"
        f"coating2: {c2}\n"
        f"source flux density: {fluxInitial:.6e} photons/(s 0.1%bw)\n"
        f"source total photons/s: {NphotonsI:.6e}\n"
        f"total created rays: {createdRays}\n"
        f"intensity: {intensity:.2e}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
        f"created/accepted: {createdRays/nrays}\n"
        f"final photons/s: {NphotonsF:.6e}\n"
        f"final photon density (p/(s mm\u00b2)): {photonDensityF:.6e}\n"
        f"{fwhmHString}"
        f"{fwhmVString}"
        f"{fwhmEstring}")
        df.loc[n] = [energy,harmonic,runMirrors,m1angle,m2angle,secondCrystalRot,c1,c2,fluxInitial,NphotonsI,
                    NphotonsF,photonDensityF,createdRays,intensity,createdRays/nrays,fwhmH,fwhmV,fwhmE]
            
        df.to_csv(dfFileName,sep='\t')
        f = open(resultsFile,'a')
        f.write(string)
        f.close()

    df['photonRatio'] = df['finalPhotons/s']/np.max(df['finalPhotons/s'])
    df['harmonicRatio'] = 1
    hrs = df.loc[df['harmonic'] == True]['finalPhotons/s'].values/df.loc[df['harmonic'] == False]['finalPhotons/s'].values
    for index,hr in zip(df.loc[df['harmonic'] == True].index,hrs):
        df.loc[index,'harmonicRatio'] = hr
    for index,hr in zip(df.loc[df['harmonic'] == False].index,hrs):
        df.loc[index,'harmonicRatio'] = hr
    df.to_csv(dfFileName,sep='\t')
    f = open(resultsFile,'r')
    print(f.read())
    f.close()
    print(df)
    for n in beams:
        if plot:
            Shadow.ShadowTools.plotxy(beams[n],1,3,nbins=101,nolost=1,title="Real space")
            Shadow.ShadowTools.histo1(beams[n],11,nbins = 201, nolost=  1, ref = 23)
    
    scatter = plt.scatter(energies,intPlot,c = harmonics,cmap = 'bwr') #plot of (final photons/s)/(fwhmH * fwhmV)
    plt.xlabel('energy (eV)')
    plt.ylabel('photon density (photons/(s mm$^2$)')
    plt.legend(*scatter.legend_elements(),title = 'harmonic')
    plt.show()
    
    return df
if __name__ == '__main__':
    main(plot)