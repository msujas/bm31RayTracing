#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
np = numpy
from srxraylib.sources import srfunc
import matplotlib.pyplot as plt
from bm31_oasys import dctToFile, readConfig, whereStart, fluxEnergy, fluxDensity, direc, readCreatedRays, writeCreatedRays, initialPhotons, finalPhotons
import os

energy = 47000*3
focalEnergy = 47000
meridionalDist = 10000
nrays = 500000
eRange = 200
harmonic = True
autoStart = True #set to False for first use


f1 = 3032.8
f2 = 1811.6

def mradSurface_to_degNorm(mrad):
    return 90-mrad*180/(numpy.pi*1000)

def energyToTheta(energy):
    dspacing = 3.13379
    wavelength = 12398.47/energy
    return np.arcsin(wavelength/(2*dspacing))*180/np.pi

def saggitalRadius(energy,f1=f1,f2=f2):
    theta = energyToTheta(energy)*np.pi/180
    return (2*f1*f2*np.sin(theta)/(f1+f2))

sr = saggitalRadius(focalEnergy,f1,f2)

def srTof2(energy, harmonic,sr=sr,f1=f1):
    if harmonic:
        energy /=3
    theta = energyToTheta(energy)*np.pi/180
    return f1*sr/(2*f1*np.sin(theta) - sr)

def meridionalRadius(f1,f2,focalEnergy):
    theta = energyToTheta(focalEnergy)*np.pi/180
    return ((f1+f2)**2 - 4*f1*f2*np.sin(theta)**2)**0.5/(np.sin(2*theta))

def run(energy = 49000,  focalEnergy = 49000, meridionalDist = 1000000, writeBeam = True, 
        nrays = 1000000, traceStart = 0, autoStart = False, eRange = eRange,imageDist = f2, harmonic = False):
    
    configFile = 'config/xrdConfig.dat'
    createdRaysLog = 'config/createdRaysXRD.log'

    file333 = bytes(f'{direc}/bragg333.dat',encoding = 'utf-8')
    file111 = bytes(f'{direc}/Si5_55.111',encoding = 'utf-8')

    if harmonic:
        braggFile = file333
    else:
        braggFile = file111
    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    beamFile = 'config/starXRD'
    config = {'energy':energy,  'nrays':nrays, 'focalEnergy': focalEnergy, 'meridionalDist': meridionalDist,
              'f2':f2,'eRange':eRange, 'harmonic':harmonic}
    startDct = {'energy':0,'nrays':0,'focalEnergy':2, 'meridionalDist': 2, 'eRange':0,'harmonic':1}
    if os.path.exists(configFile):
        oldConfig = readConfig(configFile)
        autoTraceStart = whereStart(config,oldConfig,startDct)
        if autoStart:
            traceStart = autoTraceStart
    if traceStart > 0:
        loadfile = f'{beamFile}.{traceStart-1:02d}'
        print(f'loading {loadfile}')
        beam.load(loadfile)

    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    #generating source trajectory files
    traj, pars = srfunc.wiggler_trajectory(
    b_from            =1,
    inData            =f"{direc}/SW_2PA.txt",
    nPer              =1,
    nTrajPoints       =501,
    ener_gev          =6.0,
    per               =0.2,
    kValue            =7.85,
    trajFile          ="tmp.traj",
    shift_x_flag      =4,
    shift_x_value     =0.0,
    shift_betax_flag  =5,
    shift_betax_value =0.0055)
    srfunc.wiggler_cdf(traj,
    enerMin        = energy-eRange/2,
    enerMax        = energy+eRange/2,
    enerPoints     = 2001,
    outFile        =bytes(f'{direc}/xshwig.sha',encoding = 'utf-8'),
    elliptical     =False)

    #oe0 - wiggler source
    oe0.BENER = 6.0
    oe0.CONV_FACT = 100.0
    oe0.EPSI_DX = 89.4
    oe0.EPSI_DZ = -104.8
    oe0.EPSI_X = 2.16e-08
    oe0.EPSI_Z = 5e-10
    oe0.FDISTR = 0
    oe0.FILE_BOUND = bytes(f'{direc}/myslit.dat',encoding = 'utf-8')
    oe0.FILE_TRAJ = bytes(f'{direc}/xshwig.sha', encoding = 'utf-8')
    oe0.FSOUR = 0
    oe0.FSOURCE_DEPTH = 0
    oe0.F_BOUND_SOUR = 2
    oe0.F_COLOR = 0
    oe0.F_PHOT = 0
    oe0.F_WIGGLER = 1
    oe0.HDIV1 = 1.0
    oe0.HDIV2 = 1.0
    oe0.IDO_VX = 0
    oe0.IDO_VZ = 0
    oe0.IDO_X_S = 0
    oe0.IDO_Y_S = 0
    oe0.IDO_Z_S = 0
    oe0.ISTAR1 = 5676561 #seed value, takes any odd number from 1000 to 1000000
    oe0.NCOL = 0
    oe0.NPOINT = nrays
    oe0.NTOTALPOINT = 0
    oe0.N_COLOR = 0
    oe0.PH1 = energy - eRange/2
    oe0.PH2 = energy + eRange/2
    oe0.POL_DEG = 0.0
    oe0.SIGMAX = 0.0008757
    oe0.SIGMAY = 0.0
    oe0.SIGMAZ = 0.0001647
    oe0.VDIV1 = 1.0
    oe0.VDIV2 = 1.0
    oe0.WXSOU = 0.0
    oe0.WYSOU = 0.0
    oe0.WZSOU = 0.0

    print(oe0.PH1)
    #first crystal
    oe1.DUMMY = 1.0
    oe1.FHIT_C = 1
    oe1.FILE_REFL = braggFile
    oe1.FWRITE = 1
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_MOVE = 1
    oe1.PHOT_CENT = energy
    oe1.RLEN1 = 2.5
    oe1.RLEN2 = 2.5
    oe1.RWIDX1 = 3.5
    oe1.RWIDX2 = 3.5
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 5.0
    oe1.T_SOURCE = f1

    #second crystal
    oe2.ALPHA = 180.0
    oe2.DUMMY = 1.0
    oe2.FHIT_C = 1
    oe2.FILE_REFL = braggFile
    oe2.FMIRR = 3
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.F_EXT = 1
    oe2.F_MOVE = 1
    oe2.F_TORUS = 2
    oe2.PHOT_CENT = energy
    oe2.RLEN1 = 2.5
    oe2.RLEN2 = 2.5
    oe2.RWIDX1 = 3.5
    oe2.RWIDX2 = 3.5
    oe2.R_LAMBDA = 5000.0
    oe2.R_MAJ = meridionalRadius(f1,meridionalDist,focalEnergy)
    oe2.R_MIN = saggitalRadius(focalEnergy,f1,f2)
    oe2.T_IMAGE = imageDist
    oe2.T_SOURCE = 5.0
    
    #oe3.DUMMY = 1.0
    #oe3.FWRITE = 3
    #oe3.F_REFRAC = 2
    #oe3.F_SCREEN = 1
    #oe3.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    #oe3.N_SCREEN = 1
    #oe3.RX_SLIT = numpy.array([2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #oe3.RZ_SLIT = numpy.array([0.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    #oe3.T_IMAGE = 0.0
    #oe3.T_INCIDENCE = 0.0
    #oe3.T_REFLECTION = 180.0
    #oe3.T_SOURCE = 0.0
    

    if traceStart < 1:
        beam.genSource(oe0)
        createdRays = oe0.NTOTALPOINT
        writeCreatedRays(createdRays, createdRaysLog)

    if writeBeam and traceStart < 1:
        beam.write(f"{beamFile}.00")

    #run optical element 1
    if traceStart < 2:
        print("    Running optical element: %d"%(1))
        beam.traceOE(oe1,1)

    if writeBeam and traceStart < 2:
        beam.write(f"{beamFile}.01")


    #
    #run optical element 2
    #
    if traceStart < 3:
        print("    Running optical element: %d"%(2))
        beam.traceOE(oe2,2)

    if writeBeam:
        beam.write(f"{beamFile}.02")

    #
    #run optical element 3
    #
    #print("    Running optical element: %d"%(3))

    
    #beam.traceOE(oe3,3)
    #if writeBeam:
    #    beam.write(f"{beamFile}.03")
    result = beam.histo2(1,3, nbins= 101,nolost=1)
    if harmonic:
        nbins = 10001
    else:
        nbins = 201
    eResult = beam.histo1(11,nbins = nbins, nolost=  1, ref = 23)
    #Shadow.ShadowTools.plotxy(beam,1,3,block=block,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    #plt.imshow(beam.rays, aspect = 'auto')
    #plt.show()
    dctToFile(config,configFile)
    createdRays = readCreatedRays(createdRaysLog)
    return result, eResult,beam, createdRays

if __name__ == '__main__':
    result, eResult, beam, createdRays = run(energy = energy, nrays = nrays, focalEnergy=focalEnergy, 
                                meridionalDist = meridionalDist, autoStart=autoStart, harmonic=harmonic)
    intensity = result['intensity']
    fwhmH = result['fwhm_h']*10 #mm
    fwhmV = result['fwhm_v']*10 #mm
    fwhmE = eResult["fwhm"]
    intRatio = intensity/nrays
    energyIndex = np.abs(fluxEnergy-energy).argmin()
    fluxInitial = fluxDensity[energyIndex]
    NphotonsI = initialPhotons(fluxInitial,eRange,energy) #approximate
    NphotonsF = finalPhotons(NphotonsI,createdRays,intensity) #approximate
    print()
    fluxEnd = intRatio*fluxInitial #this is approximating equal flux density in the energy range
    f2 = srTof2(energy, harmonic,sr)
    string = (f"{energy} eV\n"
    f"source flux density: {fluxInitial:.6e}\n photons/(s 0.1%bw)"
    f"source total photons/s: {NphotonsI:.6e}\n"
    f"total created rays: {createdRays}\n"
    f"created/accepted rays = {createdRays/nrays:.6f}\n"
    f"meridional fdist: {meridionalDist} cm\n"
    f"mono energy: {energy} eV\n"
    f"focal distance: {f2:.1f} cm\n"
    f"intensity: {intensity:.1e}\n"
    f"final photons/s: {NphotonsF:.6e}\n"
    f"fwhm_h: {fwhmH} mm\n"
    f"fwhm_v: {fwhmV} mm\n"
    f"energy fwhm: {fwhmE} eV\n"
    f"intensity/fwhm_v: {intensity/fwhmV:.1f}\n"
    f"intensity/(fwhm_v*fwhm_h) {intensity/(fwhmV*fwhmH):.1f}\n")
    print(string)
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    if harmonic:
        nbins = 10001
    else:
        nbins = 501
    Shadow.ShadowTools.histo1(beam,11,nbins = nbins, nolost=  1, ref = 23)
    
    #plt.stairs(eResult['histogram'],eResult['bins'])
    #plt.xlabel('energy (eV)')
    #plt.ylabel('intensity')
    #plt.show()
    #print(eResult)
    #print(result)
    #plt.imshow(result['histogram'].transpose(), extent = [result['bin_h_edges'][0],result['bin_h_edges'][-1],
    #                                                      result['bin_v_edges'][0],result['bin_v_edges'][-1]], aspect = 'auto')
    #plt.show()