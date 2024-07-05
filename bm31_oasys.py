#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
np = numpy
from srxraylib.sources import srfunc
import os
import matplotlib.pyplot as plt
from enum import Enum

direc = os.path.dirname(os.path.realpath(__file__))
energy = 9000
harmonic = False

secondCrystalRot = 0.00 #0.001 for detuning 60%
firstMirrorAngle = 3
torroidalMirrorAngle = 3 #mrad from surface
torrMajor = 1067159.
torrMinor = 6.4533
coating1 = 'Si'
coating2 = 'Si'
mirror1type = 'spherical'
mirror2type = 'spherical'
slitx = 1
slitz = 1
runMirrors = True
dspacing = 3.13379
nrays = 1000000
eRange = 50

autoStart = True #set to False with first use
traceStart = 0


writeBeam = True

fluxFile = rf'{direc}/spectrum5000_150000.dat'
fluxArray = np.loadtxt(fluxFile,comments='#', unpack=True)
fluxEnergy = fluxArray[0]
fluxDensity = fluxArray[1]
powerDensity = fluxArray[2]

def readCreatedRays(createdRaysLog):
    f = open(createdRaysLog,'r')
    createdRays = int(f.read())
    f.close()
    return createdRays

def writeCreatedRays(createdRays, createdRaysLog):
    f = open(createdRaysLog,'w')
    f.write(str(createdRays))
    f.close()

def mradSurface_to_degNorm(mrad):
    return 90-mrad*180/(numpy.pi*1000)


def energy_to_incidentAngle(energy):
    wavelength = 12398.47/energy
    return 90-np.arcsin(wavelength/(2*dspacing))*180/np.pi

def dctToFile(dct,file):
    string = ''
    for key in dct:
        string += f'{key};{dct[key]}\n'
    f = open(file,'w')
    f.write(string)
    f.close()

def readConfig(file):
    dct = {}
    f = open(file,'r')
    lines = f.read().split('\n')
    f.close()
    for line in lines[:-1]:
        lsplit = line.split(';')
        dct[lsplit[0]] = lsplit[1]
    return dct

def initialPhotons(flux, eRange, energy):
    return flux*eRange/(energy/1000)

def finalPhotons(initialPhotons, totalnrays, intensity):
    return initialPhotons*intensity/totalnrays

def whereStart(config,oldConfig, startDct):
    startDctSorted = dict(sorted(startDct.items(), key = lambda x: x[1]))
    for item in startDctSorted:
        if item not in list(oldConfig.keys()):
            return startDct[item]
        elif str(config[item]) != oldConfig[item]:
            return startDct[item]
    return startDct[item]+1

def mirrorSagitalRadius(angle,f1 = 3242.7,f2 = 1601.8):
    #angle in deg
    angleRad = angle*np.pi/180
    return (2*f1*f2*np.sin(angleRad)/(f1+f2))

def mirrorMerRadius(angle, f1 = 10**26, f2= 1601.8):
    angleRad = angle*np.pi/180
    sRad = mirrorSagitalRadius(angleRad,f1,f2)
    return sRad/np.sin(angleRad)**2

class Coating(Enum):
    Rh = 'Rh'
    Pt = 'Pt'
    Si = 'Si'
    def file(self):
        if self.value == 'Rh':
            return bytes(f'{direc}/RhReflec_4p9_150.dat',encoding = 'utf-8')
        elif self.value == 'Pt':
            return bytes(f'{direc}/PtReflec_4p9_150.dat',encoding = 'utf-8')
        elif self.value == 'Si':
            return bytes(f'{direc}/SiReflec_4p9_150.dat',encoding = 'utf-8')

def mirrorType(coatingFile,angle, mtype, sdist, imgdist,colsimg, colssour,torrMajor = None, torrMinor = None):
    #dist = 200 for mirror 1, 209.9 for mirror2
    if torrMinor == None:
        torrMinor = mirrorSagitalRadius(90-angle)
    if torrMajor == None:
        torrMajor = mirrorMerRadius(90-angle)
    oe = Shadow.OE()
    if mtype == 'collimating' or mtype == 'spherical':
        oe.DUMMY = 1.0
        oe.FCYL = 1
        oe.FHIT_C = 1
        oe.FILE_REFL = coatingFile
        oe.FILE_RIP = bytes(f'{direc}/Colmirror1m_1u.dat', encoding = 'utf-8')
        oe.FMIRR = 1
        oe.F_DEFAULT = 0
        oe.F_G_S = 2
        oe.F_REFLEC = 1
        oe.F_RIPPLE = 1
        oe.RLEN1 = 50.0
        oe.RLEN2 = 50.0
        oe.RWIDX1 = 4.0
        oe.RWIDX2 = 4.0
        oe.SIMAG = colsimg #1.00000003e+16 #image side focal diastance
        oe.SSOUR = colssour #2886.3 #object side focal distance
        oe.THETA = angle
        oe.T_IMAGE = imgdist
        oe.T_INCIDENCE = angle
        oe.T_REFLECTION = angle
        oe.T_SOURCE = sdist
    elif mtype == 'torroidal':
        oe.DUMMY = 1.0
        oe.FHIT_C = 1
        oe.FILE_REFL = coatingFile
        oe.FILE_RIP = bytes(f'{direc}/Colmirror1m_1u.dat',encoding = 'utf-8')
        oe.FMIRR = 3
        oe.FWRITE = 2
        oe.F_ANGLE = 1
        oe.F_EXT = 1
        oe.F_G_S = 2
        oe.F_REFLEC = 1
        oe.F_RIPPLE = 1
        oe.RLEN1 = 50.0
        oe.RLEN2 = 50.0
        oe.RWIDX1 = 5.0
        oe.RWIDX2 = 5.0
        oe.R_MAJ = torrMajor #meridional radius. 1667159 - 2mrad, 1067159.1907 - 3mrad
        oe.R_MIN = torrMinor #saggital radius. 4.24533 - 2mrad, 6.4533 - 3mrad
        oe.T_IMAGE = imgdist
        oe.T_INCIDENCE = angle
        oe.T_REFLECTION = angle
        oe.T_SOURCE = sdist
    else:
        raise ValueError(f'"{mtype}" is an invalid mirror type. Must be "spherical"/"collimating" or "torroidal"')
    return oe

if not os.path.exists(f'{direc}/config/'):
    os.makedirs(f'{direc}/config/')

def run(energy = 9000, eRange = 100, colMirrorRad = 3, torrAnglemRad = 3, secondCrystalRot = 0, writeBeam=True, 
        nrays = 100000, traceStart = 0, autoStart = False, harmonic = False, coating1 = 'Rh', coating2 = 'Rh',mirror1type = 'spherical',
        mirror2type = 'torroidal', torrMajor = 1067159.1907, torrMinor = 6.4533, slitx = 0.4, slitz = 1, runMirrors=True):

    torrAngleDeg = mradSurface_to_degNorm(torrAnglemRad)
    colMirrorDeg = mradSurface_to_degNorm(colMirrorRad)

    configFile = 'config/mirrorsConfig.dat'
    createdRaysLog = 'config/createdRays.log'
    fname = 'output.dat'

    config = {'energy':energy,'eRange':eRange, 'colMirrorRad': colMirrorRad, 'torrAnglemRad':torrAnglemRad, 'secondCrystalRot':secondCrystalRot,
              'nrays':nrays, 'harmonic':harmonic, 'coating1':coating1, 'coating2':coating2,'torrMajor':torrMajor,'torrMinor':torrMinor, 
              'slitx':slitx, 'slitz':slitz, 'runMirrors':runMirrors}
    startDct = {'energy':0, 'eRange': 0, 'colMirrorRad':2, 'secondCrystalRot':4,'torrAnglemRad':5,'nrays':0, 'harmonic': 3, 'coating1':2,
                'coating2':5, 'torrMajor':5, 'torrMinor':5, 'slitx':6,'slitz':6,'runMirrors':2}
    if os.path.exists(configFile):
        oldConfig = readConfig(configFile)
        autoTraceStart = whereStart(config,oldConfig,startDct)
        if autoStart:
            traceStart = autoTraceStart
    #
    # initialize shadow3 source (oe0) and beam
    #

    file111 = bytes(f'{direc}/bragg111.dat',encoding = 'utf-8')
    file333 = bytes(f'{direc}/bragg333.dat',encoding = 'utf-8')
    if harmonic:
        dcmfile = file333
    else:
        dcmfile = file111
    
    coatingFile1 = Coating(coating1).file()
    coatingFile2 = Coating(coating2).file()
    
    print(f'beginning at {traceStart}')
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    #oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    oe4 = Shadow.OE()
    #oe5 = Shadow.OE()
    oe6 = Shadow.OE()
    beam = Shadow.Beam()




    beamFile = 'config/beam'
    if traceStart > 0:
        beam.load(f'{beamFile}.{traceStart-1:02d}')

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
    enerPoints     = 1001,
    outFile        =bytes(f'{direc}/xshwig.sha',encoding = 'utf-8'),
    elliptical     =False)

    calculate_spectrum = False

    if calculate_spectrum:
        e, f, w = srfunc.wiggler_spectrum(traj,
            enerMin=5000.0,
            enerMax=150000.0,
            nPoints=10000,
            electronCurrent=200*1e-3,
            outFile=fluxFile,
            elliptical=False)
        from srxraylib.plot.gol import plot
        plot(e, f, xlog=False, ylog=False,show=False,
            xtitle="Photon energy [eV]",ytitle="Flux [Photons/s/0.1%bw]",title="Flux")
        plot(e, w, xlog=False, ylog=False,show=True,
            xtitle="Photon energy [eV]",ytitle="Spectral Power [E/eV]",title="Spectral Power") 
        return        
    #oe0 - wiggler source
    oe0.BENER = 6.0
    oe0.CONV_FACT = 100.0
    oe0.EPSI_DX = 89.4
    oe0.EPSI_DZ = -104.8
    oe0.EPSI_X = 2.16e-08
    oe0.EPSI_Z = 5e-10
    oe0.FDISTR = 0
    oe0.FILE_BOUND = bytes(f'{direc}\\myslit.dat',encoding = 'utf-8')
    oe0.FILE_TRAJ = bytes(f'{direc}\\xshwig.sha',encoding = 'utf-8')
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

    #oe1 - pslits
    oe1.DUMMY = 1.0
    oe1.FWRITE = 3
    oe1.F_REFRAC = 2
    oe1.F_SCREEN = 1
    oe1.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe1.N_SCREEN = 1
    oe1.RX_SLIT = numpy.array([2.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.RZ_SLIT = numpy.array([0.3, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe1.T_IMAGE = 0.0
    oe1.T_INCIDENCE = 0.0
    oe1.T_REFLECTION = 180.0
    oe1.T_SOURCE = 2686.3

    #oe2 - first mirror (collimatring)
    oe2 = mirrorType(coatingFile1,colMirrorDeg, sdist = 200, imgdist=0, colsimg=1e16, colssour=2886.3,mtype = mirror1type)

    #oe3 - first mono crystal
    oe3.ALPHA = 0
    oe3.DUMMY = 1.0
    oe3.FHIT_C = 1
    oe3.FILE_REFL = dcmfile
    oe3.FWRITE = 1
    oe3.F_CENTRAL = 1
    oe3.F_CRYSTAL = 1
    oe3.F_MOVE = 1
    oe3.PHOT_CENT = energy
    oe3.RLEN1 = 2.5
    oe3.RLEN2 = 2.5
    oe3.RWIDX1 = 3.5
    oe3.RWIDX2 = 3.5
    oe3.R_LAMBDA = 5000.0
    oe3.FSHAPE = 1
    oe3.T_IMAGE = 0.0
    oe3.T_INCIDENCE = 0
    oe3.T_REFLECTION = 0
    oe3.T_SOURCE = 146.4

    #oe4 - second mono crystal
    oe4.ALPHA = 180.0
    oe4.DUMMY = 1.0
    oe4.FHIT_C = 1
    oe4.FILE_REFL = dcmfile
    oe4.FWRITE = 1
    oe4.F_CENTRAL = 1
    oe4.F_CRYSTAL = 1
    oe4.F_MOVE = 1
    oe4.PHOT_CENT = energy
    oe4.FSHAPE = 1
    oe4.RLEN1 = 2.5
    oe4.RLEN2 = 2.5
    oe4.RWIDX1 = 3.5
    oe4.RWIDX2 = 3.5
    oe4.R_LAMBDA = 5000.0
    oe4.T_IMAGE = -1.054
    oe4.T_INCIDENCE = 0
    oe4.T_REFLECTION = 0
    oe4.T_SOURCE = 1.166
    oe4.X_ROT = secondCrystalRot

    #oe5 - torroidal mirror
    oe5 = mirrorType(coatingFile2,torrAngleDeg, sdist = 209.9, imgdist=1601.8,colsimg=1601.8, colssour=1e16,mtype = mirror2type, torrMajor=torrMajor, torrMinor=torrMinor)

    oe6.DUMMY = 1.0
    oe6.FWRITE = 3
    oe6.F_REFRAC = 2
    oe6.F_SCREEN = 1
    oe6.I_SLIT = numpy.array([1, 0, 0, 0, 0, 0, 0, 0, 0, 0])
    oe6.N_SCREEN = 1
    oe6.RX_SLIT = numpy.array([slitx, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe6.RZ_SLIT = numpy.array([slitz, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0])
    oe6.T_IMAGE = 150.0
    oe6.T_INCIDENCE = 0.0
    oe6.T_REFLECTION = 180.0
    oe6.T_SOURCE = -150.0

    #Run SHADOW to create the source
    if traceStart < 1:
        beam.genSource(oe0)
        createdRays = oe0.NTOTALPOINT #find total generated rays (greater than nrays)
        writeCreatedRays(createdRays, createdRaysLog)
    if writeBeam and traceStart < 1:
        print(f'writing {beamFile}.00')
        beam.write(f"{beamFile}.00")
        
    #run optical element 1
    if traceStart < 2:
        print("    Running optical element: %d"%(1))
        beam.traceOE(oe1,1)
    
    if writeBeam and traceStart < 2:
        print(f'writing {beamFile}.01')
        beam.write(f"{beamFile}.01")

    #run optical element 2
    if traceStart < 3 and runMirrors:
        print("    Running optical element: %d"%(2))
        beam.traceOE(oe2,2)

    if writeBeam and traceStart < 3:
        print(f'writing {beamFile}.02')
        beam.write(f"{beamFile}.02")

    #run optical element 3
    if traceStart < 4:
        print("    Running optical element: %d"%(3))
        beam.traceOE(oe3,3)

    if writeBeam and traceStart < 4:
        print(f'writing {beamFile}.03')
        beam.write(f"{beamFile}.03")

    #run optical element 4
    if traceStart < 5:
        print("    Running optical element: %d"%(4))
        beam.traceOE(oe4,4)

    if writeBeam and traceStart < 5:
        print(f'writing {beamFile}.04')
        beam.write(f"{beamFile}.04")

    #run optical element 5
    if traceStart < 6 and runMirrors:
        print("    Running optical element: %d"%(5))
        beam.traceOE(oe5,5)

    if writeBeam and traceStart < 6:
        print(f'writing {beamFile}.05')
        beam.write(f"{beamFile}.05")

    #run optical element 6
    if traceStart < 7 and traceStart < 7:
        print("    Running optical element: %d"%(6))
        beam.traceOE(oe6,6)

    if writeBeam:
        print(f'writing {beamFile}.06')
        beam.write(f"{beamFile}.06")

    result = beam.histo2(1,3, nbins= 101,nolost=1)
    if fname != None:
        fsplit = os.path.splitext(fname)
        basename = fsplit[0]
        ext = fsplit[1]
        efname = f'{basename}_E{ext}'
    else:
        efname = None
    #energyResult = beam.histo1(11,nbins = 101, nolost=  1, write = efname)
    energyResult = beam.histo1(11,nbins = 501, nolost=  1, write = efname, ref = 23)
    #Shadow.ShadowTools.plotxy(beam,1,3,block=block,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    #plt.imshow(beam.rays, aspect = 'auto')
    #plt.show()
    dctToFile(config,configFile)
    createdRays = readCreatedRays(createdRaysLog)
    return result, energyResult,beam, createdRays


if __name__ == '__main__':
    result, eResult, beam, createdRays = run(energy = energy, eRange = eRange, colMirrorRad=firstMirrorAngle, torrAnglemRad=torroidalMirrorAngle, 
                                secondCrystalRot =  secondCrystalRot, writeBeam=writeBeam, nrays=nrays, 
                                traceStart=traceStart, autoStart=autoStart, harmonic=harmonic, coating1=coating1, coating2=coating2,
                                mirror1type=mirror1type,mirror2type=mirror2type, torrMajor=torrMajor, torrMinor=torrMinor,slitx=slitx,slitz=slitz,
                                runMirrors=runMirrors)
    #print(result.keys())
    
    intensity = result['intensity']
    energyIndex = np.abs(fluxEnergy-energy).argmin()
    fluxInitial = fluxDensity[energyIndex]
    
    print()
    NphotonsI = initialPhotons(fluxInitial,eRange,energy) #approximate
    NphotonsF = finalPhotons(NphotonsI,createdRays,intensity) #approximate
    eFWHM = eResult['fwhm']

    if eFWHM == None:
        eFWHMstring = f"energy fwhm: NA\n"
    else:
        eFWHMstring = f"energy fwhm: {eFWHM:.6f} eV\n"
    string = (f"energy: {energy}\n"
    f"harmonic: {harmonic}\n"
    f""
    f"source flux density: {fluxInitial:.6e} photons/(s 0.1%bw)\n"
    f"source total photons/s: {NphotonsI:.6e}\n"
    f"total created rays {createdRays}\n"
    f"created/accepted: {createdRays/nrays:.6f}\n"
    f"intensity: {result['intensity']:.1e}\n" #result parameters: nrays, good_rays, fwhm_h, fwhm_v, fwhm_coordinates_h, fwhm_coordinates_v. #lengths in cm
    f"final photons/s: {NphotonsF:.6e}\n"
    f"fwhm_h: {result['fwhm_h']*10:.6f} mm\n"
    f"fwhm_v: {result['fwhm_v']*10:.6f} mm\n"
    f"{eFWHMstring}")
    print(string)
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    Shadow.ShadowTools.histo1(beam,11,nbins = 201, nolost=  1, ref = 23)
    #Shadow.ShadowTools.histo1(beam,11)

