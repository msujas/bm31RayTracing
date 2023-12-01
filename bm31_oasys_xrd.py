#
# Python script to run shadow3. Created automatically with ShadowTools.make_python_script_from_list().
#
import Shadow
import numpy
np = numpy
from srxraylib.sources import srfunc
import matplotlib.pyplot as plt


energy = 49000
monoEnergy = 49000
focalEnergy = 49000
torroidalMirrorAngle = 3 #mrad from surface
nrays = 1000000
eRange = 200

def mradSurface_to_degNorm(mrad):
    return 90-mrad*180/(numpy.pi*1000)

# write (1) or not (0) SHADOW files start.xx end.xx star.xx
iwrite = 0
dspacing = 3.13379
def energyToTheta(energy):
    wavelength = 12398.47/energy
    return np.arcsin(wavelength/(2*dspacing))*180/np.pi

def saggitalRadius(f1,f2,energy):
    theta = energyToTheta(energy)*np.pi/180
    return (2*f1*f2*np.sin(theta)/(f1+f2))

def meridionalRadius(f1,f2,energy):
    theta = energyToTheta(energy)*np.pi/180
    return ((f1+f2)**2 - 4*f1*f2*np.sin(theta)**2)**0.5/(np.sin(2*theta))

f1 = 3032.8
f2 = 1811.6

def run(energy = 49000,  monoEnergy = 49000, focalEnergy = 49000, meridionalDist = 1000000, iwrite = 0, writeBeam = 1, 
        nrays = 1000000, traceStart = 0):

    #
    # initialize shadow3 source (oe0) and beam
    #
    beam = Shadow.Beam()
    oe0 = Shadow.Source()
    oe1 = Shadow.OE()
    oe2 = Shadow.OE()
    oe3 = Shadow.OE()
    if traceStart > 0:
        beam.load(f'star.{traceStart-1:02d}')

    #
    # Define variables. See meaning of variables in: 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/source.nml 
    #  https://raw.githubusercontent.com/srio/shadow3/master/docs/oe.nml
    #

    #generating source trajectory files
    traj, pars = srfunc.wiggler_trajectory(
    b_from            =1,
    inData            ="C:/Users/kenneth1a/Documents/mirrors/SW_2PA.txt",
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
    outFile        =b'C:\\Users\\kenneth1a\\Documents\\mirrors\\xshwig.sha',
    elliptical     =False)

    #oe0 - wiggler source
    oe0.BENER = 6.0
    oe0.CONV_FACT = 100.0
    oe0.EPSI_DX = 89.4
    oe0.EPSI_DZ = -104.8
    oe0.EPSI_X = 2.16e-08
    oe0.EPSI_Z = 5e-10
    oe0.FDISTR = 0
    oe0.FILE_BOUND = b'C:\\Users\\kenneth1a\\Documents\\mirrors\\myslit.dat'
    oe0.FILE_TRAJ = b'C:\\Users\\kenneth1a\\Documents\\mirrors\\xshwig.sha'
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
    oe0.ISTAR1 = 5676561
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

    oe1.DUMMY = 1.0
    oe1.FHIT_C = 1
    oe1.FILE_REFL = b'C:/Users/kenneth1a/Documents/mirrors/Si5_55.111'
    oe1.FWRITE = 1
    oe1.F_CENTRAL = 1
    oe1.F_CRYSTAL = 1
    oe1.F_MOVE = 1
    oe1.PHOT_CENT = monoEnergy
    oe1.RLEN1 = 2.5
    oe1.RLEN2 = 2.5
    oe1.RWIDX1 = 3.5
    oe1.RWIDX2 = 3.5
    oe1.R_LAMBDA = 5000.0
    oe1.T_IMAGE = 5.0
    oe1.T_SOURCE = f1

    oe2.ALPHA = 180.0
    oe2.DUMMY = 1.0
    oe2.FHIT_C = 1
    oe2.FILE_REFL = b'C:/Users/kenneth1a/Documents/mirrors/Si5_55.111'
    oe2.FMIRR = 3
    oe2.FWRITE = 1
    oe2.F_CENTRAL = 1
    oe2.F_CRYSTAL = 1
    oe2.F_EXT = 1
    oe2.F_MOVE = 1
    oe2.F_TORUS = 2
    oe2.PHOT_CENT = monoEnergy
    oe2.RLEN1 = 2.5
    oe2.RLEN2 = 2.5
    oe2.RWIDX1 = 3.5
    oe2.RWIDX2 = 3.5
    oe2.R_LAMBDA = 5000.0
    oe2.R_MAJ = meridionalRadius(f1,meridionalDist,focalEnergy)
    oe2.R_MIN = saggitalRadius(f1,f2,focalEnergy)
    oe2.T_IMAGE = f2
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
    

    if iwrite:
        oe0.write("start.00")
    if traceStart < 1:
        beam.genSource(oe0)

    if iwrite:
        oe0.write("end.00")
    if writeBeam and traceStart < 1:
        beam.write("star.00")


    #
    #run optical element 1
    #
    print("    Running optical element: %d"%(1))
    if iwrite:
        oe1.write("start.01")
    if traceStart < 2:
        beam.traceOE(oe1,1)

    if iwrite:
        oe1.write("end.01")
    if writeBeam and traceStart < 2:
        beam.write("star.01")


    #
    #run optical element 2
    #
    print("    Running optical element: %d"%(2))
    if iwrite:
        oe2.write("start.02")

    beam.traceOE(oe2,2)

    if iwrite:
        oe2.write("end.02")
    if writeBeam:
        beam.write("star.02")

    #
    #run optical element 3
    #
    print("    Running optical element: %d"%(3))
    if iwrite:
        oe3.write("start.03")

    #beam.traceOE(oe3,3)

    #if iwrite:
    #    oe3.write("end.03")
    #    beam.write("star.03")
    result = beam.histo2(1,3, nbins= 101,nolost=1)
    eResult = beam.histo1(11,nbins = 201, nolost=  1, ref = 23)
    #Shadow.ShadowTools.plotxy(beam,1,3,block=block,nbins=101,nolost=1,title="Real space")
    # Shadow.ShadowTools.plotxy(beam,1,4,nbins=101,nolost=1,title="Phase space X")
    # Shadow.ShadowTools.plotxy(beam,3,6,nbins=101,nolost=1,title="Phase space Z")
    #plt.imshow(beam.rays, aspect = 'auto')
    #plt.show()
    return result, eResult,beam

if __name__ == '__main__':
    result, eResult, beam = run(energy = energy, iwrite=iwrite, monoEnergy=monoEnergy,nrays = nrays, focalEnergy=focalEnergy)
    print('intensity:',result['intensity'])
    print('fwhm_h', result['fwhm_h']*10, 'mm')
    print('fwhm_v', result['fwhm_v']*10, 'mm')
    print('energy fwhm', eResult['fwhm'], 'eV')
    Shadow.ShadowTools.plotxy(beam,1,3,nbins=101,nolost=1,title="Real space")
    Shadow.ShadowTools.histo1(beam,11,nbins = 201, nolost=  1, ref = 23)
    #plt.stairs(eResult['histogram'],eResult['bins'])
    #plt.xlabel('energy (eV)')
    #plt.ylabel('intensity')
    #plt.show()
    #print(eResult)
    #print(result)
    #plt.imshow(result['histogram'].transpose(), extent = [result['bin_h_edges'][0],result['bin_h_edges'][-1],
    #                                                      result['bin_v_edges'][0],result['bin_v_edges'][-1]], aspect = 'auto')
    #plt.show()