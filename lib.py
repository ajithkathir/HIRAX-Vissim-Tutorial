#Author : Ajith Sampath
#Date : 20-10-2024


import matplotlib.pyplot as plt
import numpy as np
import matplotlib
#from pygdsm import GlobalSkyModel, GSMObserver
import datetime as dt
import healpy as hp
from tqdm import tqdm
import itertools


from casatools import simulator, image, table, coordsys, measures, componentlist, quanta, ctsys,vpmanager
from casatasks import importfits, tclean, ft, imhead, listobs, exportfits, flagdata, bandpass, applycal,imtrans,uvsub
from casatasks.private import simutil
import sys
import os
from astropy.io import fits
from astropy.wcs import WCS


# Instantiate all the required tools
sm = simulator()
ia = image()
tb = table()
cs = coordsys()
me = measures()
qa = quanta()
cl = componentlist()
mysu = simutil.simutil()
vp = vpmanager()


def get_index_for_antenna_pair(antenna_pair, num_antennas):
    # Generate all unique antenna pairs (combinations)
    all_pairs = list(itertools.combinations(range(num_antennas), 2))
    
    # Find the index of the given antenna pair
    try:
        index = all_pairs.index(antenna_pair)
        return index
    except ValueError:
        return "Antenna pair not found"


def makeMSFrame(array_coord,msname = 'sim_data.ms',lat=31,start_freq='400MHz',deltafreq='390.625kHz',n_chan=1024,ra='10h00m00.0s',dec='-31d00m00.0s'):
    """
    Construct an empty Measurement Set that has the desired observation setup.
    """
    #print("\nCreating measurement set....")

    os.system('rm -rf '+msname)

    ## Open the simulator
    sm.open(ms=msname);

    ## Set the antenna configuration
    
    
    E = array_coord[:,0]
    N = array_coord[:,1]
    U = array_coord[:,2]

 
    ant = np.arange(0,len(E),1)
    #ant = np.arange(0,2,1)
    antid = (np.zeros(len(ant)))

    x = -1*N*np.sin(lat) + U*np.cos(lat)
    y = E
    z = N*np.cos(lat) + U*np.sin(lat)
    dia = 6
    dia_col = np.zeros((len(x)))+dia
    sm.setconfig(telescopename='NGVLA',
                     x=x,
                     y=y,
                     z=z,
                     dishdiameter=dia_col,
                     mount=['alt-az'],
                     antname=ant,
                     coordsystem='global',
                     referencelocation=me.observatory('vla'));

    ## Set the polarization mode (this goes to the FEED subtable)
    sm.setfeed(mode='perfect R L', pol=['']);

    ## Set the spectral window and polarization (one data-description-id).
    ## Call multiple times with different names for multiple SPWs or pol setups.
    sm.setspwindow(spwname="hiraxband",
                   freq=start_freq,
                   deltafreq=deltafreq,
                   freqresolution=deltafreq,
                   nchannels=n_chan,
                   stokes='RR');

    ## Setup source/field information (i.e. where the observation phase center is)
    ## Call multiple times for different pointings or source locations.
    sm.setfield( sourcename="fake",
                 sourcedirection=me.direction(rf='J2000', v0=ra,v1=dec));

    ## Set shadow/elevation limits (if you care). These set flags.
    #sm.setlimits(shadowlimit=0.01, elevationlimit='1deg');

    ## Leave autocorrelations out of the MS.
    sm.setauto(autocorrwt=0.0);

    ## Set the integration time, and the convention to use for timerange specification
    ## Note : It is convenient to pick the hourangle mode as all times specified in sm.observe()
    ##        will be relative to when the source transits.
    sm.settimes(integrationtime='10s',
                usehourangle=True,
                referencetime=me.epoch('UTC','2019/10/4/00:00:00'));

    ## Construct MS metadata and UVW values for one scan and ddid
    ## Call multiple times for multiple scans.
    ## Call this with different sourcenames (fields) and spw/pol settings as defined above.
    ## Timesteps will be defined in intervals of 'integrationtime', between starttime and stoptime.
    sm.observe(sourcename="fake",
               spwname='hiraxband',
               starttime='-5s',
               stoptime='+5s');

    ## Close the simulator
    sm.close()

    ## Unflag everything (unless you care about elevation/shadow flags)
    flagdata(vis=msname,mode='unflag')
    #print("Done!\n")


 


def makeEmptyImage(imname_true='sim_onepoint_true.im',ra='19h00m00.0s',dec='-31d00m00.0s',imsize=1501,n_chan=1024,ang_res='0.1deg',start_freq='400MHz',deltafreq='390.625kHz'):
    ## Define the center of the image
    radir = ra
    decdir = dec
    #print("\nCreating empty image...")
    ## Make the image from a shape
    ia.close()
    ia.fromshape(imname_true,[imsize,imsize,1,n_chan],overwrite=True)

    ## Make a coordinate system
    cs=ia.coordsys()
    cs.setunits(['rad','rad','','Hz'])
    cell_rad=qa.convert(qa.quantity(ang_res),"rad")['value']
    cs.setincrement([-cell_rad,cell_rad],'direction')
    cs.setreferencevalue([qa.convert(radir,'rad')['value'],qa.convert(decdir,'rad')['value']],type="direction")
    cs.setreferencevalue(start_freq,'spectral')
    cs.setreferencepixel([0],'spectral')
    cs.setincrement(deltafreq,'spectral')

    ## Set the coordinate system in the image
    ia.setcoordsys(cs.torecord())
    ia.setbrightnessunit("Jy/pixel")
    ia.set(0.0)
    ia.close()
    #print("Done!\n")
### Note : If there is an error in this step, subsequent steps will give errors of " Invalid Table Operation : SetupNewTable.... imagename is already opened (is in the table cache)"
## The only way out of this is to restart the kernel (equivalent to exit and restart CASA).
## Any other way ?


def copyModel2modelcolumn(model_ms='sim_data_1.ms',main_ms='sim_data.ms'):
    tb.open(model_ms,nomodify=False);
    moddata = tb.getcol(columnname='DATA');
    tb.close()
    tb.open(main_ms,nomodify=False);
    tb.putcol(columnname='MODEL_DATA',value = moddata);
    #tb.putcol(columnname='CORRECTED_DATA',value=moddata);
    tb.close();


def twoD_Gaussian(x,y,amp,sx,sy,xo, yo, theta, off):
    xo = float(xo)
    yo = float(yo)
    x,y=np.meshgrid(x,y)
    a = (np.cos(theta)**2)/(2*sx**2) + (np.sin(theta)**2)/(2*sy**2)
    b = -(np.sin(2*theta))/(4*sx**2) + (np.sin(2*theta))/(4*sy**2)
    c = (np.sin(theta)**2)/(2*sx**2) + (np.cos(theta)**2)/(2*sy**2)
    return off + amp*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2))) 



def predictSim(msname='sim_data.ms',
                imname='sim_onepoint_true.im'):
    #print("\nPredicting visibility from given image...")
    ## Open an existing MS Frame
    sm.openfromms(msname)
    
    # Predict from a model image
    sm.predict( imagename = imname, incremental=False)
    # Close the tool
    sm.close()
    #print("Done!\n")


def createGSMsky(nchan,imsize,year,month,day,hour,minute,start_freq,end_freq):
    im = np.zeros((nchan,imsize,imsize))
    freq = np.linspace(start_freq,end_freq,nchan)
    counts=0
    for f in freq:
        ov = GSMObserver()
        ov.lon = obs_lon
        ov.lat = obs_lat
        ov.elev = obs_elev
        ov.date = dt.datetime(year,month,day,hour,minute)
        ov.generate(f)
        d=ov.view()
        plt.close()
        d1 = hp.orthview(d, half_sky=True,return_projected_map=True, xsize=imsize)
        d1 = np.nan_to_num(d1.data, neginf=0) 
        plt.close()
        # Convert to a regular array, replacing masked values with a fill value (e.g., 0)
        d1[d1<0] = 0
        im[counts] = d1
        counts=counts+1
    return im;



def simulate(skyimage,beam,array_coord,latitude,start_freq,chan_width,n_chan,RA,DEC,image_angular_resolution,imsize,dirname):  

    os.mkdir(dirname)
    os.chdir(dirname)

    
    sky_beam = np.zeros_like(skyimage)
    for i in range(skyimage.shape[0]):
        sky_beam[i] = skyimage[i]*beam[i]
    
    makeMSFrame(array_coord,msname ="test.ms",lat=latitude,start_freq=str(start_freq)+'MHz',deltafreq=chan_width,n_chan=n_chan,ra=RA,dec=DEC)
    
    
    ## Make an empty CASA image
    makeEmptyImage(imname_true="empty.im",ra=RA,dec=DEC,imsize=imsize,n_chan=n_chan,ang_res=image_angular_resolution,start_freq=str(start_freq)+'MHz',deltafreq=chan_width)
    
    
    exportfits(imagename="empty.im",fitsimage="empty.fits")
    
    hdu = fits.open("empty.fits")
    header = hdu[0].header
    
    hdu1=fits.PrimaryHDU(sky_beam,header=header)
    
    hdu1.writeto("im1.fits",overwrite=True)
    
    #print("Done!\n")
    
    importfits(fitsimage="im1.fits",imagename="PtSrc.im",beam=['0.35arcsec', '0.24arcsec', '25deg'])
    
    imtrans(imagename="PtSrc.im",outfile="PtSrc1.im",order='0132')
    
    predictSim(msname="test.ms",imname="PtSrc1.im")
    
    tb.open("test.ms")
    vis = tb.getcol("DATA")[0]
    tb.close()

    #print("\n Removing unwanted files and cleaning...")
    #os.system("rm -rf *")
    os.chdir("..")
    os.system("rm -rf "+dirname)

    return vis;

def calculate_baselines(num_antennas):
    # Formula for number of baselines
    return num_antennas * (num_antennas - 1) // 2


def simulate_drift(data,beam,array_coord,latitude,start_freq,
               chan_width,n_chan,RA,DEC,image_angular_resolution,
               imsize,dirname):
    bl = calculate_baselines(array_coord.shape[0])
    Vis = np.zeros((data.shape[0],n_chan,bl),dtype=np.complex128)
    for i in tqdm(range(data.shape[0]), desc="computing..."):
        Vis[i] = simulate(data[i],beam,array_coord,latitude,start_freq,
                   chan_width,n_chan,RA,DEC,image_angular_resolution,
                   imsize,dirname)

    return Vis;
   
def plot_waterfall(V1,antenna_pair,array_coord):
    ind = get_index_for_antenna_pair(antenna_pair, array_coord.shape[0])
    
    print("Plotting baseline index..",ind)
    plt.pcolor(np.log10(np.abs(V1[:,:,ind].T)))
    plt.xlabel("Ntime")
    plt.ylabel("Nchan")
    plt.colorbar()

