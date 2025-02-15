#ATE60cm-make_fits
# dir /b/s ..\calibratedA\LFB\*.out | sort | ATE60cm-make_fits.py
#
# Linux 
# find ../calibratedA -name "*.out" -print | sort | python ATE60cm-make_fits_new_vel_grid.py
#
#this script is to make fits files from calibratedA file(*.out)
# 
# Written by Lin Zhenhui(zhlin@pmo.ac.cn)
# Modified on Feb 24, 2024
# V1.0
#
# V1.1 modified by Yiping Ao, tested under Mac OSX
# python ATE60cm-make_fits.py
# script in directory scripts, .out data in directory calibratedA, and .fits in directory fitsData
# V2.0 add the velocity correction. this can be used to import to class files directly. Yan Gong
# V3.0 using CO (4-3) compact sources to correct the pointing error.


from astropy.io import fits
import numpy as np
import pylab as pl
import os,sys,re
from time import gmtime,strftime
from datetime import datetime
from astropy import units as u
from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_body_barycentric, get_body, solar_system_ephemeris
from astropy.time import Time
from astropy.constants import c
from astropy.constants import R_earth
from astropy.coordinates import FK5


def aztoj2000(az,el, obstime):
    # this function is used for ATE coordinate conversion. 
    # Observer's location (replace with your actual latitude, longitude, and elevation)
    observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093.5)
    #observer_location = EarthLocation(lat='37d58m36.2s', lon='96d35m08.6s', height=4820)
    # Example: AZ, EL coordinates of a target
    #azimuth = 24.9857*u.deg#-8./60*u.deg# replace with your azimuth angle in degrees
    #elevation = 37.4345*u.deg  # replace with your elevation angle in degrees
    azimuth = az*u.deg
    elevation = el*u.deg
    # Current UTC time
    #current_time_utc = Time.now()
    #specific_utc_time = Time('2024-01-09T03:58:6.865')
    specific_utc_time = Time(obstime)
    # Create AltAz coordinate object
    #altaz_coord = AltAz(az=azimuth, alt=elevation, location=observer_location, obstime=current_time_utc)
    altaz_coord = AltAz(az=azimuth, alt=elevation, location=observer_location, obstime=specific_utc_time)
    # Convert AltAz coordinates to Equatorial coordinates (RA, DEC)
    #print(altaz_coord.transform_to(Galactic())) # to galactic coordinates
    equatorial_coord = altaz_coord.transform_to(FK5())
    print(equatorial_coord)
    print(equatorial_coord.ra.hms)
    print(equatorial_coord.dec.dms)
    ra = equatorial_coord.ra.value
    dec = equatorial_coord.dec.value
    return ra, dec

def j2000toaz(ra,dec, obstime):
    observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093.5)
    #observer_location = EarthLocation(lat='37d58m36.2s', lon='96d35m08.6s', height=4820)
    # Create a SkyCoord object for the target coordinates
    # target_coord = SkyCoord("10:43:23.075", "-59:29:35.69", unit=(u.hourangle, u.deg), frame='icrs')
    target_coord = SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='icrs')
    specific_utc_time = Time(obstime)
    altaz_frame = AltAz(obstime=specific_utc_time, location=observer_location)
    altaz_coord = target_coord.transform_to(altaz_frame)
    el = altaz_coord.alt.value
    az = altaz_coord.az.value
    return az, el

def pointcorr(ra,dec, obstime, offaz, offel):
    az, el = j2000toaz(ra,dec, obstime)
    newaz = az + offaz
    newel = el + offel
    newra, newdec = aztoj2000(newaz, newel, obstime)
    return newra, newdec


### this is used to estimate the velocity shift by earth rotation, etc.
def cal_vel_corr(myra='+17h20m53.35s', mydec='-35d47m1.50s', mylon='77d06m38.82s', mylat='-80d25m01.6289s', myalt=4093.5, utctime='2025-01-06T19:30:38.96'):
    '''
    myra:     source ra coordinate
    mydec:    source dec coordinate
    mylon:    longtidue of the site
    mylat:    latitude of the site
    myalt:    the elevation of the site
    utctime:  observing utc time
    shift.value  return the correction velocity in units of km/s
    '''
    # Define the coordinates of the observed object
    sc = SkyCoord(ra=myra, dec=mydec, unit=(u.hourangle, u.deg), frame='icrs')
    gl = sc.galactic.l.value
    gb = sc.galactic.b.value
    # Define the location of the observer (telescope)
    ate60 = EarthLocation(lat=mylat, lon=mylon, height=myalt*u.m)
    
    # Define the time of observation
    obs_time = Time(utctime)

    
    # Calculate the Altitude-Azimuth coordinates of the observed object
    barycorr = sc.radial_velocity_correction(obstime=obs_time, location=ate60)  
    barycorr.to(u.km/u.s)

    heliocorr = sc.radial_velocity_correction('heliocentric', obstime=obs_time, location=ate60) 
    heliocorr.to(u.km/u.s) 

    # heliocentric to LSR (Reid+2009)
    # in units of km/s
    Us = 10.27#10.3
    Vs = 15.32#15.3
    Ws = 7.74#7.7

    ## velocity to be corrected 
    cor_val = (Us*np.cos(np.deg2rad(gl))+ Vs*np.sin(np.deg2rad(gl)))*np.cos(np.deg2rad(gb)) + Ws *np.sin(np.deg2rad(gb))
    shift = heliocorr.to(u.km/u.s)+cor_val*u.km/u.s
    return shift.value




lightspeed = 299792458  #m/s

def dms_to_decimal(dms_str):
    # Regular expression to match the DMS format (e.g., " +160d24m34.26s")
    match = re.match(r"([+-]?\d+)d(\d+)m([\d.]+)s", dms_str.strip())
    if match:
        degrees = int(match.group(1))
        minutes = int(match.group(2))
        seconds = float(match.group(3))
        
        # Convert DMS to decimal degrees
        decimal_degrees = degrees + (minutes / 60.0) + (seconds / 3600.0)
        return decimal_degrees
    else:
        raise ValueError(f"Invalid DMS format: {dms_str}")


def initial_Headers():
    hdr = fits.Header()
    hdr.set('SIMPLE',True)
    hdr.set('BITPIX',-32)
    hdr.set('NAXIS',3)
    hdr.set('NAXIS1',32768)
    hdr.set('NAXIS2',1)
    hdr.set('NAXIS3',1)  
    #hdr.set('NAXIS4',1)   
    
    #hdr.set('BLANK',2147483647,'Blanking value')
    #hdr.set('BSCALE',1.0,'Biase scale of intensity')
    #hdr.set('BZERO',0.0,'Bias offset of intensity')
    hdr.set('DATAMIN',-1000.,'min value of data')
    hdr.set('DATAMAX',1000.,'max value of data')
    
    hdr.set('BUNIT','K (Ta*)','unit of spectral line')
    
    hdr.set('CTYPE1','FREQ')
    hdr.set('CUNIT1','Hz','unit of Hertz')
    hdr.set('CRVAL1', 0.,'offset frequency,unit or Hz')
    hdr.set('CDELT1', 0.732421875E+005,'Frequency resolution,unit or Hz')
    hdr.set('CRPIX1',16383,'Number of center channel')
    
    hdr.set('EQUINOX', 2000.,'J2000')
    hdr.set('RADESYS','FK5' ) 
    hdr.set('CTYPE2','RA---GLS','Right ascension of target,unit of degree')
    hdr.set('CUNIT2','deg','unit of degree')
    hdr.set('CRVAL2', 0.,'Reference of RA,unit of degree')
    hdr.set('CDELT2', 0.,'DeltaRA from center pixel,unit of degree')
    hdr.set('CRPIX2',0.,'Pixel of Center' )
    hdr.set('CTYPE3','DEC--GLS','Declination of target,unit of degree')
    hdr.set('CUNIT3','deg','unit of degree')
    hdr.set('CRVAL3', 0.,'Declination,unit of degree')
    hdr.set('CDELT3', 0.,'DeltaDEC from center pixel,unit of degree')
    hdr.set('CRPIX3',0.,'Pixel of Center')

    hdr.set('SITE-LON',0.,'Lontitude of Observation,unit of degree')
    hdr.set('SITE-LAT',0.,'Lattitude of Observation,unit of degree')    
    hdr.set('SITE-ALT',0.,'height of Observation,unit of m')
    
    hdr.set('TELESCOP', 'ATE60','name of telescope')
    hdr.set('OBJECT', 'ORIONA','name of target')
    hdr.set('RA',0.,'Right ascension,unit of degree')
    hdr.set('DEC', 0.,'Declination,unit of degree')
    hdr.set('REF-RA', 0.,'Reference of RA,unit of degree')
    hdr.set('REF-DEC', 0.,'Reference of DEC,unit of degree')
    hdr.set('DELTARA', 0.,'DeltaRA from center pixel,unit of degree')
    hdr.set('DELTADEC', 0.,'DeltaDEC from center pixel,unit of degree')  
    hdr.set('GRID-RA', 0.,'Grid RA')
    hdr.set('GRID-DEC', 0.,'Grid DEC')   
    hdr.set('DELTA-X', 0.,'detla X from main pixel,unit of arcminute')
    hdr.set('DELTA-Y', 0.,'detla Y from main pixel,unit of arcminute')      
    hdr.set('LINE','s','name of spectral line')
    hdr.set('RESTFREQ',0.,'Rest frequency of source,unit of Hz')
    hdr.set('SKY-FREQ',0.,'Sky frequency of source,unit of Hz')  
    hdr.set('LO-FREQ',0.,'LO frequency of Rx,unit of Hz')      
    hdr.set('VELO-LSR',0.,'Velocity of reference channel')
    hdr.set('DELTAV', 0.,'Velocity spacing of channels')
    hdr.set('IMAGFREQ',0.,'Image frequency of signal,unit of Hz')
    hdr.set('TSYS', 0.,'Temperature of receiver system,unit of K')
    hdr.set('TCAL', 0.,'Temperature of calibrator,unit of K')
    hdr.set('OBSTIME',10.,'Integration time on Source,unit of sec')
    hdr.set('SCAN-NUM',1,' Scan number')
    hdr.set('TAU-ATM',0.,'Atmospheric opacity')
    hdr.set('GAINIMAG',1.,'Image sideband gain ratio')
    hdr.set('BEAMEFF',1.,'Beam efficiency')
    hdr.set('FORWEFF', 1.,'Image sideband gain ratio')
    hdr.set('TIMESYS','UTC')
    hdr.set('DATE-OBS','2024-01-01T00:00:00','Date of Observed')
    hdr.set('LST','00:00:00','Sideral time of Observed')
    
    hdr.set('ELEVATIO',0.,'Telescope elevation,unit of deg')
    hdr.set('AZIMUTH',0.,'Telescope azimuth,unit of deg')
    hdr.set('ACT-EL',0.,'Telescope elevation,actually,unit of deg')
    hdr.set('ACT-AZ',0.,'Telescope azimuth,actually,unit of deg') 
    hdr.set('CAL-ACT',0.,'Telescope CAL_ACT,unit of deg')
   
    hdr.set('VCM-PRES',0.,'Pressure of VCM,unit of mbar')  
    hdr.set('5K-TEMP',0.,'Temperature of 5K stage,unit of K') 
    hdr.set('80K-TEMP',0.,'Temperature of 80K stage,unit of K') 
    hdr.set('CPS-TEMP',0.,'Temperature of CPS,unit of K')
    hdr.set('CHD-TEMP',0.,'Temperature of CHD,unit of K')
    hdr.set('MIXER-I', 0.,'Current of mixer bias,unit of uA')
    hdr.set('MIXER-V', 0.,'Voltage of mixer bias,unit of mV')
    hdr.set('LNA-I', 0., 'Current of Low noise amplifier,unit of mA')
    hdr.set('LNA-V', 0.,'Voltage of Low noise amplifier,unit of V')
    hdr.set('IF-TP', 0. ,'Voltage of IF detection output,unit of V')
    hdr.set('LO-V15V', 0.,'Voltage of LO_V15V,unit of V')
    hdr.set('LO-I15V', 0.,'Current of LO_15V,unit of mA')    
    hdr.set('LO-I12V',0.,'Current of LO_12V,unit of mA')
    hdr.set('LO-TEMP',0.,'Temperature of LO,unit of K') 
    hdr.set('LO-YIGV',0.,'Volt of LO YIG,unit of V')     
    hdr.set('DEW-TEMP',0.,'Temperature of Dew,unit of K')     
    hdr.set('HTR-TEMP',0.,'Temperature of HTR,unit of K')  
    
    hdr.set('AMB-TEMP',0.,'Temperature of Ambient,unit of K')  
    hdr.set('RH',0.,'RH of Ambient,unit of %')  
    hdr.set('PRESS',0.,'Pressure of Ambient,unit of hPa') 
    hdr.set('WS',0.,'Wind speed of Ambient,unit of m/s')    
    hdr.set('WD',0.,'Wind direction of Ambient,unit of degree') 
    
    hdr.set('DATE-RED','2024-01-01T00:00:00','Date of reduced')
    hdr.set('DATE-WRT','2024-02-20T00:00:00.000','Date written')   
    hdr.set('ORIGIN','CLASS-Grenoble ')
    hdr.set('OBSERVER', 'NAME1','Observer\'s name')
    hdr.set('ENGINEER', 'NAME2','Engineer\'s name')
    hdr.set('OPERATOR', 'NAME3','Operator\'s name')
   
    return hdr
    
 
def write_fitshdr(file,hdr):
    fp = open(file,'r')
    data = fp.readlines()
    fp.close()
    
    for line in data:
        if re.match('#siteLon*',line):
            P = line.split(':')[1]
            hdr['SITE-LON'] = float(P)
        if re.match('#siteLat*',line):
            P = line.split(':')[1]
            hdr['SITE-LAT'] = float(P)
        if re.match('#siteAlt*',line):
            P = line.split(':')[1]
            hdr['SITE-ALT'] = float(P)   
        if re.match('#srcName*',line):
            P = line.split(': ')[1].split('\n')[0]
            hdr['OBJECT'] = P  
        if re.match('#srcVel*',line):   
            P = line.split(':')[1]
            hdr['VELO-LSR'] = float(P)*1e3 #m/s
        if re.match('#J2000RaOn:*',line):   
            P = line.split(':')[1]
            hdr['RA'] = float(P)
            hdr['CRVAL2'] = float(P)
        if re.match('#J2000DecOn:*',line):   
            P = line.split(': ')[1]
            hdr['DEC'] = float(P)
            hdr['CRVAL3'] = float(P)
        if re.match('#J2000RaOff:*',line):   
            P = line.split(':')[1]
            hdr['REF-RA'] = float(P)            
        if re.match('#J2000DecOff:*',line):   
            P = line.split(':')[1]
            hdr['REF-DEC'] = float(P)           
        if re.match('#deltaRa*',line):   
            P = line.split(':')[1]
            hdr['DELTARA'] = float(P)/3600.  
            hdr['CDELT2'] = float(P)/3600.            
        if re.match('#deltaDec*',line):   
            P = line.split(':')[1]
            hdr['DELTADEC'] = float(P)/3600.   
            hdr['CDELT3'] = float(P)/3600.  
        if re.match('#gridRa*',line): 
            P = line.split(':')[1]
            hdr['GRID-RA'] = int(P)  
        if re.match('#gridDec*',line):   
            P = line.split(':')[1]
            hdr['GRID-DEC'] = int(P)              
        if re.match('#deltaX*',line):   
            P = line.split(':')[1]
            hdr['DELTA-X'] = float(P)                
        if re.match('#deltaY*',line):   
            P = line.split(':')[1]
            hdr['DELTA-Y'] = float(P)               
        if re.match('#srcFreq*',line):   
            P = line.split(':')[1]
            hdr['RESTFREQ'] = float(P)*1e9  
            srcFreq = float(P)*1e9 #+23.*1000./lightspeed*float(P)*1e9
        if re.match('#slName*',line):   
            P = line.split(': ')[1].split('\n')[0]
            hdr['LINE'] = P               
        if re.match('#skyFreq*',line):   
            P = line.split(':')[1]
            hdr['SKY-FREQ'] = float(P)*1e9  
            skyFreq = float(P)*1e9
        if re.match('#loFreq*',line):   
            P = line.split(':')[1]
            hdr['LO-FREQ'] = float(P)*1e9  
            loFreq = float(P)*1e9  
        if re.match('#instChn*',line):   
            P = line.split(':')[1]
            hdr['NAXIS1'] = int(P)             
        if re.match('#instRbw*',line):   
            P = line.split(':')[1]
            hdr['CDELT1'] = float(P)*1e3  #frequency resolution
            chn_bw = float(P)*1e3
        if re.match('#utcDateTime*',line):   
            P = line.split(': ')[1].split('\n')[0].replace(' ','T')
            date_obj = datetime.strptime(P, "%Y-%m-%dT%H:%M:%S.%f")
            hdr['DATE-OBS'] = date_obj.strftime("%Y-%m-%dT%H:%M:%S.%f")            
        if re.match('#lstTime*',line):   
            P = line.split(': ')[1].split('\n')[0]
            hdr['LST'] = P             
        if re.match('#intTime*',line):   
            P = line.split(':')[1]
            hdr['OBSTIME'] = float(P) 
        if re.match('#Scan*',line):   
            P = line.split(':')[1]
            hdr['SCAN-NUM'] = int(P)    
            
        if re.match('#cmdAz*',line):   
            P = line.split(':')[1]
            hdr['AZIMUTH'] = float(P) #dms_to_decimal(P) #float(P)             
        if re.match('#cmdEl*',line):   
            P = line.split(':')[1]
            hdr['ELEVATIO'] = float(P) #dms_to_decimal(P) #float(P)              
        if re.match('#actAz*',line):   
            P = line.split(':')[1]
            hdr['ACT-AZ'] = float(P) #dms_to_decimal(P) # float(P)              
        if re.match('#actEl*',line):   
            P = line.split(':')[1]
            hdr['ACT-EL'] = float(P) #dms_to_decimal(P)# float(P)  
        if re.match('#calAct*',line):   
            P = line.split(':')[1]
            hdr['CAL-ACT'] = float(P)    
            
        if re.match('#Pvcm*',line):   
            P = line.split(':')[1]
            hdr['VCM-PRES'] = float(P)               
        if re.match('#T5K*',line):   
            P = line.split(':')[1]
            hdr['5K-TEMP'] = float(P)               
        if re.match('#T80K*',line):   
            P = line.split(':')[1]
            hdr['80K-TEMP'] = float(P)               
        if re.match('#Tcps*',line):   
            P = line.split(':')[1]
            hdr['CPS-TEMP'] = float(P)              
        if re.match('#Tchd*',line):   
            P = line.split(':')[1]
            hdr['CHD-TEMP'] = float(P)               
        if re.match('#Tcal*',line):   
            P = line.split(':')[1]
            hdr['TCAL'] = float(P)                
        if re.match('#MIX_Vbias*',line):   
            P = line.split(':')[1]
            hdr['MIXER-V'] = float(P)   
        if re.match('#MIX_Ibias*',line):   
            P = line.split(':')[1]
            hdr['MIXER-I'] = float(P)   
        if re.match('#LNA_Vds*',line):   
            P = line.split(':')[1]
            hdr['LNA-V'] = float(P)   
        if re.match('#LNA_Ids*',line):   
            P = line.split(':')[1]
            hdr['LNA-I'] = float(P)   
        if re.match('#IF_Tp*',line):   
            P = line.split(':')[1]
            hdr['IF-TP'] = float(P)   
        if re.match('#LO_V15V*',line):   
            P = line.split(':')[1]
            hdr['LO-V15V'] = float(P) 
        if re.match('#LO_I15V*',line):   
            P = line.split(':')[1]
            hdr['LO-I15V'] = float(P) 
        if re.match('#LO_I12V*',line):   
            P = line.split(':')[1]
            hdr['LO-I12V'] = float(P) 
        if re.match('#LO_Temp*',line):   
            P = line.split(':')[1]
            hdr['LO-TEMP'] = float(P) 
        if re.match('#LO_YIGV*',line):   
            P = line.split(':')[1]
            hdr['LO-YIGV'] = float(P)
        if re.match('#Tdew*',line):   
            P = line.split(':')[1]
            hdr['DEW-TEMP'] = float(P) 
        if re.match('#Thtr*',line):   
            P = line.split(':')[1]
            hdr['HTR-TEMP'] = float(P) 
        if re.match('#Thtr*',line):   
            P = line.split(':')[1]
            hdr['HTR-TEMP'] = float(P) 
        if re.match('#Tamb*',line):   
            P = line.split(':')[1]
            hdr['AMB-TEMP'] = float(P)            
        if re.match('#RH[%]*',line):   
            P = line.split(':')[1]
            hdr['RH'] = float(P) 
        if re.match('#Pressure*',line):   
            P = line.split(':')[1]
            hdr['PRESS'] = float(P) 
        if re.match('#WS[m/s]*',line):   
            P = line.split(':')[1]
            hdr['WS'] = float(P) 
        if re.match('#WD[deg]*',line):   
            P = line.split(':')[1]
            hdr['WD'] = float(P) 
        if re.match('#Spectrum*',line):  
            break

    
    #image frequency from lofreq and srcfreq       
    hdr['IMAGFREQ'] = 2 * loFreq - skyFreq
    #compute I0(CRPIX1)
    # for the OTF, the RA and DEC needs to be recomputed
    hdr["DEC"] = hdr["DEC"] + ((hdr['GRID-DEC']-1)*hdr['deltaDec'])

    hdr["RA"] = hdr["RA"] + ((hdr['GRID-RA']-1)* hdr['deltaRa']/np.cos(np.deg2rad(hdr["DEC"])))

    # let's do a test with +4' +4' 
    # apply a pointing correction here 
    myra, mydec = pointcorr(hdr["RA"], hdr["DEC"], hdr["DATE-OBS"], -1.7/60., 1.6/60.)

    #myra, mydec = aztoj2000(hdr['ACT-AZ'],hdr['ACT-EL'], hdr["DATE-OBS"])

    hdr["RA"]  = myra
    hdr["DEC"] = mydec

    hdr['CRVAL2'] = hdr["RA"] 
    hdr['CRVAL3'] = hdr["DEC"]


    #IF_Sig = abs(skyFreq - loFreq)  # Hz
    #Chn_Sig = IF_Sig/chn_bw        # Hz/Hz
    #hdr['CRPIX1'] = int(Chn_Sig + 0.5) # Signal Channel (I0) # someting wrong with the reference number and rest frequency
    # frequency shift for CO frequency
    #cofreq = 461.0407682*1e9
    cc = SkyCoord(ra=hdr["RA"]*u.degree, dec=hdr["DEC"]*u.degree, frame='icrs')
    
    shift_velo = cal_vel_corr(myra=cc.ra.to_string(unit='hour'), mydec=cc.dec.to_string(unit='deg'), mylon='77d06m38.82s', mylat='-80d25m01.6289s', myalt=4093.5, utctime=hdr["DATE-OBS"])
    print(shift_velo)
    
    shift_freq = shift_velo*1000./lightspeed*hdr['restfreq']
    
    hdr['CRPIX1'] = 1
    hdr['restfreq'] = loFreq - shift_freq
    # Compute Dv and df
    df = 0.
    if(srcFreq >= loFreq):
        df = chn_bw  # Hz
    else:
        df = -1. * chn_bw  #Hz
    
    Dv = -1.* lightspeed * df/srcFreq     #m/s
    hdr['DELTAV'] = Dv
    hdr['CDELT1'] = df
    
    
    
    
    time_str = strftime("%Y-%m-%dT%H:%M:%S", gmtime())
    hdr['DATE-RED'] = time_str 
    hdr['DATE-WRT'] = time_str 
    
    #hdr['NAXIS'] = 3
    #hdr['NAXIS1'] = 32768
    #hdr['NAXIS2'] = 1
    #hdr['NAXIS3'] =1  
    
if __name__ == '__main__':

    hdr = initial_Headers()
    #file_path = '../calibratedA/'
    #files = os.listdir(file_path)
    
    files = []
    for line in sys.stdin:
        P = line.split('\n')[0]
        files.append(P)
    
    files = np.array(files)
    N = len(files)

    for i in range(N):
    
        source_dir = os.path.dirname(files[i])
        dest_dir = source_dir.replace('calibratedA', 'fitsData')
        if not os.path.isdir(dest_dir):
            os.makedirs(dest_dir)
        basename = os.path.basename(files[i])
        newfile = basename.replace('out', 'fits')
        outf = os.path.join(dest_dir, newfile)
        #file_name = os.path.join(file_path,file_name)

        data = np.loadtxt(files[i])
        sp = data[:, 1]
        sp = sp.astype(np.float32)
        maxV = np.max(sp)
        minV = np.min(sp)
        if maxV > 1e9:
            maxV = 1e9
        if minV < -1e9:
            minV = -1e-9
        sp = sp.reshape(1,1,len(sp))
        hdu = fits.PrimaryHDU(sp)
#       hdr = hdu.header
        write_fitshdr(files[i],hdr)
        hdr['DATAMAX'] = maxV
        hdr['DATAMIN'] = minV
        hdu.header = hdr
        
        #print(hdu.header)

        hdu.writeto(outf, overwrite=True)
        #print(outf)

    
    
    
    
    
    
    
    
    
