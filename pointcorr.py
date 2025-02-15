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


az1, el1 = j2000toaz(258.0663326, -38.5434273, obstime="2025-01-14T15:00:10")
az2, el2 = j2000toaz(258.0956096, -38.5148204, obstime="2025-01-14T15:00:10")

print("pointing correction in AZ and EL (arcmin)")
print((az1-az2)*60, (el1-el2)*60)
