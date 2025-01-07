'''
This is desgined for ATE30/ATE60 to show the source distribution in sky.
contact: ygong@pmo.ac.cn, gongyan2444@gmail.com
update: 26-02-2024
add the Equatorial plane: 12-03-2024
update: 10-11-2024 
add the optical guide stars and sources to be observed
'''
import sys
from PyQt5.QtWidgets import (QApplication, QWidget, QLabel, QLineEdit, QPushButton, QHBoxLayout, QVBoxLayout, QMenuBar, QFileDialog)
from PyQt5.QtGui import QFont
from PyQt5.QtCore import QSize
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
import matplotlib.pyplot as plt
from astropy.coordinates import EarthLocation, AltAz, get_sun, get_body
from astropy.time import Time
from astropy import units as u
from astropy.coordinates import AltAz, SkyCoord
import numpy as np



def ate30_aztoj2000(az,el, obstime):
    # this function is used for ATE30 coordinate conversion. 
    # Observer's location (replace with your actual latitude, longitude, and elevation)
    observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093)
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
    print(altaz_coord.transform_to(Galactic())) # to galactic coordinates
    equatorial_coord = altaz_coord.transform_to(FK5())
    print(equatorial_coord)
    print(equatorial_coord.ra.hms)
    print(equatorial_coord.dec.dms)
    ra = equatorial_coord.ra.value
    dec = equatorial_coord.dec.value
    return ra, dec

def ate30_j2000toaz(ra,dec, obstime):
    observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093)
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

def milkyway(obstime):
    gl = np.arange(360)
    gb = gl * 0 
    target_coord = SkyCoord(gl*u.deg, gb*u.deg, frame='galactic')
    target_coord_eq = target_coord.transform_to('icrs')
    az, el = ate30_j2000toaz(target_coord_eq.ra.value, target_coord_eq.dec.value, obstime)
    return az, el

#def equatorial_plane(obstime):
#   ra = np.arange(241)*0.1
#    dec = ra * 0 
#    target_coord = SkyCoord(ra, dec,  unit=(u.hourangle, u.deg), frame='icrs')
#    az, el = ate30_j2000toaz(target_coord.ra.deg, target_coord.dec.deg, obstime)
#    #az, el = ate30_j2000toaz(ra, dec, obstime)
#    return az, el
   


def pltsource(ax, sounm,ra, dec, times):
    '''
    black: current position
    red:   position in 1 hour
    '''
    target_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    az, el = ate30_j2000toaz(target_coord.ra.deg, target_coord.dec.deg, times)
    az2, el2 = ate30_j2000toaz(target_coord.ra.deg, target_coord.dec.deg, times+1*u.hour)
    if (el > 0):
        ax.scatter(az,el,marker='o',color='k',s=12)
        #plt.scatter(az2,el2,marker='o',color='r',s=12)
        #plt.plot([az,az2], [el, el2],ls='dotted', color='gray')
        ax.text(az,el,sounm, ha='left')


def pltstar(ax, sounm,ra, dec, times):
    '''
    plot the optical guide stars
    green: current position
    '''
    target_coord = SkyCoord(ra, dec, unit=(u.hourangle, u.deg), frame='icrs')
    az, el = ate30_j2000toaz(target_coord.ra.deg, target_coord.dec.deg, times)
    az2, el2 = ate30_j2000toaz(target_coord.ra.deg, target_coord.dec.deg, times+1*u.hour)
    if (el > 0):
        ax.scatter(az,el,marker='o',color='tab:green',s=12)
        #plt.scatter(az2,el2,marker='o',color='r',s=12)
        #plt.plot([az,az2], [el, el2],ls='dotted', color='gray')
        ax.text(az,el,sounm, ha='left')
        

def pltplant(ax, sounm, times):
    '''
    black: current position
    red:   position in 1 hour
    sounm: planets like Jupiter
    '''
    #observer_location = EarthLocation(lat='37d58m36.2s', lon='96d35m08.6s', height=4820)
    observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093)
    planet_position = get_body(sounm,times)
    planet_altaz = planet_position.transform_to(AltAz(obstime=times, location=observer_location))
    if (planet_altaz.alt.deg > 0):
        ax.scatter(planet_altaz.az.deg,planet_altaz.alt.deg,marker='o', color='tab:purple',s=20)
        ax.text(planet_altaz.az.deg, planet_altaz.alt.deg, sounm, ha='left',color='tab:purple')

class ATE60VIEW(QWidget):
    def __init__(self):
        super().__init__()

        self.initUI()

    def initUI(self):
        # Create widgets
        lbl_d = QLabel("Delay (hrs):", self)
        self.txt_d = QLineEdit(self)
        self.txt_d.setText("0")
        lbl_d.setFont(QFont("Arial", 16))

        hbox_d = QHBoxLayout()
        hbox_d.addWidget(lbl_d)
        hbox_d.addWidget(self.txt_d)
##################################################################################
        # File input to specify star file
        lbl_file = QLabel("Guide star:", self)
        self.txt_file = QLineEdit(self)
        self.txt_file.setText("star.txt")  # Default file name
        lbl_file.setFont(QFont("Arial", 16))


        # Button to browse files
        self.browse_button = QPushButton("Browse", self)
        self.browse_button.setFont(QFont("Arial", 16))
        self.browse_button.setFixedSize(100, 30)


        # File input to specify target file
        lbl_target = QLabel("Targets:", self)
        self.txt_target = QLineEdit(self)
        self.txt_target.setText("sou.txt")  # Default file name
        lbl_target.setFont(QFont("Arial", 16))

        # Button to browse files
        self.browse_button_target = QPushButton("Browse", self)
        self.browse_button_target.setFont(QFont("Arial", 16))
        self.browse_button_target.setFixedSize(100, 30)
        
        #
        hbox_file = QHBoxLayout()
        hbox_file.addWidget(lbl_file)
        hbox_file.addWidget(self.txt_file)
        hbox_file.addWidget(self.browse_button)

        hbox_file.addWidget(lbl_target)
        hbox_file.addWidget(self.txt_target)
        hbox_file.addWidget(self.browse_button_target)

        #
        #hbox_file = QHBoxLayout()
        #hbox_file.addWidget(lbl_file)
        #hbox_file.addWidget(self.txt_file)
        #hbox_file.addWidget(self.browse_button)
##################################################################################
        
 ##################################################################################       
        # plot button
        self.plot_button = QPushButton("Refresh Skyview", self)
        self.plot_button.setFont(QFont("Arial", 16))

        vbox = QVBoxLayout()
        vbox.addLayout(hbox_d)
        vbox.addLayout(hbox_file)
        vbox.addWidget(self.plot_button)

##################################################################################        
        # Matplotlib Figure and Canvas
        self.figure = plt.figure()
        self.canvas = FigureCanvas(self.figure)
        vbox.addWidget(self.canvas)

        # Connect signals and slots
        self.plot_button.clicked.connect(self.plot_skyview)
        self.browse_button.clicked.connect(self.browse_file)
        self.browse_button_target.clicked.connect(self.browse_target)

        # Create a layout for the canvas and toolbar
        self.toolbar = NavigationToolbar(self.canvas, self)
        self.toolbar.setIconSize(QSize(16, 16))
        vbox.addWidget(self.toolbar)

        self.plot_skyview()

        # Set window properties
        self.setLayout(vbox)
        self.setWindowTitle('ATE60 Skyview')
        self.setGeometry(100, 100, 1000, 700)

    # define browse_file function

    def plot_skyview(self):
        self.figure.clear()  # Clear the previous plot
        #observer_location = EarthLocation(lat='37d58m36.2s', lon='96d35m08.6s', height=4820)
        observer_location = EarthLocation(lat='-80d25m01.6289s', lon='77d06m38.82s', height=4093)
        delay = 0 
        #print(self.txt_d.text().lstrip('-'))
        #if self.txt_d.text().lstrip('-').isdigit():
        ##self.txt_d.text().lstrip('-').isascii():
        #    delay = float(self.txt_d.text())
        #else:
        #    print("invalid input values")
        try:
            delay = float(self.txt_d.text())
        except ValueError:
            print("Invalid input values")


        times = Time.now()+delay*u.hour
        sun_position = get_sun(times)
        moon_position = get_body('moon',times)

        sun_altaz = sun_position.transform_to(AltAz(obstime=times, location=observer_location))
        moon_altaz = moon_position.transform_to(AltAz(obstime=times, location=observer_location))
           
        ax = self.figure.add_subplot(111)
        if (sun_altaz.alt.deg > 0):
            ax.scatter(sun_altaz.az.deg,sun_altaz.alt.deg,marker='o',color='r',s=60, label='Sun')
        if (moon_altaz.alt.deg > 0):
            ax.scatter(moon_altaz.az.deg,moon_altaz.alt.deg,marker='o', color='b',s=20,label='Moon')

        pltplant(ax, "Jupiter", times)
        pltplant(ax, "Saturn", times)
        pltplant(ax, "Mercury", times)
        pltplant(ax, "Mars", times)
        pltplant(ax, "Venus", times)
        pltplant(ax, "Uranus", times)
        pltplant(ax, "Neptune", times)

        # Additional plotting logic here...
        #pltsource(ax, "NGC253"   , "00:47:33.134","-25:17:19.68", times)
        #pltsource(ax, "N159W", "05:39:36.0", "-69:45:35.0", times) # circumpolar
        #pltsource(ax, "N113", "05:13:17.4", "-69:22:22.0", times)  # circumpolar
        #pltsource(ax, "30Dor", "05:38:38.0", "-69:06:47.9", times)  # circumpolar

        #high-mass
        #targetfile = 'sou.txt'
        targetfile = self.txt_target.text()
        targetname, target_ra, target_dec = np.loadtxt(targetfile, dtype=str, unpack=True, usecols=(0,1,2))
        for i in np.arange(len(targetname)):
            pltsource(ax, targetname[i], target_ra[i], target_dec[i], times)
        #pltsource(ax, "NGC6334I", "17:20:53.35","-35:47:01.5", times)
        #pltsource(ax, "Carina"  , "10:43:23.075","-59:29:35.69", times)  # circumpolar
        #pltsource(ax, "SgrA"    , "17:45:40.036",  "-29:00:28.17", times)
        #pltsource(ax, "SgrB2"    , "17:47:19.8",  "-28:22:17.0", times)
        #pltsource(ax, "RCW79"  , "13:40:18",  "-61:44:12", times)
        #pltsource(ax, "RCW120", "17:12:18", "-38:27:43", times)
        #pltsource(ax, "G305", "13:13:43", "-62:30:38", times)

        starfile = self.txt_file.text()

        # optical guide stars
        #data = np.loadtxt('star.txt', dtype=str)
        starname, star_ra, star_dec = np.loadtxt(starfile, dtype=str, unpack=True, usecols=(0,1,2))
        for i in np.arange(len(starname)):
            pltstar(ax, starname[i], star_ra[i], star_dec[i], times)

        #pltsource(ax, "SgrC"     , "17:44:47.0",  "-29:28:24.8", times)

        gaz, gel = milkyway(times)
        ax.scatter(gaz, gel, marker='.', s=1, color='tab:orange', label="Galactic plane")
        #e_az, e_el = equatorial_plane(times)
        #ax.scatter(e_az, e_el, marker='.', s=1, color='tab:pink')

        # Plotting logic here...
        ax.axhline(y=30, xmin=-30, xmax=390, ls='--')
        ax.legend()
        ax.set_xlabel('Azimuth (degree)')
        ax.set_ylabel('Elevation (degree)')
        ax.set_title("ATE60 UTC:"+times.iso)
        ax.set_xlim(-30, 390)
        ax.set_ylim(0, 90)

        # Draw the plot
        self.canvas.draw()

    def browse_file(self):
        # Open a file dialog to let the user choose a file
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", "Text Files (*.txt);;All Files (*)")
        
        if file:
            # Set the selected file path in the text input
            self.txt_file.setText(file)
            print(f"Selected file: {file}")

    def browse_target(self):
        # Open a file dialog to let the user choose a file
        file, _ = QFileDialog.getOpenFileName(self, "Select File", "", "Text Files (*.txt);;All Files (*)")
        
        if file:
            # Set the selected file path in the text input
            self.txt_target.setText(file)
            print(f"Selected file: {file}")


if __name__ == '__main__':
    app = QApplication(sys.argv)
    ex = ATE60VIEW()
    ex.show()
    sys.exit(app.exec_())
