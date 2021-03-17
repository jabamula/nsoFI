#!/usr/bin/env python
# coding: utf-8
# Jari Backman Slooh Astronomer 27.5.2020
# astronomical calculations based on Paul Schlyter's web site <a href="https://stjarnhimlen.se/comp/ppcomp.html
# ver 1.2
# - added Messier, Caldwell and Variable Stars combobox
# ver 1.1
# - added the choice of saving the image
# - corrected the date explanation on x-axis
#
from PyQt5 import QtCore, QtGui, QtWidgets
from PyQt5.QtCore import QDate, QDateTime
import math
import re
from datetime import datetime, timedelta
import datetime as dt
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 
from matplotlib.dates import DateFormatter
import csv

# arcsin and arccos
def fnasin(x):
    return math.atan(float(x)/math.sqrt(1-float(x)*float(x)));

def fnacos(x): 
    return (pi/2 - fnasin(float(x)));

# Trig. functions in degrees
def fnsind(x):
    return math.sin(float(x)/radeg);

def fncosd(x):
    return math.cos(float(x)/radeg);

def fntand(x):
    return math.tan(float(x)/radeg);

def fnasind(x):
    if x == 1:
        return pi / 2;
    elif x == -1:
        return - pi/2;
    else:
        return radeg*math.atan(float(x)/math.sqrt(1-float(x)*float(x)));

def fnacosd(x):
    return 90 - fnasind(float(x));

def fnatand(x):
    return radeg*math.atan(float(x));
#
# arctan in all four quadrants
def fnatan2(y,x):
    if x == 0.0:
        if y == 0.0:
            return 99.9;
        elif y > 0.0:
            return pi / 2;
        else:
            return - pi / 2;
    else:
        if x > 0.0:
            return math.atan(float(y/x));
        elif x < 0.0:
            if y >= 0.0:
                return math.atan(float(y/x)) + pi;
            else:
                return math.atan(float(y/x)) - pi;

def fnatan2d(y,x):
    return radeg*fnatan2((y/radeg),(x/radeg));

# Normalize an angle between 0 and 360 degrees
# Use Double Precision, if possible
def fnrev(x):
    rv = x - math.trunc(float(x/360.0))*360.0;
    if rv < 0.0:
        return rv + 360;
    return rv;

# Cube Root (needed for parabolic orbits)
def fncbrt(x):
    if x > 0:
        return math.exp( math.log(x)/3 );
    elif x < 0:
        return - math.exp( math.log(-x)/3 );
    else:
        return 0;

# spherical to rectangular coordinates
def sphe2rect(r,RA,Decl):
    x = r * fncosd(RA) * fncosd(Decl); 
    y = r * fnsind(RA) * fncosd(Decl);
    z = r * fnsind(Decl);
    return x,y,z;

# rectangular to spherical coordinates
def rect2sphe(x,y,z):
    r    = math.sqrt( x*x + y*y + z*z );
    RA   = fnatan2d( y, x );
    Decl = fnatan2d( z, math.sqrt( x*x + y*y ) );
    return r,RA,Decl;

# ecliptic to equatorial coordinates
def eclip2equat(x,y,z,o):
    xequat = x;
    yequat = y * fncosd(o) - z * fnsind(o);
    zequat = y * fnsind(o) + z * fncosd(o);
    return xequat, yequat, zequat;

# equatorial to ecliptic coordinates
def equat2eclip(x,y,z,o):
    xeclip = x;
    yeclip = y * fncosd(o) + z * fnsind(o);
    zeclip = -y * fnsind(o) + z * fncosd(o);
    return xeclip, yeclip, zeclip;

# function for converting degrees, minutes, seconds (DMS) coordinates to decimal degrees (DD)
def dms_to_dd(d, m, s):
    dd = float(d) + float(m)/60 + float(s)/3600;
    return dd

# function for converting decimal degrees (DD) to degrees, minutes, seconds (DMS) coordinates
def dd_to_dms(dd):
    mnt,sec = divmod(dd*3600,60)
    deg,mnt = divmod(mnt,60)
    deg = int(deg);
    mnt = int(mnt);
    return deg,mnt,sec;

def sol(Y, M, D, UT, LON, lat):
    # date calculation
    d = 367*Y - (7*(Y + ((M+9)//12)))//4 + (275*M)//9 + D - 730530;

    # The Sun's position
    #=========================
    w = 282.9404 + 4.70935E-5   * d;         #    (longitude of perihelion)
    a = 1.000000;                            #    (mean distance, a.u.)
    e = 0.016709 - 1.151E-9 * d;             #    (eccentricity)
    M = 356.0470 + 0.9856002585 * d;         #    (mean anomaly)

    # obliquity of the ecliptic
    oblecl = 23.4393 - 3.563E-7 * d;

    # mean anomaly reduction to normal degrees, M
    M = fnrev(M);

    # Sun's mean longitude, L
    L = w + M;
    L = fnrev(L);

    # eccentric anomaly, E
    E = M + radeg * e * fnsind(M) * (1 + e * fncosd(M));

    # Sun's rectangular coordinates in the plane of the ecliptic, where the X axis points towards the perihelion
    x = fncosd(E) - e;
    y = fnsind(E) * math.sqrt(1 - e*e);

    # Convert to distance and true anomaly
    r = math.sqrt(x*x + y*y);
    v = fnatan2d(y,x);

    # longitude of the Sun at time
    lon = fnrev(v + w);

    #Sun's ecliptic rectangular coordinates
    x = r * fncosd(lon);
    y = r * fnsind(lon);
    z = 0.0;

    #rotate these to equatorial coordinates
    xequat, yequat, zequat = eclip2equat(x, y, z, oblecl);

    # rectangular to spherical
    r,RA,Decl = rect2sphe(xequat, yequat, zequat);
    # convert RA to hours
    RA = RA/15;

    # Sidereal Time at the Greenwich
    GMST0 = ( L + 180 ) / 15;

    # Sidereal Time for the time meridian of Central Europe at UTC 0.00 (noe lon = 15°)
    # ask UT
    # ask longitude         
    # Longitude to hours => /15     My calibration
    SIDTIME = GMST0 + UT + LON/15 - 0.0002;
    if SIDTIME > 24:
        SIDTIME = SIDTIME - 24;
    elif SIDTIME < 0:
        SIDTIME = SIDTIME + 24;
    
    #The Hour Angle in hours
    HA = SIDTIME - RA;
    #The Hour Angle in degrees
    HA = HA * 15;

    # Convert Sun's HA and Decl to a rectangular (x,y,z) coordinate system 
    x = fncosd(HA) * fncosd(Decl);
    y = fnsind(HA) * fncosd(Decl);
    z = fnsind(Decl);

    # Rotate this x,y,z system along an axis going east-west, L is +15deg, lat +60deg
    # ask latitude
    xhor = x * fnsind(lat) - z * fncosd(lat);
    yhor = y;
    zhor = x * fncosd(lat) + z * fnsind(lat);

    #azimuth, altitude
    azimuth  = fnatan2d( yhor, xhor ) + 180;
    altitude = fnatan2d( zhor, math.sqrt(xhor*xhor+yhor*yhor) );
    return azimuth, altitude;


def cobj(Y, M, D, UT, RA, Decl, LON, lat):
    # date calculation
    d = 367*Y - (7*(Y + ((M+9)//12)))//4 + (275*M)//9 + D - 730530;

    # obliquity of the ecliptic
    oblecl = 23.4393 - 3.563E-7 * d;

    # Sidereal Time for the time meridian of Central Europe at UTC 0.00 (noe lon = 15°)
    # first in degrees, fit 0 to 360 and lastly divide with 15 to hours   My calibration
    SIDTIME = fnrev(100.4606184 + 0.9856473662862 * d + UT*15 + LON) / 15 - 0.0984;
    #
    if SIDTIME > 24:
        SIDTIME = SIDTIME - 24;
    elif SIDTIME < 0:
        SIDTIME = SIDTIME + 24;

    #The Hour Angle in hours
    HA = SIDTIME - RA;
    #The Hour Angle in degrees
    HA = HA * 15;
    
    # Convert HA and Decl to a rectangular (x,y,z) coordinate system 
    x = fncosd(HA) * fncosd(Decl);
    y = fnsind(HA) * fncosd(Decl);
    z = fnsind(Decl);

    # Rotate this x,y,z system along an axis going east-west, L is +15deg, lat +60deg
    # ask latitude
    xhor = x * fnsind(lat) - z * fncosd(lat);
    yhor = y;
    zhor = x * fncosd(lat) + z * fnsind(lat);

    #azimuth, altitude
    azimuth  = fnatan2d( yhor, xhor ) + 180;
    altitude = fnatan2d( zhor, math.sqrt(xhor*xhor+yhor*yhor) );
    return azimuth, altitude;

def input_data(c, x, y, z):         # data inquiry part
    global start_date;
    # Observatory data
    df2 = pd.read_csv('Observatory.csv');
    obs = df2.set_index('Name').index.get_loc(x)
    name = x
    obslimit = df2['obslimit'][obs]

    # Observatory hours based on location
    start_date = datetime.strptime(str(z), "%Y-%m-%d")
    start_date = start_date + dt.timedelta(hours=int(df2['starthour'][obs])) + dt.timedelta(minutes=int(df2['startminute'][obs]));
    end_date = start_date + dt.timedelta(hours=int(df2['durahour'][obs])) + dt.timedelta(minutes=int(df2['duraminute'][obs]));
    
    # Get longitude and latitude
    LAT = df2['Lat'][obs];
    LON = df2['Lon'][obs];

    if c == 'MessierObjects':
        df = pd.read_csv('MessierObjects.csv');
    elif c == 'CaldwellObjects':
        df = pd.read_csv('CaldwellObjects.csv');
    else:
        df = pd.read_csv('VariableStarObjects.csv');

    obj = df.set_index('Object').index.get_loc(y)
    
    id_mo = df['Object'][obj];   # Object identification

    # get object position
    RA = dms_to_dd(int(df['RA'][obj][0:2]),int(df['RA'][obj][3:5]),float(df['RA'][obj][6:]));
    Decl = dms_to_dd(int(df['Dec'][obj][1:3]),int(df['Dec'][obj][4:6]),float(df['Dec'][obj][7:]));
    
    # Southern hemisphere correction
    if df['Dec'][obj][0] == '-':
        Decl = -1.0 * Decl;    
    return RA, Decl, start_date, end_date, obslimit, id_mo, name, LAT, LON;

def loop(start_date, end_date, obslimit, info_text, RA, Decl, LON, LAT):    #going trough data
    delta_time = dt.timedelta(minutes=1);      # CHANGE this for loop steps

    # loop through time 
    while start_date <= end_date:
        Y = start_date.year;
        M = start_date.month;
        D = start_date.day;
        UT = start_date.hour + start_date.minute / 60;

        # azimuthal position on the observation place    
        azimuth, altitude = cobj(Y, M, D, UT, RA, Decl, LON, LAT);

        # check sun's position, append appropriate data
        sol_azimuth, sol_altitude = sol(Y, M, D, UT, LON, LAT);

        if sol_altitude < -6 and altitude > obslimit:
            # Note rising time
            if len(arvo) == 0:
                text = 'Nousee ' + datetime.time(start_date).strftime("%H:%M") + ' korkeudella ' + f'{altitude:.2f}°';
                info_text.append(text);
                
            arvo.append(altitude);
            pvm.append(start_date);
        else:
            altitude = altitude;

        # increase the date value and get back to make a new
        start_date += delta_time;

    # print transit and setting times
    if len(arvo) > 1:
        max_index = arvo.index(max(arvo))
        if max_index == 0 or max_index + 1 == len(arvo):
            text = '      --     Ei huippua      --          ';
        else:
            text = ' Huippu ' + datetime.time(pvm[max_index]).strftime("%H:%M") + ' korkeudella ' + f'{max(arvo):.2f}°';
        info_text.append(text);
        text = ' Laskee ' + datetime.time(pvm[-1]).strftime("%H:%M") + ' korkeudella ' + f'{arvo[-1]:.2f}°';
        info_text.append(text);
    return arvo, pvm, info_text

def figure_mo(arvo, pvm, obslimit, info_text, id_mo, RA, Decl, LON, LAT):
    global fig_name, fs
    # Only on subplot, screen division (rows, columns), screen size in inches
    fig, ax = plt.subplots(figsize=(16,10));

    # data frame
    df3 = pd.DataFrame(arvo, columns = ['arvo'], index = pvm);
    
    # graphic plotting
    ax.plot(df3);
    ax.set_title('Tähtitaivaan kohde ' + id_mo + ' RA / Dec ' + "{:.4f}".format(RA) + ', ' + "{:.4f}".format(Decl) + ' paikassa: ' + name + ' Lat / Lon '
                 + "{:.4f}".format(LAT) + '°' + "{:.4f}".format(LON) + '° pvm: ' + start_date.strftime("%Y-%m-%d"));
    ax.set_xlabel('aika (hh-mm pp)');                      # labels
    ax.set_ylabel('Kohteen korkeus (°)');
    ax.grid(linestyle='--', linewidth=1);           # grid specifics
    ax.format_xdata = mdates.DateFormatter('%m');   # date formatting to x-axis
    plt.setp(plt.xticks()[1], rotation=30, ha='right')
    fig.autofmt_xdate();
    ax.set_ylim(obslimit-10, 100, 10);                        # set y-axis
    ax.spines['bottom'].set_position(('data', obslimit)); # lift x-axis to critical altitude
    if len(info_text) > 0:
        ax.text(.01, 0, info_text[0]+info_text[1]+info_text[2],
            horizontalalignment='left',
            verticalalignment='bottom',
            transform=ax.transAxes,
            bbox=dict(boxstyle="round",
            ec=(1., 0, 0),
            fc=(1., 0.6, 0.05),
            ));
    sig = "NSO ver 1.2 JB 2020";
    ax.text(1,0.,sig,horizontalalignment='right', verticalalignment='bottom', transform=ax.transAxes);

    # Save the figure
    if fs:
        fig_name = id_mo + '_' + start_date.strftime("%Y-%m-%d") + '_' + name + '.jpg'
    else:
        fig_name = 'sivutallennus.jpg'
    fig.savefig(fig_name)
    
def find_day_curve(c, x, y, z):
    # Main program:
    # initial parameters
    global pi, radeg, fii, arvo, pvm, info_text, name, end_date, RA, Decl
    global LON, LAT, obslimit, id_mo
    pi = 3.14159265359;
    radeg = 180/pi;
    fii = 90;
    arvo = []; # altitudes
    pvm = [];  # dates
    info_text = [];

    # Get Input data
    RA, Decl, start_date, end_date, obslimit, id_mo, name, LAT, LON = input_data(c, x, y, z);

    # Find all the data with the object during the day
    loop(start_date, end_date, obslimit, info_text, RA, Decl, LON, LAT);

    # Plot the image
    figure_mo(arvo, pvm, obslimit, info_text, id_mo, RA, Decl, LON, LAT);

class Ui_MainWindow(object):
    global currentYear, currentMonth, currentDay, x, y, z, fig_name, fs, df1, d
    currentMonth = datetime.now().month
    currentYear = datetime.now().year
    currentDay = datetime.now().day
    fs = False
    fig_name = "Mallikuva.jpg"
    d = {'Catalog':['MessierObjects','CaldwellObjects','VariableStarObjects']}
    #df1 = pd.DataFrame(data=d)
    
    def setupUi(self, MainWindow):
        MainWindow.setObjectName("MainWindow")
        MainWindow.resize(1152, 1020)
        self.centralwidget = QtWidgets.QWidget(MainWindow)
        self.centralwidget.setObjectName("centralwidget")
        # Calendar
        self.calendar = QtWidgets.QCalendarWidget(self.centralwidget)
        self.calendar.setGeometry(576, 0, 500, 300)
        font = QtGui.QFont()
        font.setPointSize(12)
        self.calendar.setFont(font)
        self.calendar.setGridVisible(True)
        self.calendar.setMinimumDate(QDate(currentYear-20, 1, 1))
        self.calendar.setMaximumDate(QDate(currentYear+20, 12, 31))
        self.calendar.setSelectedDate(QDate(currentYear, currentMonth, currentDay))
        self.calendar.setObjectName("calendar")
        # Catalogue
        self.comboZ = QtWidgets.QComboBox(self.centralwidget)
        self.comboZ.setGeometry(QtCore.QRect(126, 30, 300, 30))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.comboZ.setFont(font)
        self.comboZ.setObjectName("comboZ")
        df1 = pd.DataFrame(data=d)
        self.comboZ.addItems(df1['Catalog'])
        self.comboZ.activated[str].connect(self.comboBoxChanged)
        # Observatory
        self.comboX = QtWidgets.QComboBox(self.centralwidget)
        self.comboX.setGeometry(QtCore.QRect(126, 70, 300, 30))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.comboX.setFont(font)
        self.comboX.setObjectName("comboX")
        df2 = pd.read_csv('Observatory.csv')
        self.comboX.addItems(df2['Name'])
        # Objects
        self.comboY = QtWidgets.QComboBox(self.centralwidget)
        self.comboY.setGeometry(QtCore.QRect(126, 110, 300, 30))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.comboY.setFont(font)
        self.comboY.setObjectName("comboY")
        df = pd.read_csv(self.comboZ.currentText() + '.csv')
        self.comboY.addItems(df['Object'])
        # Button
        self.submit = QtWidgets.QPushButton(self.centralwidget)
        self.submit.setGeometry(QtCore.QRect(176, 250, 220, 40))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.submit.setFont(font)
        self.submit.setObjectName("submit")
        # Label
        self.label = QtWidgets.QLabel(self.centralwidget)
        self.label.setGeometry(QtCore.QRect(150, 150, 300, 30))
        font = QtGui.QFont()
        font.setPointSize(16)
        self.label.setFont(font)
        self.label.setObjectName("label")
        # Check Box
        self.checkbox = QtWidgets.QCheckBox(self.centralwidget)
        self.checkbox.setGeometry(QtCore.QRect(126, 150, 21, 30))
        self.checkbox.setObjectName("checkbox")
        self.checkbox.stateChanged.connect(self.checkboxChangedAction)
        # Photo        
        self.photo = QtWidgets.QLabel(self.centralwidget)
        self.photo.setGeometry(QtCore.QRect(0, 300, 1152, 720))
        self.photo.setText("")
        self.photo.setPixmap(QtGui.QPixmap(fig_name))
        self.photo.setScaledContents(True)
        self.photo.setObjectName("photo")
        # Main Window
        MainWindow.setCentralWidget(self.centralwidget)
        self.menubar = QtWidgets.QMenuBar(MainWindow)
        self.menubar.setGeometry(QtCore.QRect(0, 0, 800, 21))
        self.menubar.setObjectName("menubar")
        MainWindow.setMenuBar(self.menubar)
        self.statusbar = QtWidgets.QStatusBar(MainWindow)
        self.statusbar.setObjectName("statusbar")
        MainWindow.setStatusBar(self.statusbar)

        self.submit.clicked.connect(self.pressed)

        self.retranslateUi(MainWindow)
        QtCore.QMetaObject.connectSlotsByName(MainWindow)

    def retranslateUi(self, MainWindow):
        _translate = QtCore.QCoreApplication.translate
        MainWindow.setWindowTitle(_translate("MainWindow", "Kohteiden näkyvyys Slooh:n observatorioissa"))
        self.submit.setText(_translate("MainWindow", "Hae tiedot"))
        self.label.setText(_translate("MainWindow", "Kuvaa ei tallenneta"))
        
    def comboBoxChanged(self):
        self.comboY.clear()
        df = pd.read_csv(self.comboZ.currentText() + '.csv')
        self.comboY.addItems(df['Object'])

    def checkboxChangedAction(self, state):
        global fs
        if (QtCore.Qt.Checked == state):
            self.label.setText("Kuva tallennetaan")
            fs = True
        else:
            self.label.setText("Kuvaa ei tallenneta")
            fs = False

    def pressed(self):
        c = self.comboZ.currentText()
        x = self.comboX.currentText()
        y = self.comboY.currentText()
        z = self.calendar.selectedDate()
        z = z.toPyDate()

        find_day_curve(c, x, y, z)
        self.photo.setPixmap(QtGui.QPixmap(fig_name))

if __name__ == "__main__":
    import sys
    app = QtWidgets.QApplication(sys.argv)
    MainWindow = QtWidgets.QMainWindow()
    ui = Ui_MainWindow()
    ui.setupUi(MainWindow)
    MainWindow.show()
    sys.exit(app.exec_())