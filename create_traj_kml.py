#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# pylint: disable=trailing-whitespace,bad-whitespace,invalid-name
# pylint: disable=anomalous-backslash-in-string,bad-continuation
# pylint: disable=multiple-statements,redefined-outer-name,global-statement
"""
FILE: create_traj_kml.py
DATE: 16 Sep 2020
AUTH: G. E. Deschaines
PROG: Creates trajectory KML file for display with Google Earth.
DESC: Given traj3D version and name of input namelist file,
      this program reads associated traj3D case input and
      output files to create KML trajectory file containing
      following sections.
      
         HEADER:
             - Open kml tag
             - Open Document tag
             - Label, Line and Poly Style tags
             - Open Folder tag
             - Open Placemark for trajectory LineString
             - Open coordinates tag
         COORDS:
             - lines of CSV for longitude, latitude & height 
         FOOTER:
             - Close coordinates tag
             - Close Placemark for trajectory LineString
             - Placemark for trajectory launch point
             - Placemark for trajectory maximum altitude point
             - Placemark for trajectory impact point
             - Placemark for great circle ground path
             - Close Folder tag
             - Close Document tag
             - CLose kml tag
     
     The following traj3D case namelist input and output files
     are assumed to exist for the given version 'vers' and
     namelist input .txt file with 'name'.
     
       - namelist input:        ./txt/traj3D{vers}_{name}.txt
       - standard output:       ./out/traj3D{vers}_{name}.out
       - geodesic data output:  ./dat/traj3D{vers}_geod_{name}.dat
       
     The created KML file will be:
         
       ./kml/traj3D{vers}_traj_{name}.kml

     using following trajectory KML header template, folder
     description and footer template files.

       - ./kml/traj_header.kml
       - ./kml/traj_folder_desc_{name}.kml
       - ./kml/traj_footer.kml
"""
import sys
from math import pi, atan2, sin, cos

try:
    import numpy as np
except ImportError:
    print("* Error: NumPy required.")
    print("         Suggest installing the SciPy stack.")
    sys.exit()
    
# Units of Measure Conversion Constants

RPD  = pi/180.0   # radians per degree
DPR  = 180.0/pi   # degrees per radian
MPFT = 0.3048     # meters per foot
MPNM = 1853.184   # meters per nautical mile

def aviation_bearing(latA,lonA,latB,lonB):
    """ Assumes geocentric latitudes and longitudes are given in radians.
    """
    dlon = lonB - lonA
    azd = atan2(sin(dlon)*cos(latB),\
                cos(latA)*sin(latB)-sin(latA)*cos(latB)*cos(dlon))*DPR
    azd = np.mod(azd+360.0, 360.0)
    return azd

if __name__ == "__main__":
    
    # Check arguments for version and input namelist name.
    
    name = ''
    if len(sys.argv) > 2:
        vers = sys.argv[1]
        name = sys.argv[2]
    else:
        print("usage: create_kml [0|1] name")
        print("where: [0|1] - traj3D executable version.") 
        print("       name  - name of input namelist .txt file.")
        exit()

    xname = 'traj3D' + vers
    
    # Open traj3D namelist input text file.
  
    ipath = './txt/' + name + '.txt'
    try:
        txt_file = open(ipath, 'r')
    except OSError:
        print('open error for', ipath)
        exit()
      
    # Read case description from namelist file.
        
    lines = txt_file.readlines()
    txt_file.close()
    
    sdesc = lines[0].lstrip().rstrip()
  
    # Read TFINAL value to determine trajectory mode (Observed or Inertial).
    
    tfinal_idx = lines[2].find('TFINAL')
    if tfinal_idx > -1:
        equal_idx = lines[2].find('=', tfinal_idx)+1
        comma_idx = lines[2].find(',', tfinal_idx)
        tfinal = float(lines[2][equal_idx:comma_idx])
        if tfinal > -1.0:
            if tfinal > 0.0:
                smode = 'Inertial'
                sevnt = 'Impact Time - __TOFHRS__ hours'
                color = 'Red'
            else:
                smode = 'Relative'
                sevnt = 'Launch and Impact Time'
                color = 'Purple'
        else:
            smode = 'Observed'
            sevnt = 'Launch Time + __TOFHRS__ hours'
            color = 'Green'
    else:
        smode = 'Observed'
        sevnt = 'Launch Time + __TOFHRS__ hours'
        color = 'Green'

    # Read ROTOPT value to determine trajectory color (Observed, Inertial or Non-rotating).
    
    rotopt_idx = lines[6].find('ROTOPT')
    if rotopt_idx > -1:
        equal_idx = lines[6].find('=', rotopt_idx)+1
        comma_idx = lines[6].find(',', rotopt_idx)
        rotopt = float(lines[6][equal_idx:comma_idx])
        if rotopt == 0.0:
            smode = 'Non-rotating'
            sevnt = 'Launch and Impact Time'
            color = 'Blue'

    # Read GAMREL value.
    
    gamrel_idx = lines[3].find('GAMREL')
    if gamrel_idx > -1:
        equal_idx = lines[3].find('=', gamrel_idx)+1
        comma_idx = lines[3].find(',', gamrel_idx)
        el0d = float(lines[3][equal_idx:comma_idx])
    else:
        el0d = 0.0
    
    # Open traj3D standard output file.
    
    opath = './out/' + xname + '_' + name + '.out'
    try:
        out_file = open(opath, 'r')
    except OSError:
        print('open error for', opath)
        exit()
        
    # Read trajectory final state values from standard output file.
          
    lines = out_file.readlines()
    out_file.close()
    
    [stofsec, svrfps, saltft, slatdeg, _] = lines[-6][1:].split(None,4)
    [  sgsec,  smach, srngnm, slondeg, _] = lines[-5][1:].split(None,4)
    [          sppsf, sazdeg,  sgamma, _] = lines[-4][1:].split(None,3)
    
    tofsec = float(stofsec[0:])  # time of flight (seconds)
    tofhrs = tofsec/3600.0       # time of flight (hours)
    vrfps  = float(svrfps[0:])   # relative velocity (feet per second)
    vrmps  = vrfps*MPFT          # relative velocity (meters per second)
    altft  = float(saltft[0:])   # altitude (feet)
    altm   = altft*MPFT          # altitude (meters)
    rngnm  = float(srngnm[0:])   # ground range (nautical miles)
    rngkm  = rngnm*MPNM/1000.0   # ground range (kilometers)
    azdeg  = float(sazdeg[0:])   # flight path heading (degrees)
    gamma  = float(sgamma[0:])   # flight path pitch angle (degrees)

    period_idx = slatdeg.find('.')
    latdeg = float(slatdeg[0:period_idx+3])  # latitude (degrees)  
    period_idx = slondeg.find('.')
    londeg = float(slondeg[0:period_idx+3])  # longitude (degrees)
    if londeg >  180.0: londeg = londeg - 360.0
    if londeg < -180.0: londeg = londeg + 360.0
    
    # Open traj3D geodesic data output file.
    
    dpath = './dat/' + xname + '_geod_' + name + '.dat'
    try:
        dat_file = open(dpath, 'r')
    except OSError:
        print('open error for', dpath)
        exit()

    # Read lines of geodesic data and store in lists
    # of latitude, longitude and height values.
    
    lines = dat_file.readlines()
    dat_file.close()
    
    n = len(lines)
    flon = np.ndarray(n, dtype=np.float)
    flat = np.ndarray(n, dtype=np.float)
    ihgt = np.ndarray(n, dtype=np.int)
    for i in range(n):
        [slon, slat, shgt] = lines[i].split(',')
        flon[i] = float(slon)
        flat[i] = float(slat)
        ihgt[i] = int(shgt)
    m = np.argmax(ihgt)
  
    az0d = aviation_bearing(flat[0]*RPD,flon[0]*RPD,flat[1]*RPD,flon[1]*RPD)
    azMd = aviation_bearing(flat[0]*RPD,flon[0]*RPD,flat[m]*RPD,flon[m]*RPD)
    azNd = aviation_bearing(flat[-1]*RPD,flon[-1]*RPD,flat[-2]*RPD,flon[-2]*RPD)
    
    # Read trajectory header KML file.
    
    kml_file = open('./kml/traj_folder_desc_' + name + '.kml', 'r')
    folder_kml = kml_file.read()
    kml_file.close()
    
    # Read trajectory header KML file.
    
    kml_file = open('./kml/traj_header.kml', 'r')
    header_kml = kml_file.read()
    kml_file.close()
    
    # Read trajectory footer KML file.
  
    kml_file = open('./kml/traj_footer.kml', 'r')
    footer_kml = kml_file.read()
    kml_file.close()
  
    # Replace header and footer key values with those read
    # from traj3D input namelist text file, standard output
    # and geodesic data output files.
  
    srngkm  = "%.3f" % rngkm
    stofhrs = "%.4f" % tofhrs
    slon0   = "%.6f" % flon[0]
    slat0   = "%.6f" % flat[0]
    shgt0   = "%d"   % ihgt[0]
    saz0d   = "%.2f" % az0d
    sel0d   = "%.2f" % el0d
    slonM   = "%.6f" % flon[m]
    slatM   = "%.6f" % flat[m]
    shgtM   = "%d"   % ihgt[m]
    shgtMkm = "%.3f" % (float(ihgt[m])/1000.0)
    sazMd   = "%.2f" % azMd
    slonN   = "%.6f" % flon[n-1]
    slatN   = "%.6f" % flat[n-1]
    shgtN   = "%d"   % ihgt[n-1]
    sazNd   = "%.2f" % azNd
  
    header_kml = header_kml.replace('__FLDRDESC__', folder_kml)
    
    header_kml = header_kml.replace('__VERS__',   vers)
    header_kml = header_kml.replace('__DESC__',   sdesc)
    header_kml = header_kml.replace('__MODE__',   smode)
    header_kml = header_kml.replace('__RNGKM__',  srngkm)
    header_kml = header_kml.replace('__TOFHRS__', stofhrs)
    header_kml = header_kml.replace('__LON0__',   slon0)
    header_kml = header_kml.replace('__LAT0__',   slat0)
    header_kml = header_kml.replace('__AZ0__',    saz0d)
    header_kml = header_kml.replace('__EL0__',    sel0d)
    header_kml = header_kml.replace('__COLOR__',  color)
    
    footer_kml = footer_kml.replace('__DESC__',   sdesc)
    footer_kml = footer_kml.replace('__MODE__',   smode)
    footer_kml = footer_kml.replace('__EVNT__',   sevnt)
    footer_kml = footer_kml.replace('__TOFHRS__', stofhrs)
    footer_kml = footer_kml.replace('__LON0__',   slon0)
    footer_kml = footer_kml.replace('__LAT0__',   slat0)
    footer_kml = footer_kml.replace('__HGT0__',   shgt0)
    footer_kml = footer_kml.replace('__AZ0__',    saz0d)
    
    footer_kml = footer_kml.replace('__HGTMKM__', shgtMkm)
    footer_kml = footer_kml.replace('__LONM__',   slonM)
    footer_kml = footer_kml.replace('__LATM__',   slatM)
    footer_kml = footer_kml.replace('__HGTM__',   shgtM)
    footer_kml = footer_kml.replace('__AZM__',    sazMd)
    
    footer_kml = footer_kml.replace('__LONN__',   slonN)
    footer_kml = footer_kml.replace('__LATN__',   slatN)
    footer_kml = footer_kml.replace('__HGTN__',   shgtN)
    footer_kml = footer_kml.replace('__AZN__',    sazNd)
    footer_kml = footer_kml.replace('__COLOR__',  color)
   
    # Create trajectory KML file.
    
    kpath = './kml/' + xname + '_traj_' + name + '.kml'
    kml_file = open(kpath, 'w')
    # Write header section
    kml_file.write(header_kml)
    # Write coords section.
    for i in range(n):
        line = "        %.6f,%.6f,%d\n" % (flon[i],flat[i],ihgt[i])
        kml_file.write(line)
    # Write footer section.
    kml_file.write(footer_kml)
    kml_file.close()
  
