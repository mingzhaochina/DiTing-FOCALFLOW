#! /usr/bin/python
"""stereonet: Prepares structural data for plotting on a stereonet
using GMT's psxy.

Provided under an MIT-style license.

Joe Kington <jkington@wisc.edu>
10/2008
"""
__license__   = "MIT License <http://http://www.opensource.org/licenses/mit-license.php>"
__copyright__ = "2008, Joe Kington IV"
__author__    = "Joe Kington <jkington@wisc.edu>"


#-------------------------------------------------------------------------
#--Try to import modules--------------------------------------------------
#-------------------------------------------------------------------------
import sys, copy
from math import radians,degrees,cos,sin,asin,atan2,atan,sqrt
try:
    from optparse import OptionParser, OptionGroup
except ImportError:
    message = """%s requires the optparse module, included by 
    default in python 2.3 or greater.  If you are using an older version of 
    python, either install these modules or upgrade your python interpreter.  
    You are currently using python ver. %s""" % (sys.argv[0], sys.version)
    message = message.replace('   ','') #For better display
    sys.exit(message)



#----------------------------------------------------------------------------------
#--Option Handling-----------------------------------------------------------------
#----------------------------------------------------------------------------------

def main(argv):
    """Parse arguments and execute routine based on them"""


    #---------------------------------------------------------------
    #--Set up Options-----------------------------------------------
    #---------------------------------------------------------------
    usage = "usage: %prog [infile] [outfile] [-p|-P|-L|-R] [-I] [-H] [-C] [-i] [-:]"
    description = """Formats structural data for plotting on a 
    stereonet using GMT's psxy. Outputs tab delimited x,y pairs 
    representing a feature based on a strike and dip, plunge and
    bearing, or rake along a plane. To use this program to plot
    data on a stereonet, use psbasemap with the -JA0/0 option, and
    set the extent to be -90/90/-90/90 (-Rd works fine). If no infile
    or outfile are specified, the program reads from stdin and writes 
    to stdout. Infile should contain one structural measurement per 
    line. Lines in an infile starting with '#' will be treated as 
    comments. Strikes may be in azimuth or quadrant form.  If no dip
    direction is given, the right hand rule is assumed. See the 
    examples given in the options for details on how data may be
    formatted.  Optionally, an 'H' and the S/D of an originally 
    horizontal plane (e.g. bedding) may be appended after a 
    measurement to rotate the measurement so that the given plane
    is horizontal. (e.g. if you have a paleocurrent measurement of
    20/315N that has been tilted so that bedding in the area dips
    at 076/89, "20/315N H 076/89" would rotate the output point so
    that bedding is horizontal and the paleocurrent measurement is
    in its original orientation.)"""
    description = description.replace('    ','') #Strip out spaces for better display

    examples = """To plot the plane 052/36SW and its pole on a stereonet:
    psbasemap -JA0/0/6i -Rd -B -K > output.ps                              
    echo 052/36SW | stereonet --planes | psxy -JA0/0/6i -Rd -O -K -W2p/red >> output.ps
    echo 052/36SW | stereonet --poles | psxy -JA0/0/6i -Rd -O -G10p/blue >>output.ps"""

    parser = OptionParser(usage=usage, description=description)
    
    parser.add_option("-p", "--planes", dest="PlotType",\
            help="Output planes as lines given a S/D measurement (e.g 330/42W, 150/42, 75/10SW, N30E/20NW). This is the default.", \
            action="store_const", const='Planes')
    parser.add_option("-P", "--poles", dest="PlotType",\
            help="Output poles to planes as points given a S/D measurement (e.g 330/42W, 150/42, 75/10SW).", \
            action="store_const", const='Poles')
    parser.add_option("-L", "--lines", dest="PlotType",\
            help="Output lines as points given a P/B measurement (e.g 25/200N, 20/025, 75/225SW, 23/E30N SW).", \
            action="store_const", const='Lines')
    parser.add_option("-R", "--rakes", dest="PlotType",\
            help="Output rakes along a plane as points given a S/D and rake measurement (e.g 330/42W 22N, 150/42 30, 075/38SE 32SW).", \
            action="store_const", const='Rakes')
    parser.add_option("-H", "--horizontal", dest="Flatten",\
            help="Rotate all measurements in input so that the specified plane is horizontal. FLATTEN should be in the form of a S/D measurement, e.g. 234/64NW.", \
            action="store", type="string")
    parser.add_option("-I", "--invert", dest="Invert",\
            help="Convert each long,lat pair in input to a plunge/bearing or the strike/dip of a plane defined as the pole to the long,lat pair. Useful for converting the output an analysis (e.g. 'fitcircle') back to a more readable format.",\
            action="store", choices=('Lines','Line','line','lines','Planes','Plane','plane','planes','Pole','pole','Poles','poles') )
    parser.add_option("-C", "--clean", dest="Clean",\
            help="'Clean-up' input data. Outputs S/D, P/B, or rakes in azimuth format following the right hand rule. Use --planes, --lines, or --rakes to set the type of measurement. Defaults to planes.", \
            action="store_true")
    parser.add_option("-:", "--reverse_xy", dest="ReverseXY",\
            help="Output y,x pairs instead of x,y pairs or expect y,x from input if used with -I. This does not change the format of S/D's, P/B's or rakes in either input or output, only the ordering of longitude,latitude pairs.",\
            action="store_true")
    parser.add_option("--parse", dest="Parse",\
            help="Parse the last word of each input line unchanged",\
            action="store_true")
    parser.add_option("-i", "--increment", dest="inc", \
            help="Increment to insert vertices at along a line representing a plane. Default: 10 degrees", \
            action="store", type="int")

    parser.set_defaults(PlotType="Planes", inc=10, ReverseXY=False)

    #Bit of a hack to add examples.  Adds an empty option group with them.
    #Need to write a new formatter that leaves in newlines in some cases
    #   Too much trouble, really... Would need to extend IndentedHelpFormatter 
    #   Or just go back to getopt... Very verbose, either way...
    examples = OptionGroup(parser, 'Examples', description=examples)
    parser.add_option_group(examples)


    #---------------------------------------------------------------------
    #--Parse Options------------------------------------------------------
    #---------------------------------------------------------------------
    (options, args) = parser.parse_args(args=argv[1:])
    try:
        #How many files are we working with?
        if len(args) == 0:   #None, read from stdin, write to stdout
            infile = sys.stdin
            outfile = sys.stdout
        elif len(args) == 1: #One, read from file, write to stdout
            infile = file(args[0], 'r')
            outfile = sys.stdout
        elif len(args) ==2:  #Two, read from first, write to second
            infile = file(args[0], 'r')
            outfile = file(args[1], 'w')

        #More, raise an options parser error and print message
        else: parser.error("Only one input file and one output file are allowed")

    except IOError (errno, strerror):
        #If opening a file fails...
        sys.exit("Cannot access file!\nI/O error(%s): %s" % (errno, strerror)) 

    #-----------------------------------------------------------------------
    #--Read input file and output properly formatted data-------------------
    #-----------------------------------------------------------------------
    for line in infile:
        line = line.strip()

        #--Skip Comments and Blank Lines-------------
        if line.startswith('#') or len(line)==0: continue

        #--W.B.: Parse any last word--
        if options.Parse:
            line, parsant = line.split()[:-1], line.split()[-1]
            line = ' '.join(line)

        try:
            #--Which function are we preforming?---------
            if options.Invert:    data = InvertGeographic(line,options) #Invert long,lat to a S/D or P/B
            elif options.Clean:   data = CleanInput(line,options)       #Output data in azimuths following the RHR     
            else:                 data = OutputXY(line,options)         #Output long,lat pairs corresponding to measurements

            if options.Parse: data = data.strip() + ' ' + parsant + '\n'

            #--Write to output---------------------------
            try: outfile.write(data)
            except: sys.exit('Data could not be written to output!')

        #--If the data wasn't properly formatted, print error and continue
        except InputError as message:
            print(sys.stderr, 'Invalid Input: %s\n' % line, message, '\nSkipping this line...')

#---------------------------------------------------------------------------------------    
#--Decide what to do--------------------------------------------------------------------
#---------------------------------------------------------------------------------------    

def OutputXY(input,options):
    """Calculates the strike and dip based on the input 
    string and returns x,y pairs for plotting in a 
    stereographic projection"""

    #--Basic Concepts-----------------------------------
    """A stereonet in <long,lat> coordinates:
                      <0,90>
                       ***
                    *       *
           <-90,0> *         *<90,0>
                   *         *
                    *       *
                       ***
                     <0,-90>                  
       If strike=0, plotting lines, rakes, planes or
       poles to planes is simple.  For a plane, it's
       a line of constant longitude at long=90-dip.
       For a line, it's a point at long=0,lat=90-dip.
       For a rake, it's a point at long=90-dip,
       lat=90-rake.  These points can then be rotated
       to the proper strike. (A rotation matrix around
       the X-axis is much simpler than the trig 
       otherwise necessary!)"""
    #---------------------------------------------------

    #--Options------------------------------------------
    PlotType = options.PlotType.capitalize()
    inc = options.inc
    Flatten = options.Flatten

    #--Is data going to be "flattened"?-----------------
    #   This checks for the line-by-line case denoted
    #   by an H in the infile, e.g. 340/20 H 035/76
    if 'H' in input.upper():
        portions = [item.strip() for item in input.upper().split('H')]
        if len(portions) == 2:  input,Flatten = portions
        else: raise InputError("Too many H's!")
        
    #---------------------------------------------------
    #--Make Data with strike=North----------------------
    #---------------------------------------------------
    header = ''
    if PlotType == 'Planes':
        strike,dip = ParsePlanes(input) #Returns S/D following RHR

        #If strike=north, planes are lines of constant longitude.
        header = '> %s\n' % (input) #Multisegment GMT format, annotated with S/D
        x,y = [],[]
        for lat in range(-90,90+inc,inc):
            x.append(90-dip)
            y.append(float(lat))

        #If data needs to be 'unfolded', rotate the pole, instead of the plane
        #  This prevents problems with length changes if x,y coordinates from
        #  a plane are rotated instead of a point.
        if Flatten: x,y = [-dip], [0]

    elif PlotType == 'Poles':
        strike,dip = ParsePlanes(input) 

        #If strike=north, the pole to a plane will be at lat=0, long=-dip 
        x,y = [-dip], [0]

    elif PlotType == 'Lines':
        #Returns bearing/plunge with bearing at the end the plunge direction is measured from
        strike,dip = ParseLines(input) 

        #For Lines, plot lat=codip, long=0 and rotate later
        x,y = [0], [90-dip]

    elif PlotType == 'Rakes':
        #Returns a negative rake if the rake angle is measured from the "south" end of the plane
        strike,dip,rake = ParseRakes(input)
        if rake>0: coRake = 90-rake
        else:      coRake = 90+rake

        #For Rakes, Lat = 90-rake, Long = 90-dip
        x,y = [90-dip], [coRake]

    else:  #Shouldn't Happen
        sys.exit("Invalid Plot Type: %s (This shouldn't happen!) Programming error!" % PlotType)

    #---------------------------------------------------
    #--Rotate and Format Data---------------------------
    #---------------------------------------------------

    #--Rotate Data to proper strike----------------------
    X,Y = Rotate(x,y,strike)

    #--Do we need to "Flatten" the data?-----------------
    #   i.e. Rotate to horizontal based on another plane
    if Flatten:
        horizStrike,horizDip = ParsePlanes(Flatten)
        #Rotate to horizStrike=north and make horizDip horizontal
        X,Y = Rotate(X,Y,-horizStrike,-horizDip)
        #Unrotate back to the original strike
        X,Y = Rotate(X,Y,horizStrike)

        #If a pole to a plane was rotated, go back and create the plane
        #  This needs to be explained more clearly... Also, it's a bit hackish...
        if PlotType == 'Planes':
            tmpOpts = copy.copy(options)
            tmpOpts.Flatten = False
            tmpOpts.Invert = 'plane'
            StrikeDip = InvertGeographic('%.5f\t%.5f'%(X[0],Y[0]), tmpOpts)
            return OutputXY(StrikeDip, tmpOpts)

    #--Return as a string--------------------------------
    coords = header  #In case it's a line, use GMT multisegment format
    outputFormat = '%.2f\t%.2f\n' #Formatting string for coordinates
 
    for lon,lat in zip(X,Y): 

        #If point is in the upper hemisphere, get the opposite end
        if lon>90: lon-=180; lat=-lat
        if lon<-90: lon+=180; lat=-lat

        #Is -: set? If so, output lat-long, otherwise output long-lat
        if options.ReverseXY:  coords += outputFormat % (lat,lon) 
        else:                  coords += outputFormat % (lon,lat)

    return coords


def InvertGeographic(input,options):
    """Converts a string containing a long, lat pair into a plunge/bearing
    or strike/dip of the plane perpendictular to it."""
    
    #--Split input into long and lat and convert to floats-----------------
    input = input.split()
    if len(input) != 2:
        raise InputError("Too many or two few fields in input. (Expecting X and Y seperated by whitespace.)")

    #Is -: set? If so, expect lat-long instead of long-lat
    if options.ReverseXY: lat,long = input
    else:                 long,lat = input

    #Convert to floats
    try: long,lat = float(long), float(lat)
    except ValueError: raise InputError("Can't convert %s or %s to a number!" % (long,lat) ) 

    #--Sanity check -------------------------------------------------------
    if (abs(lat) > 90) or (abs(long) > 360) or (long < -180):
        raise InputError("(%.1f, %.1f) is not a valid lat, long pair." % (lat, long))

    #--If using 0<long<360, convert to -180<long<180-----------------------
    if long > 180: long -= 360

    #--Make sure it's in the right hemisphere, if not get the opposite end of the line
    if long > 90:    long -= 180; lat=-lat
    elif long < -90: long += 180; lat=-lat
    
    #--Convert back to plunge and bearing----------------------------------
    x,y,z = sph2cart(long,lat)           #Cartesian Coords
    bearing = atan2(z,y)                 #Bearing will be in y-z plane
    plunge = atan( x/sqrt(y**2 + z**2) ) #Plunge is the angle btw the line and the y-z plane
    plunge,bearing = degrees(plunge), degrees(bearing)
    
    #--Rotate so that 0 is north, not east---------------------------------
    bearing = 90-bearing
    if bearing<0: bearing += 360

    #--Calculate S/D of plane to which the measurement is the pole---------
    strike = bearing+90
    dip = 90-plunge
    if strike>360: strike-=360

    #--Return P/B of line or S/D of plane
    outputFormat = '%.2f/%.2f\n'
    if options.Invert.lower() in ['line','lines']:    return outputFormat % (plunge,bearing)
    elif options.Invert.lower() in ['plane','planes','pole','poles']: return outputFormat % (strike,dip)

def CleanInput(input, options):
    """Takes a line with a S/D, P/B, or rake measurement 
    in any acceptable form and returns a measurement in
    azimuth format following the RHR."""

    #--Direction letters for each quad----------------------
    LineDir = {'I':'NE', 'II':'SE', 'III':'SW', 'IV':'NW'} #Same quad as end of line
    DipDir =  {'I':'SE', 'II':'SW', 'III':'NW', 'IV':'NE'} #Dip dir following RHR

    #--Output Format----------------------------------------
    outFormat = '%.0f/%.0f%s'           #E.g. Strike, dip, and direction
    rakeFormat = outFormat + ' %.0f%s'  #Strike/Dip of plane + Rake angle and direction 

    #--Strike/Dip-------------------------------------------
    if options.PlotType in ['Planes', 'Poles']:
        strike,dip = ParsePlanes(input)
        quad = FindQuadrant(strike)
        output = outFormat % (strike,dip, DipDir[quad])

    #--Plunge/Bearing---------------------------------------
    elif options.PlotType == 'Lines':
        bearing,plunge = ParseLines(input)
        quad = FindQuadrant(bearing)
        output = outFormat % (plunge,bearing, LineDir[quad])

    #--Rakes------------------------------------------------
    elif options.PlotType == 'Rakes':
        strike,dip,rake = ParseRakes(input)
        quad = FindQuadrant(strike)
        end = LineDir[quad]
        #Is the rake measured from the non-RHR end of the plane?
        if rake<0:         #rake is negative if it is...
            rake = -rake
            otherend = strike-180
            if otherend<0: otherend+=360
            end = LineDir[FindQuadrant(otherend)]
        dir = DipDir[quad]
        output = rakeFormat % (strike,dip,dir,rake,end)

    else: #Shouldn't Happen...
        sys.exit("Invalid Plot Type: %s (This shouldn't happen!) Programming error!" % options.PlotType)

    return output+'\n'

#---------------------------------------------------------------------------------------    
#--Parsing Functions--------------------------------------------------------------------
#---------------------------------------------------------------------------------------    

def ParsePlanes(measurement):
    """Reads a string containing strikes and 
    dips in the form xxx/yyDIR (e.g. 315/38SW 
    or 78/9S) and returns a strike and a dip 
    tuple following the right hand rule"""

    #--Read input and split into strike and dip---------------
    #   This part is repeated in other functions, need to refactor...
    pair = measurement.split('/')
    if len(pair) != 2:
        raise InputError("Too many or too few /'s!")
    strikestring,dipstring = pair

    #Get Letters from End Indicating Dip Direction
    (dip,dipdir) = GetDir(dipstring)

    #Parse strike portion
    if strikestring[0].isalpha(): #Then it's a quadrant measurement
        strike = ParseQuadrant(strikestring)
    else:
        try: strike = float(strikestring)
        except: raise InputError("Can't convert %s to a number!" % strikestring)


    #--Does this measurement follow the Right Hand Rule?-------

    #Directions that follow RHR for each quadrant
    #Remember that this is the _dip_ direction, not the strike itself
    Directions = {'I':['S','E','SE'], \
                 'II':['S','W','SW'], \
                'III':['N','W','NW'], \
                 'IV':['N','E','NE']}

    quad = FindQuadrant(strike)
    if dipdir not in Directions[quad] and dipdir: #In case of empty string, assume RHR
        #Not RHR, so get the opposite end
        strike = strike-180
        if strike < 0:
            strike = strike+360

    return (strike,dip)

def ParseLines(input):
    """Takes an input string of the form plunge/bearing 
    (e.g. 10/330S, 20/95, 25/95W) and returns the bearing
    and plunge of the measurement as floats."""

    #--Read input and split into plunge and bearing----------
    portions = input.split('/')
    if len(portions) != 2: raise InputError("Too many or too few /'s!")
    plungeStr, bearingStr = portions


    #--Decide what type of measurement we're dealing with----
    
    #Remove trailing whitespace and internal spaces 
    #  (e.g. make "N30E NW\n" --> N30ENW)
    bearingStr = bearingStr.strip().split()
    bearingStr = ''.join(bearingStr)

    #Do we start with a letter? (if so it's a quadrant measurement)
    if bearingStr[0].isalpha():    
        try:  #length of a quadrant measurement must be 3 or 4!
            if bearingStr[2].isalpha(): #e.g. N9E
                bearing = ParseQuadrant(bearingStr[0:3]) 
                plungeDir = bearingStr[3:]
            else: #e.g. N39E
                bearing = ParseQuadrant(bearingStr[0:4]) 
                plungeDir = bearingStr[4:]
        except IndexError:
            raise InputError("%s is not a valid quadrant measurement!" % bearingStr)

    #Not a quadrant, deal with it normally.
    else: (bearing, plungeDir) = GetDir(bearingStr)

    #--Float conversions-------------------------------------
    try: plunge = float(plungeStr)
    except: raise InputError("Can't convert %s to a number!" % plungeStr)

    #--Sanity Check------------------------------------------
    if (bearing>360) or (bearing<0) or (plunge>90) or (plunge<0):
        raise InputError("Bearing not btw. 0 and 360 or plunge not btw. 0 and 90!")

    #--Do we need the opposite end of the bearing?-----------
    if not isSameEnd(bearing,plungeDir):
        #Opposite direction given, get other end
        bearing = bearing-180
        if bearing < 0: bearing += 360

    return (bearing,plunge)

def ParseRakes(input):
    """Takes an input string of the from "strike/dip plunge"
    (e.g. 330/10W 25S) and returns the strike, dip, and rake
    angle of the line represented by the measurement. Returns
    a negative rake angle if the rake is not from the end of
    the plane the would result from following the RHR based
    on the dip of the plane."""

    #--Read input and split into portions--------------------
    #Seperate S/D from rake angle (seperated by whitespace)
    portions = input.strip().split()
    if len(portions)!= 2: #Needs to be changed.  Get last characters instead of splitting by whitespace
        raise InputError("Too many or too few parts of %s seperated by whitespace" % input)
    SD,rakeStr = portions

    #Treat the S/D of plane like any other
    strike,dip = ParsePlanes(SD)

    #Get the rake angle and direction
    rake,rakeDir = GetDir(rakeStr)

    #--Sanity check-----------------------------------------
    if (rake > 90) | (rake < 0):
        raise InputError('Rake angle must be between 0 and 90!')
    
    #--Do we need the opposite end?-------------------------
    if not isSameEnd(strike,rakeDir): rake = -rake
    
    return (strike,dip,rake)

def ParseQuadrant(strikestring):
    """Takes an input string with a quadrant strike 
    measurement (e.g. N40W instead of 320) and returns
    an azimuth strike as a float."""

    #--Seperate into portions--------------------------
    first_dir= strikestring[0].upper() #First letter
    sec_dir = strikestring[-1].upper() #Last letter
    angle = strikestring[1:-1]         #Everything in between

    #--Try converting angle to a float----------------
    try: angle=float(angle)
    except: raise InputError("%s is not a valid strike. Can't convert %s to a number!" % (strikestring, angle) )

    #--Sanity checks----------------------------------
    valid_dir = ['NE','NW','SE','SW','EN','WN','WS','ES']  
    letters = first_dir + sec_dir

    if (angle<0) or (angle>90) or (letters not in valid_dir):
        raise InputError('%s is not a valid quadrant-format strike!' % strikestring)
    
    #--Starting angles for each direction-------------
    start = {'N':0, 'S':180, 'E':90, 'W':270}

    #--Do we need to add or subtract angle from first_dir?
    if letters in ['NE','SW','WN','ES']: 
        strike = start[first_dir] + angle
    else: 
        strike = start[first_dir] - angle
    #--Make positive for NxxW measurements-----------
    if strike<0: strike += 360

    return strike

#---------------------------------------------------------------------------------------    
#--Various Utility Functions------------------------------------------------------------
#---------------------------------------------------------------------------------------    

def Rotate(longs,lats,strike,dip=0):
    """Returns a two lists of rotated x and y coordinates, 
    given a list of x coords, a list of y coords, and a strike
    (in degrees) for the coords to be rotated to. Optionally, a
    dip angle may be specified. (usually to rotate measurements
    back to horizontal)"""

    #--Simplified rotation matrices---------------------
    def XAxisRotate(x,y,z,theta):
        X = x
        Y = y*cos(theta) + z*sin(theta)
        Z = -y*sin(theta) + z*cos(theta)
        return X,Y,Z
    def ZAxisRotate(x,y,z,omega):
        X = x*cos(omega) + y*sin(omega)
        Y = -x*sin(omega) + y*cos(omega)
        Z = z
        return X,Y,Z

    #--Rotate data--------------------------------------
    #-Convert to radians
    theta = radians(strike)
    omega = radians(dip)

    #-Convert to cartesian and rotate
    X,Y = [],[]
    for lon,lat in zip(longs,lats):
        #Rotate around center: lon=0,lat=0 (X-axis)
        x,y,z = sph2cart(lon,lat)
        x,y,z = XAxisRotate(x,y,z,theta)

        #Are we rotating back to horizontal? ("unfolding")
        #if so, rotate around north pole (Z-axis)
        if dip != 0:  x,y,z = ZAxisRotate(x,y,z,omega)

        #Back to lat,long and add to coordinate list
        rotLon,rotLat = cart2sph(x,y,z)
        X.append(rotLon)
        Y.append(rotLat)

    return X,Y

#--Conversions btw cartesian & spherical-------------
def sph2cart(lon, lat): 
    """Converts a long, lat pair in degrees to cartesian 
    coordinates <x,y,z> assuming a radius of 1."""
    lat, lon = radians(lat), radians(lon)
    x = cos(lat)*cos(lon)
    y = cos(lat)*sin(lon)
    z = sin(lat)
    return x,y,z
def cart2sph(x,y,z): 
    """Converts a <x,y,z> triplet to spherical coordinates.
    Returns long, lat in degrees."""
    r = sqrt(x*x + y*y + z*z)
    lat = asin(z/r)
    lon = atan2(y,x)
    return degrees(lon),degrees(lat)


def FindQuadrant(strike):
    """Returns the quadrant of a given strike value 
    between 0 and 360 degrees.  Quadrant value 
    returned is one of: 'I', 'II', 'III', or 'IV'. """ 

    #1st Quadrant
    if (strike >= 0) & (strike <=90):
        quad = 'I'
    #2nd Quadrant
    elif (strike > 90) & (strike <= 180):
        quad = 'II'
    #3rd Quadrant
    elif (strike > 180) & (strike <= 270):
        quad = 'III'
    #4th Quadrant
    elif (strike > 270) & (strike <=360):
        quad = 'IV'
    else:
        raise InputError('%.2f is not a valid strike. Strike must be between 0 and 360!' % strike)
    return quad


def FindEnd(strike):
    """Returns the longitude and latitude of the 
    starting point for the plane defined by the 
    given strike measurement (as a float, 0-360)"""

    #Which point on the "globe" are we working with?
    quad = FindQuadrant(strike)
    if quad in ['I', 'II']:
        pole_lon = 90
        pole_lat = 90-strike
    elif quad in ['III','IV']:
        pole_lon = -90
        pole_lat = strike-270

    return (pole_lon,pole_lat)

def GetDir(dipstring):
    """Given a string with a number and a direction,
    (e.g. 45N or 52NW) return the string indicating
    the direction and the string indicating the angle."""

    #Note, this is used for both dips and strikes.
    #Don't put a sanity check 0<dip<90!

    #--Strip trailing chars
    dipstring = dipstring.strip()

    #--Get Letters from End Indicating Dip Direction
    if dipstring[-2:].isalpha(): #One Letter
        #e.g. 340/45NE
        dipdir = dipstring[-2:]
        dip = dipstring[0:-2].strip()

    elif dipstring[-1:].isalpha(): #Two Letters
        #e.g. 340/45E or 340/45N
        dipdir = dipstring[-1:]
        dip = dipstring[0:-1].strip()

    else: #Just dip, no direction
        #Assume RHR (e.g. 340/45)
        dipdir = ''
        dip = dipstring.strip()
    
    try: dip = float(dip)
    except: raise InputError("%s is not a valid dip.  Can't convert %s to a number!" % (dipstring, dip) )
   
    #--Captialize dipdir for easier checking
    dipdir = dipdir.upper()

    return (dip, dipdir)

def isSameEnd(bearing,plungeDir):
    """Checks to see if the plungeDir matches the bearing
    given.  plungeDir is a string, bearing is a float.
    Returns True if the bearing matches the direction,
    and False if it doesn't."""

    quad = FindQuadrant(bearing)

    #Correct end of line for each quadrant
    Directions = {'I':['N','E','NE'], \
                 'II':['S','E','SE'],\
                'III':['S','W','SW'],\
                 'IV':['N','W','NW']}
    #Is the direction given opposite the bearing?
    if plungeDir:
        same_end = plungeDir in Directions[quad]
    else:
        #Empty string, assume they're the same. This is deliberate
        same_end = True

    return same_end


#--Error Classes----------------------------------------------------------------
class InputError(Exception):
    """Exception raised for errors in input format."""
    def __init__(self, message):   self.message = message
    def __str__(self):             return self.message

#-------------------------------------------------------------------------------
#--Unless imported, execute main function---------------------------------------
#-------------------------------------------------------------------------------
if __name__ == "__main__":
    try:
        main(sys.argv)

    #gracefully exit on ctrl-c
    except KeyboardInterrupt: 
        pass


