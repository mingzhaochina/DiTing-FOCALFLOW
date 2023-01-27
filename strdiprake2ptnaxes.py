#! /usr/bin/python3
"""
Computes Dip and Azimuth of P-, T-, and N-axis of a focal mechanism given as strike, dip, rake

usage: echo "0 0 0" | ./strdiprake2ptnaxes.py
returns: 45/180 45/000 (DipP/AzimuthP DipT/AzimuthT DipN/AzimuthN)
"""

from sys import stdin, stdout, path
from numpy import arctan2, sqrt, pi, sin, cos, matrix, cross, argmin, argmax
from numpy.linalg import eig
from argparse import ArgumentParser

parser = ArgumentParser(description='Reads a list of strike, dip and rake from stdin and returns orientation of the P-, T-, and N-axes')
parser.add_argument('-u', '--upper', help='Return dip/azimuth of P-, T-, and N-axis on the upper hemisphere of a polar plot', action='store_true')
parser.add_argument('-l', '--lower', help='Return dip/azimuth of P-, T-, and N-axis on the lower hemisphere of a polar plot', action='store_true')
parser.add_argument('-s', '--sin2', help='Return the squared sinuses of the dip of the P-, T-, and N-axes', action='store_true')
parser.add_argument('-p', '--profile', help='Return dip/length of P-, T-, and N-axis projeted on an east-west profile, in the northern sector', action='store_true')
parser.add_argument('-o', '--rotate', help='Rotate moment tensor by this many degree around N-axis', type=float, default=0)
parser.add_argument('-r', '--rake', help='Return the rake as a last argument', action='store_true')
parser.add_argument('-a', '--parse', help='Parse last word unchanged as last argument, overrides --rake', action='store_true')

def rotate2upperhemisphere(azi, dip):
    """ Roates azimuth and dip from lower to upper hemisphere"""
    azi = azi + 180
    if dip < 0:
        dip = -dip
        azi = azi + 180
    return azi, dip

def cleanazi(azi):
    while azi < 0:
        azi = azi + 360
    while azi > 360:
        azi = azi - 360
    return azi

def rotatedptn(s, d, r, a): 
    """ 
    Return the North-, East-, and Down- component
    of the P, T, and N-axis of a double couple source
    defined strike (s), dip (d) and rake (r) in a coordinate system
    rotated by a degree around the N axis
    """
    s *= pi/180
    d *= pi/180
    r *= pi/180
    a *= 2*pi/180

    # Moment tensor (see Aki & Richards)
    mxx = -(sin(d)*cos(r)*sin(2*s) + sin(2*d)*sin(r)*sin(s)**2) #north
    mxy =   sin(d)*cos(r)*cos(2*s) + 0.5*sin(2*d)*sin(r)*sin(2*s)
    mxz = -(cos(d)*cos(r)*cos(s)   + cos(2*d)*sin(r)*sin(s))
    myy =   sin(d)*cos(r)*sin(2*s) - sin(2*d)*sin(r)*cos(s)**2 #east
    myz = -(cos(d)*cos(r)*sin(s)   - cos(2*d)*sin(r)*cos(s))
    mzz =   sin(2*d)*sin(r) #down

    M = matrix([[mxx, mxy, mxz], [mxy, myy, myz], [mxz, myz, mzz]])
    R = matrix([[1, 0, 0], [0, cos(a), -sin(a)], [0, sin(a), cos(a)]])
    Mr = R*M #rotated moment tensor
    w, v = eig(Mr)

    P = v[:, argmin(w)].A1
    T = v[:, argmax(w)].A1
    N = cross(P, T)

    Pn, Pe, Pd = P[0], P[1], P[2]
    Tn, Te, Td = T[0], T[1], T[2]
    Nn, Ne, Nd = N[0], N[1], N[2]

    return Pn, Pe, Pd, Tn, Te, Td, Nn, Ne, Nd

def ned2dipazi(n, e, d):
    """
    Transforms vector given in north, east, down component to azimuth and dip
    """
    dip = (arctan2(sqrt(n**2 + e**2), -abs(d)) - pi/2) * 180/pi
    azi = 360 - arctan2(e, -n) * 180/pi

    if d > 0:
        azi -= 180
    azi = cleanazi(azi)
    return dip, azi

def ned2EWlendir(n, e, d):
    """
    Returns length and direction of the vector given by North-, East-, and Down- component, projected on the EW-plane
    """
    direction = arctan2(d, e) * 180/pi
    length = sqrt(d**2 + e**2) * 999 # scale factor to fulfill prerequists of output format
    if d > 0:
        direction -= 180 
    while direction < 0:
        direction += 360
    return length, direction

def rotate2map(azi, dip):
    length = abs(cos(dip * np.pi/180))*100
    return azi, length

def sin2(angle):
    """returns sin^2 of angle in degree"""
    return sin(angle*pi/180)**2

infile = stdin
outfile = stdout
args = parser.parse_args()

for line in infile:
    line = line.strip()
    if args.parse:
        line, parsant = line.split()[:-1], line.split()[-1]
        line = ' '.join(line)
    strk, dip, rak = map(float, line.split())
    Pn, Pe, Pd, Tn, Te, Td, Nn, Ne, Nd = rotatedptn(strk, dip, rak, args.rotate)
    dip_p, azi_p = ned2dipazi(Pn, Pe, Pd)
    dip_t, azi_t = ned2dipazi(Tn, Te, Td)
    dip_n, azi_n = ned2dipazi(Nn, Ne, Nd)
    if args.upper:
        azi_p, dip_p = rotate2upperhemisphere(azi_p, dip_p)
        azi_t, dip_t = rotate2upperhemisphere(azi_t, dip_t)
        azi_n, dip_n = rotate2upperhemisphere(azi_n, dip_n)
    azi_p = cleanazi(azi_p)
    azi_t = cleanazi(azi_t)
    azi_n = cleanazi(azi_n)
    if args.profile: 
        azi_p, dip_p = ned2EWlendir(Pn, Pe, Pd)
        azi_t, dip_t = ned2EWlendir(Tn, Te, Td)
        azi_n, dip_n = ned2EWlendir(Nn, Ne, Nd) #azi is length, dip is direction
    if args.rake:
        parsant = rak
    if args.parse or args.rake and not args.sin2:
        data = '{:02.0f}/{:03.0f} {:02.0f}/{:03.0f} {:02.0f}/{:03.0f} {:} \n'.format(dip_p, azi_p, dip_t, azi_t, dip_n, azi_n, parsant)
    elif not args.sin2:
        data = '{:02.0f}/{:03.0f} {:02.0f}/{:03.0f} {:02.0f}/{:03.0f}\n'.format(dip_p, azi_p, dip_t, azi_t, dip_n, azi_n)
    if args.sin2:
        sin2p, sin2t, sin2n = sin2(dip_p), sin2(dip_t), sin2(dip_n)
        data = '{:.4f} {:.4f} {:.4f}\n'.format(sin2p, sin2t, sin2n)
    outfile.write(data)
