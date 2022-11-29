#! /usr/bin/python3
from subprocess import Popen, PIPE
from struct import unpack
from glob import glob
import logging
from warnings import warn, filterwarnings
from copy import deepcopy
import numpy as np
from obspy import read_events, read, UTCDateTime
from obspy.signal.util import util_geo_km
from obspy.signal.trigger import recursive_sta_lta, pk_baer

def print_config_tipps():
    """ Print meaning of configuration parameters to screen """
    print('')
    print('----------------------------------')
    print('Control parameters of hashpy2 are:')
    print('')
    print('CATALOG:')
    print('    Regular expression pattern for NonLinLoc .hyp files holding')
    print('    hypocenter, location error and phase information.')
    print('')
    print('RESULTS:')
    print('    Directory in which the results will be stored. Must be present.')
    print('')
    print('WAVEFORMS:')
    print('    Directory in which the waveform files are stored. Directory')
    print('    structure is:')
    print('    WAVEFORMS/YYYY/YYMMDDhhmmss, where')
    print('    - YYYY is four digit year,')
    print('    - YY, MM, DD, hh, mm, ss are two digit year, month, day,')
    print('      hour, minute, and second.')
    print('    Make sure that all three components are present and that they')
    print('    are present only once.')
    print('')
    print('VELOCITIES:')
    print('    List of velocity model files holding two columns:')
    print('    - depth (km), velocity (km/s)')
    print('    Velocities are interpolated linearly. Use zero thickness layers')
    print('    to model discontinuities.')
    print('')
    print('LOGGING:')
    print('    Choose verbosity of program output.')
    print('')
    print('-------------------------------')
    print('Control parameters of HASH are:')
    print('')
    print('dang:')
    print('    Angle increment for grid search.')
    print('')
    print('nmc:')
    print('    Number of perutbations of take-off angles for different source')
    print('    depths and velocity models.')
    print('')
    print('maxout:')
    print('    Maximum number focal mechanisms that match misfit critria to')
    print('    return.')
    print('')
    print('cangle:')
    print('    Angular distance between different families of focal mechanisms.')
    print('')
    print('prob_max:')
    print('    Fraction of focal mechanisms that need to be within cangle')
    print('    to make up a new famlily of focal mechanisms.')
    print('')
    print('qbadfac:')
    print('    log10 of uncertainty factor for s/p ratios.')
    print('')


def write_default_config(filename):
    """
    Write default configuration parameters to file.
    filename: (str) Name of the file
    """
    with open(filename, 'w') as f:
        f.write('CATALOG: ./hyp/*.hyp\n')
        f.write('RESULTS: ./results\n')
        f.write('WAVEFORMS: ./events\n')
        f.write('STATIONS: stations.nll\n')
        f.write('VELOCITIES:\n')
        f.write('    - vmodel.zv\n')
        f.write('LOGGING: WARNING\n')
        f.write('dang: 1\n')
        f.write('nmc: 50\n')
        f.write('maxout: 100\n')
        f.write('cangle: 45\n')
        f.write('prob_max: 0.2\n')
        f.write('qbadfac: 0.2\n')


# Run HASH
################################

def RunHASH(controlfile):
    """Call an instance of HASH"""
    with open(controlfile, 'r') as cf:
        for n, line in enumerate(cf):
            if n == 2:
                resultfile = line.strip()  
                out, err, ret = runBash("./hash_hashpy1D < " + controlfile)
    if ret != 0:
        msg = 'HASH endend with an error:\n' + str(err)
        logging.warning(msg)
    else:
        logging.info('HASH endded successfully')
    logging.info('HASH output:\n' + out.decode('utf8'))
    ang_rms = np.inf
    qual = 'Z'  # worst quality ever
    print (resultfile)
    with open(resultfile, 'r') as rf:
        for nn, line in enumerate(rf):
            if line.split()[4] <= qual and float(line.split()[3]) <= ang_rms:
                strike, dip, rake, ang_rms = map(int, line.split()[0:4])
                qual = line.split()[4]
    if nn > 0:
       logging.warning('Found multiple solutions. Reporting best one.')

    return strike, dip, rake, ang_rms


def runBash(cmd):
    """Run cmd as a shell comand, return stdout, stderr and return code"""
    p = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
    out = p.stdout.read().strip()
    err = p.stderr.read().strip()
    p.poll()
    ret = p.returncode
    return out, err, ret

def PlotMechanism(resultdir, hashfile, ID):
    allfpsf = resultdir + ID + '.all.fps'
    angf = resultdir + ID + '.rays' 
    fpsf = resultdir + ID + '.fps'
    plotfile = resultdir + ID + '.ps'
    cmd = 'sh plot_mechanism.sh {:} {:} {:} {:} {:}'.format(
            hashfile, allfpsf, fpsf, angf, plotfile)
    out, _, _ = runBash(cmd)  # run plotting script
    logging.info('The plotting script returned:\n' + out.decode('utf-8'))
  #  _, _, _ = runBash("gv " + plotfile)  # run gv

