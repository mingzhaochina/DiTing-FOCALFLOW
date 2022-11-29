#! /usr/bin/env python3

"""
hashpy2 -- an interactive HASH wrapper for meassuring first motion P wave
polarities, P to S amplitude ratios and invert for focal mechanisms
"""

import logging
import yaml
import argparse
import sys
from glob import glob
from os.path import isfile
from numpy import cos, pi
import hashpy_utilities as hu


logging.basicConfig(level=logging.INFO)

par = argparse.ArgumentParser(prog='hashpy2',
                              description='An interactive HASH wrapper',
                              epilog=('Authors: Wasja Bloch ' +
                                      '(wasja@gfz-potsdam.de) ' +
                                      'and Lukas Lehmann ' +
                                      '(luklehma@uni-potsdam.de)'))

par.add_argument('ID', type=str,
                 nargs='?',
                 help="Origin time of event in format: YYYYMMDDhhmmss",
                 default=None)
par.add_argument('--config', type=str,
                 help="Configuration file [config.yaml]",
                 default="config.yaml")
par.add_argument('--setup', action="store_true",
                 help=('Print configuration tipps. ' +
                       'Produce default_config.yaml.'))

print("This is hashpy2")

args = par.parse_args()
ID = args.ID

if not ID:
    par.print_help()
    sys.exit()

#  Read config file
configfile = args.config
with open(configfile, 'r') as stream:
    params = yaml.safe_load(stream)
logging.basicConfig(level=params['LOGGING'])

#  Set directories
resultdir = params['RESULTS'] + '/'
hashfile = resultdir + ID + '.pol.hash'
ctrlfile = resultdir + ID + '.inp'
wvdir = params['WAVEFORMS'] + '/'
stationfile = params['STATIONS']
hypfiles = glob(params['CATALOG'])

# Run HASH #
############
ans = 'y'
if ans == 'y':
    strike, dip, rake, misfit = hu.RunHASH(ctrlfile)
    print('Strike: Dip: Rake: RMS: (deg)')
    print('{: 6.0f}  {: 3.0f}  {: 4.0f}  {: 3.0f}'.format(
        strike, dip, rake, misfit))
    hu.PlotMechanism(resultdir, hashfile, ID)


# the terminal will be open/activated the whole time
# has to be killed manually with the except clause!!
#    try:
#        cmd = ('ID=$(xdotool getactivewindow); while true; do sleep 0.2; ' +
#               'xdotool windowfocus $ID; done;')
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        main(sID)
#    except Exception as e:
#        print(e)
#        cmd = 'pkill -f xdotool'
#        p   = Popen(cmd, shell=True, stdout=PIPE, stderr=PIPE)
#        exit()
