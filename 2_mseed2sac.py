#!/usr/bin/env python
#
# Convert waveform data from miniseed to SAC
#

import os
import obspy
from obspy import read, read_inventory, read_events
from obspy.io.sac.util import utcdatetime_to_sac_nztimes
from obspy.taup import TauPyModel
from obspy.geodetics import gps2dist_azimuth, kilometers2degrees

def obspy_to_sac_header(stream, inventory):
    for tr in stream:
        # Add stats for SAC format
        tr.stats.sac = dict()

        # Add station and channel information
        metadata = inventory.get_channel_metadata(tr.id)
        tr.stats.sac.stla = metadata["latitude"]
        tr.stats.sac.stlo = metadata["longitude"]
        tr.stats.sac.stel = metadata["elevation"]
        tr.stats.sac.stdp = metadata["local_depth"]
        tr.stats.sac.cmpaz = metadata["azimuth"]
        tr.stats.sac.cmpinc = metadata["dip"] + 90 # different definitions

        # Add event information
        tr.stats.sac.o = 0.0

st = read("./mseed/*.mseed")
inv = read_inventory("stations/*.xml")
st.attach_response(inv)
obspy_to_sac_header(st, inv)
#st.remove_response(inventory=inv)
st.rotate(method="->ZNE", inventory=inv)

os.makedirs("sac", exist_ok=True)
for tr in st:
    #tr.interpolate(sampling_rate=0.5)
    tr.write(f"sac/{tr.id}", format="SAC")
