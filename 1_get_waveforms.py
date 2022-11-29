#!/usr/bin/env python
#
# Download waveform data of events
#
import obspy
from obspy import read_events
from obspy.clients.fdsn.mass_downloader import (
    CircularDomain,
    Restrictions,
    MassDownloader,
)
# download waveform data of each event
latitude = 35.71
longitude = -117.50
starttime = obspy.UTCDateTime(2019, 7, 4)
endtime = obspy.UTCDateTime(2019, 7, 5)

domain = CircularDomain(
        latitude=latitude, longitude=longitude, minradius=0.0, maxradius=1.0
)

restrictions = Restrictions(
    starttime=starttime,
    endtime=endtime,
    reject_channels_with_gaps=False,
    #channel="LH?",
    #network="CU,GT,IC,II,IU,G", 
    #station="A*", 
    #location="00", 
    minimum_length = 0.95,
    sanitize=False,
    #minimum_interstation_distance_in_m=100e3,
)

mdl = MassDownloader()
mdl.download(domain, restrictions, mseed_storage="mseed", stationxml_storage="stations")
