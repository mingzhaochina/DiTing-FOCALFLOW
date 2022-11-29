#!/usr/bin/env python
import hashpy_utilities as hu
import os,re,glob
import argparse as ap
from obspy import UTCDateTime,read,Stream
import matplotlib.pyplot as plt
import numpy as np
import logging

def dist_calc(lon1, lat1, lon2, lat2):
    lon1 = np.radians(lon1)
    lon2 = np.radians(lon2)
    lat1 = np.radians(lat1)
    lat2 = np.radians(lat2)

    # Haversine formula
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = np.sin(dlat / 2.0) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin(dlon / 2.0) ** 2
    r = 6371.0  # km
    distance = r * 2.0 * np.arcsin(np.sqrt(a))

    return distance

def find_latlon_from_stname(stname,hypfile):
    """Get latlon from hypfile"""
    with open(hypfile, "r") as hf:
        for line in hf:
            if stname in line:
                f = line.split()
                stlo= f[4]
                stla= f[3]
    return stlo,stla

def get_spratio(traces, picktime, dist):
    logging.debug("This is get_spratio")
    """
    Consumes seismic traces, station, pick and distance information
    Automatically dertmines the amplitude of P and S amplitudes
    Asks the user for confirmation or input
    """
    cuttime = 3.0
    starttime = picktime - 5
    endtime = starttime + 50 + dist / 8.0

    tr = (
        traces.slice(starttime=starttime - 10, endtime=endtime + 10)
        .detrend("demean")
        .integrate()
        .detrend("linear")
        .slice(starttime=starttime, endtime=endtime)
    )
    df = tr[0].stats.sampling_rate
    pickindex = picktime - starttime - cuttime
    try:
        sumtrace = np.sqrt(tr[0].data ** 2 + tr[1].data ** 2 + tr[2].data ** 2)
    except (IndexError, ValueError):
        logging.warning(
            "Could not Sum trace. Maybe one trace is missing or " + "corrupt."
        )
        return 0
    sumtrace = sumtrace / max(sumtrace) * 100


    trfil = tr.filter("highpass", freq=1).slice(
        starttime=starttime + cuttime, endtime=endtime
    )
    try:
        sumtracefil = np.sqrt(
            trfil[0].data ** 2 + trfil[1].data ** 2 + trfil[2].data ** 2
        )
    except (IndexError, ValueError):
        logging.warning(
            "Could not Sum trace. Maybe one trace is missing or " + "corrupt."
        )
        return 0
    sumtracefil = sumtracefil / max(sumtracefil) * 100

    p, s = AmplitudePicker(sumtracefil, pickindex, df, dist)

    sp = s / p
    return sp

def AmplitudePicker(trace, pickindex, df, dist):
    logging.debug("This is AmplitudePicker")
    dt = (dist / 3.5) - (dist / 6)  # in seconds
    p_wind = 2  # seconds, P window
    s_wind = 10  # seconds, s window
    p = max(trace[int((pickindex - 1) * df) : int((pickindex + p_wind) * df)])
    s = max(trace[int((pickindex + dt) * df) : int((pickindex + dt + s_wind) * df)])
    return p, s



if __name__ == "__main__":
    parser = ap.ArgumentParser(
        prog='CEA_Motion_test.py',
        description='Automatic picking of seismic waves using'
                    'Generalized Phase Detection')
    parser.add_argument(
        '-In',
        type=str,
        default='./pick_motion.txt',
        help='real phase association file')
    parser.add_argument(
        '-O',
        type=str,
        default='./result',
        help='Output dir')
    parser.add_argument(
        '-P',
        type=str,
        default=False,
        help='Plot figures')
    args = parser.parse_args()

#generate hash format file
    if not os.path.exists(args.O):
        os.mkdir(args.O)

    with open(args.In, "r") as csvfile2:
        f_csv = csvfile2.readlines()
        for i in range(len(f_csv)):
            line = f_csv[i].strip('\n')
            lines = re.split('\s+', str(line))
            if lines[0] == "#":
                print ("processing: "+str(i)+": "+line)
                if i>1 and args.P:
                    plt.title(event_id, fontsize=12, color='r')
                    values=np.arange(0,idx,1)
                    plt.yticks(values*5,['%d' % val for val in values])
                    plt.xlabel("sample points", fontsize=10)
                    plt.ylabel("trace numbers", fontsize=10)
                    fig_save_name = args.O + "/" + event_id + "/" + event_id + ".motion.png"
                    plt.savefig(fig_save_name, dpi=300)
                    cmd2 = "cp " + fig_save_name + " ."
                    print (cmd2)
                    os.system(cmd2)
                    plt.close()
                if args.P:
                    plt.figure(figsize=(8, 12))
                    idx=0
                event_id= lines[1].split("/")[-1]
                if not os.path.exists(args.O+"/"+event_id):
                    os.mkdir(args.O+"/"+event_id)

                t_year = event_id[0:4]
                t_month = event_id[4:6]
                t_day = event_id[6:8]
                t_hour = event_id[8:10]
                t_min = event_id[10:12]
                t_sec = event_id[12:14]
                evla=lines[2]
                evlo = lines[3]
                evdp = lines[4]
                filepath=lines[1]
                header = (
                        "{:04d} {:02d} {:02d} {:02d} {:02d} {:05.2f} "
                        + "{:010.6f} {:010.6f} {:09.6f} {:05.2f} "
                        + "{:05.2f} {:}\n"
                ).format(
                    int(t_year),
                    int(t_month),
                    int(t_day),
                    int(t_hour),
                    int(t_min),
                    float(t_sec),
                    float(evla),
                    float(evlo),
                          float(evdp),
                    0.0,   # horizontal uncertainty,not applied so far
                    0.0,    # depth uncertainty,not applied so far
                    event_id
                )
                with open(args.O + "/" + event_id + "/" + event_id + "." + "pol.hash", "w") as hf:
                    hf.writelines(header)
                with open(args.O + "/" + event_id + "/" + event_id + ".inp", "w") as cf:
                    cf.write(args.O + "/" + event_id + "/" + event_id + ".pol.hash\n")
                    cf.write("./stations.nll\n")
                    cf.write(args.O + "/" + event_id + "/" + event_id + ".fps\n")
                    cf.write(args.O + "/" + event_id + "/" + event_id + ".all.fps\n")
                    cf.write(args.O + "/" + event_id + "/" + event_id + ".rays\n")
                    cf.write("1\n")   #dang:  Angle increment for grid search
                    cf.write("50\n")   #nmc: 50     Number of perutbations of take-off angles for different source,number of possible azimuth/takeoff angle pairs given for each station 
                    cf.write("20\n")   #maxout: 20  Maximum number focal mechanisms that match misfit critria to
                    cf.write("100\n")   #  delmax(maximum distance)  
                    cf.write("4\n")     #badpol,  number of polarities assumed bad
                    cf.write("0.2\n")   #qbadfac: 0.2,log10 of uncertainty factor for s/p ratios.
                    cf.write("45\n")   #cangle: 45,  mechanisms are “close” if less than this angle apart (degrees)
                    cf.write("0.2\n")   #prob_max: 0.2,  probability threshold for multiples (e.g., 0.1)
                    cf.write("1\n")
                    cf.write("ca.forhash\n")
            else:
                net =lines[0]
                stname = lines[1]
                polar = lines[4]
                p_pick = lines[2]
                p_prob = lines[3]
                sharpness = lines[5]
                utc_p_time = UTCDateTime(p_pick)
                if polar == 'D':
                    pol="-"
                elif polar == 'U':
                    pol = "+"
                else:
                    pol = "x"
                sharpness=lines[5]
                if sharpness == 'I':
                    qp = 0
                else:
                    qp = 1

                # S/P ratio:traces
                trace = read(filepath + "/" + net + "." + stname + "*", format='SAC').copy()
                try:
                    dist = trace[0].stats.sac.dist
                except:
                    stfile="./stations.nll"
                    stlo,stla= find_latlon_from_stname(stname,stfile)
                    dist =dist_calc(float(evlo), float(evla), float(stlo), float(stla))
                try:
                    sp = get_spratio(trace, utc_p_time, dist)
                except:
                    sp = 0
                #sp=0
                outline = "{:5s} {:1s} {:1d} {:08.3f}\n".format(
                    stname, pol, qp, sp
                )
                with open(args.O + "/" + event_id + "/" + event_id + "." + "pol.hash", "a") as hf:
                    hf.writelines(outline)

                if args.P:
                    
                    print (filepath + "/" + net + "." + stname)
                    try:
                        ztrace = read(filepath + "/" + net + "." + stname + "*Z.sac", format='SAC')
                    except:
                        continue
                    ztrace = ztrace.slice(utc_p_time - 0.64, utc_p_time + 0.64)
                    
                    if len(ztrace) ==0:
                        continue
                    input_data = ztrace[0].data[:128]
                    sigPower = sum(np.square(input_data[64:79])) / 15
                    noisePower = sum(np.square(input_data[14:64])) / 50
                    
                    snr = sigPower / noisePower
                    
                    if np.isnan(snr):
                        continue
                    input_data -= np.mean(input_data)
                    max_norm = np.std(input_data)
                    if max_norm == 0:
                        max_norm = 1
                    input_data /= max_norm

                    plt.plot(input_data[:] + idx * 5, linewidth=0.5, color='k')
                    plt.text(0, idx * 5, '{}'.format(polar), fontsize=12, color='b')
                    plt.text(30, idx * 5, 'SNR {:.3f}'.format(snr), fontsize=12, color='b')
                    plt.text(15, idx * 5, '{}'.format(sharpness), fontsize=12, color='b')
                    plt.text(134, idx * 5,
                            net + "." + stname, fontsize=10)
                    
                    if i == len(f_csv)-1:    
                        values=np.arange(0,idx,1)
                        plt.title(event_id, fontsize=12, color='r')
                        plt.yticks(values*5,['%d' % val for val in values])      
                        plt.xlabel("sample points", fontsize=10)
                        plt.ylabel("trace numbers", fontsize=10)
                        fig_save_name = args.O + "/" + event_id + "/" + event_id + ".motion.png"
                        plt.savefig(fig_save_name, dpi=300)
                        cmd2 = "cp " + fig_save_name +  " ."
                        os.system(cmd2)
                        plt.close()

                    #ax.set_title(event_id, fontsize=12, color='r')

                    idx = idx + 1


    # generate stationfile,it is NonLinLoc station file format

    # Run HASH #
    ############
    filelist = glob.glob(args.O + "/*")
    for file in filelist:
        print(file)
        ID=file.split("/")[-1]
        ctrlfile=file+"/"+ID+".inp"
        hashfile=file+"/"+ID+".pol.hash"
        print (ctrlfile,hashfile,ID)
        try:
            strike, dip, rake, misfit = hu.RunHASH(ctrlfile)
        except:
            continue
        print('Strike: Dip: Rake: RMS: (deg)')
        print('{: 6.0f}  {: 3.0f}  {: 4.0f}  {: 3.0f}'.format(
             strike, dip, rake, misfit))
        hu.PlotMechanism(file+"/", hashfile, ID)
        if os.path.exists(file+"/"+ID+".ps"):
            #cmd="ps2epsi "+file+"/"+ID+".ps"
            #os.system(cmd)
            cmd1 = "convert " + file + "/" + ID + ".ps "+ID + ".focal.png"
            #os.system(cmd1)
