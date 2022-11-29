#! /usr/bin/python3
from obspy import UTCDateTime,read
from obspy.io.sac import SACTrace
import re,glob,os
import numpy as np
fp=open("./SearchResults_example.txt","r"); fplines=fp.readlines(); fp.close()
hf=open("./hash_phase_ridgcrest_manual.dat", "w")
f=open("./catalog_ridgcrest_manual.dat", "w")

ev_id=0
nn = 0
nums=0
#dirlist= glob.glob("../2019070*")
mindex=0
ddirwa="./sac"
outdir ="./cut"
for index in range(len(fplines)):
    fpline = fplines[index].strip('\n')
    #print (len(fpline))
    if (len(fpline) >= 164):
        nn = nn + 1            
        if fpline[18] == ' ':  # N
            lat_manual = (float(fpline[16:18]) + float(fpline[19:23]) / 6000)
        else:
            lat_manual = float(fpline[16:18]) + float(fpline[19:23]) / 6000 * (-1)

        if fpline[26] == 'E':
            lon_manual = (float(fpline[23:26]) + float(fpline[27:31]) / 6000)
        else:
            lon_manual = (float(fpline[23:26]) + float(fpline[27:31]) / 6000) * (-1)
        mag_num=re.sub(u"([^\u0030-\u0039])", "", fpline[123:130])
        print (mag_num)
        mag = float(mag_num) / 100
        t_year_manual = fpline[0:4]
        t_month_manual = fpline[4:6]
        t_day_manual = fpline[6:8]
        t_hour_manual = fpline[8:10]
        t_min_manual = fpline[10:12]
        t_sec_manual = int(fpline[12:16]) / 100
        print (int(t_year_manual), int(t_month_manual), int(t_day_manual), int(t_hour_manual), int(t_min_manual), float(t_sec_manual))
        t0 = UTCDateTime(int(t_year_manual), int(t_month_manual), int(t_day_manual), int(t_hour_manual),
                             int(t_min_manual), float(t_sec_manual))
        t_manual=str(t0)
        dep_manual = float(fpline[31:36])/100
        RMS = float(fpline[48:52]) / 100
        gap = int(fpline[42:45])
        EZ = float(fpline[89:93]) / 100
        EH = float(fpline[85:89]) / 100
        num_ps = int(fpline[39:42])
        num_s = int(fpline[82:85])
        num_p_polar = int(fpline[93:96])
        evid=int(fpline[136:146])
        print (t_year_manual, t_month_manual, t_day_manual, t_hour_manual, t_min_manual, t_sec_manual,
        lat_manual, lon_manual, dep_manual, mag, EH, EZ, RMS,nn,gap, num_ps,num_s,num_p_polar)
        hf.write(
            '# {:4} {:2} {:2} {:2} {:2} {:5.2f}  {:7.4f} {:9.4f} {:5.2f} {:5.2f} {:5.2f} {:5.2f} {:5.2f} {:3d} {:3d} {:3d} {:3d} {:9d}\n'.format(
                t_year_manual, t_month_manual, t_day_manual, t_hour_manual, t_min_manual, t_sec_manual,
                lat_manual, lon_manual, dep_manual, mag, EH, EZ, RMS,gap, num_ps,num_s,num_p_polar,nn))
        f.write(
            '{:4} {:2} {:2} {:2} {:2} {:5.2f}  {:7.4f} {:9.4f} {:5.2f} {:5.2f} {:5.2f} {:5.2f} {:5.2f} {:3d} {:3d} {:3d} {:3d} {:9d}\n'.format(
                t_year_manual, t_month_manual, t_day_manual, t_hour_manual, t_min_manual, t_sec_manual,
                lat_manual, lon_manual, dep_manual, mag, EH, EZ, RMS,gap, num_ps,num_s,num_p_polar,nn))
        for j in range(index+1,index+1+num_ps):
            pline = fplines[j].strip('\n')
            station = pline[0:5].strip()
            net = pline[5:7]
            if pline[14:15] == 'P':
                p_weight = float(pline[38:41])
                p_residual = abs(int(pline[34:38]) / 100)
                p_pick = UTCDateTime(pline[17:29]) + int(pline[29:34]) / 100
                t_diff=p_pick-t0
                p_polarity = pline[15:16]
                distance = float(pline[74:78]) / 10
                if pline[13:14] == 'e':
                    sharpness_P = "E"
                elif pline[13:14] == 'i':
                    sharpness_P = "I"
                else:
                    sharpness_P = "x"
                hf.write(
                    '{:4} {:2} {:7.4f} {} {:4.2f} {:2} {:2} {:5.2f}\n'.format(
                        station, net, t_diff,pline[14:15], p_residual, p_polarity,sharpness_P,distance))
                try:
                    st_three =read(ddirwa  + '/'+'*.'+station+'.*').\
                                    slice(starttime=t0-5, endtime=p_pick + 40)
                except:
                    continue
                print (st_three)
                st_three.sort(reverse=True)
                stz=st_three[0]
                try:
                    stz= stz.detrend('linear')
                except:
                    stz= stz.detrend('constant')
                sac = SACTrace.from_obspy_trace(stz)

                sac.evdp = dep_manual
                sac.evla = lat_manual

                sac.evlo = lon_manual
                sac.mag = mag
                sac.kevnm =t_year_manual+t_month_manual+t_day_manual+t_hour_manual+t_min_manual+str(t_sec_manual)
                sac.t0 = t0
                sac.t1 =p_pick
                fname = stz.id+'.sac'
                tr = sac.to_obspy_trace()
                evdir=outdir +"/"+t_year_manual+t_month_manual+t_day_manual+t_hour_manual+t_min_manual+str(t_sec_manual)
                if not os.path.exists(evdir):
                    os.makedirs(evdir)          
                tr.write(evdir+"/"+fname)
                try:
                    st_three[1].write(evdir+"/"+st_three[1].id + '.sac')
                except:
                    continue
                try:
                    st_three[2].write(evdir+"/"+st_three[2].id + '.sac')
                except:
                    continue
            if pline[47:48] == 'S':
                s_weight = float(pline[63:66])
                s_residual = abs(int(pline[50:54]) / 100)
                s_pick = UTCDateTime(pline[17:29]) + int(pline[41:46]) / 100
                distance = float(pline[74:78]) / 10
                t_diff_s = s_pick-t0
                if pline[46:47] == 'e':
                    sharpness_S = "E"
                elif pline[46:47] == 'i':
                    sharpness_S = "I"
                else:
                    sharpness_S = "x"
                hf.write(
                    '{:4} {:2} {:7.4f} {} {:4.2f} {:2} {:2} {:5.2f}\n'.format(
                        station, net, t_diff_s, pline[47:48],  s_residual, " ", sharpness_S, distance))

hf.close()
f.close()

