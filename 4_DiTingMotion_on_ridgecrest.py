#!/usr/bin/env python
# Automatic picking of seismic waves using U-shaped convolutional neural network
# See https://github.com/mingzhaochina for more info
#
# ZHAO Ming, CHEN Shi, FANG LiHua et al.2019.
#Earthquake phase arrival auto-picking based on U-shaped convolutional neural network,
# Chinese Journal of Geophysics(in Chinese),62(8): 3034-3042,doi: 10.6038/cjg2019M0495
#
# Author: Ming Zhao (2021)
# Contact: mzhao@cea-igp.ac.cn
# Website: http://www.neobji.ac.cn/article/542
import os,re,glob
import tensorflow as tf
import sys
import pandas as pd
os.environ['CUDA_VISIBLE_DEVICES'] = '2'
import numpy as np
import argparse as ap
from obspy import UTCDateTime,read,Stream
from obspy.signal import polarization, util
#sys.path.append('./SmartMotionV1.11/smartMotion')
#sys.path.append('./SmartMotionV2.00/smartMotion')
#from model import SmartMotion
#from model_V2 import SmartMotion_V2
from geographiclib.geodesic import Geodesic
from obspy.imaging.beachball import beachball
from src.PostProcessing import postprocesser
import warnings
warnings.filterwarnings("ignore")
#tf.config.list_physical_devices('GPU')`
#tf.test.is_gpu_available()
from keras import backend as K
K.clear_session()


if __name__ == "__main__":
    parser = ap.ArgumentParser(
        prog='CEA_Motion_test.py',
        description='Automatic picking of seismic waves using'
                    'Generalized Phase Detection')
    parser.add_argument(
        '-In',
        type=str,
        default='./hash_phase_ridgcrest_manual.dat',
        help='real phase association file')
    parser.add_argument(
        '-OO',
        type=str,
        default='./pick_motion.txt',
        help='Output dir')
    parser.add_argument(
        '-M',
        type=str,
        default='./models/DiTingMotionJul.hdf5',
        help='Predict')
    args = parser.parse_args()

motion_model =tf.keras.models.load_model(args.M,compile=False)
#data = pd.read_csv(args.In,sep=',')

counter = 0

fo = open(args.OO,'w')
jinmark = False

with open(args.In, "r") as csvfile2:
    f_csv = csvfile2.readlines()
    for i in range(len(f_csv)):
        # for line in f_csv:
        line = f_csv[i].strip('\n')
        lines = re.split('\s+', str(line))
        #  line_next =f_csv[i+1].strip('\n')
        #  lines_next=re.split('\s+', str(line_next))
        if lines[0] == "#":
            t_year = lines[1]
            t_month = lines[2]
            t_day = lines[3]
            t_hour = lines[4]
            t_min = lines[5]
            t_sec = lines[6]
            print ( t_year, t_month, t_day, t_hour, t_min, t_sec)
            t = UTCDateTime(int(t_year), int(t_month), int(t_day), int(t_hour), int(t_min), float(t_sec))
            evla = lines[7]
            evlo = lines[8]
            evdep = lines[9]
            mag = lines[10]
            nums=lines[-1]
            evdir="./cut/"+t_year+t_month+t_day+t_hour+t_min+str(t_sec)           
            if os.path.exists(evdir):
               jinmark = True
               fo.write("{} {} {} {} {} {} {}\n".format("# ", evdir, evla, evlo, evdep, mag, nums))
            else:
               jinmark = False
        else:
            if not jinmark:
                continue  
            stname = lines[0]
            net = lines[1]
            print  (stname, t)
            m2 = re.match('P', lines[3])
            if m2:
                p_timestamp = t.timestamp + float(lines[2])
                p_pick = UTCDateTime(p_timestamp)
                tp_res =float(lines[4])
                t_st=Stream()
                try:
                    t_st = read(evdir+ "/*." + stname + "*Z.sac")
       
                except:
                    continue
                t_st.detrend('demean')
                try:
                    t_st.detrend(type='linear')
                except:
                    t_st.detrend(type='constant')
                t_st = t_st.taper(0.001)
                t_st=t_st.slice(p_pick - 0.64, p_pick + 0.64)
                motion_input = np.zeros([1, 128, 2])
                      
                motion_input[0, :, 0] = t_st[0].data[0:128]
                        #print(trace.data[0:128])
                if np.max(motion_input[0, :, 0]) == 0:
                     pass
                else:
                     motion_input[0, :, 0] -= np.mean(motion_input[0, :, 0])
                     norm_factor = np.std(motion_input[0, :, 0])

                     if norm_factor == 0:
                         pass
                     else:
                         motion_input[0, :, 0] /= norm_factor
                         diff_data = np.diff(motion_input[0, 64:, 0])
                         diff_sign_data = np.sign(diff_data)
                         motion_input[0, 65:, 1] = diff_sign_data[:]
                #print (motion_input)
                pred_res = motion_model.predict(motion_input)
                pred_fmp = (pred_res['T0D0'] + pred_res['T0D1'] + pred_res['T0D2'] + pred_res['T0D3']) / 4
                pred_cla = (pred_res['T1D0'] + pred_res['T1D1'] + pred_res['T1D2'] + pred_res['T1D3']) / 4
                print(pred_fmp, pred_cla)
                if np.argmax(pred_fmp[0, :]) == 1:
                    polarity = 'D'
                    if np.argmax(pred_cla[0, :]) == 0:
                        sharpness = 'I'
                    elif np.argmax(pred_cla[0, :]) == 1:
                        sharpness = 'E'
                    else:
                        sharpness = 'x'
                    fo.write("{} {} {} {} {} {}\n".format(net, stname,
                                                                    p_pick,  tp_res, polarity, sharpness))
                elif np.argmax(pred_fmp[0, :]) == 0:
                    polarity = 'U'
                    if np.argmax(pred_cla[0, :]) == 0:
                        sharpness = 'I'
                    elif np.argmax(pred_cla[0, :]) == 1:
                        sharpness = 'E'
                    else:
                        sharpness = 'x'
                    fo.write("{} {} {} {} {} {}\n".format(net, stname,
                                                                    p_pick, tp_res, polarity, sharpness))
                else:
                    polarity = 'x'
                    if np.argmax(pred_cla[0, :]) == 0:
                        sharpness = 'I'
                    elif np.argmax(pred_cla[0, :]) == 1:
                        sharpness = 'E'
                    else:
                        sharpness = 'x'
                    fo.write("{} {} {} {} {} {}\n".format(net, stname,
                                                          p_pick, tp_res, polarity, sharpness))
                counter += 1
                # break
                if counter % 100 == 0:
                    print('On {}'.format(counter))

fo.close()
# In[ ]:




