from obspy.signal.trigger import trigger_onset
import numpy as np

def _detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising', kpsh=False, valley=False):

    """
    modified from https://github.com/smousavi05/EQTransformer
    Detect peaks in data based on their amplitude and other features.
    Parameters
    ----------
    x : 1D array_like
        data.
        
    mph : {None, number}, default=None
        detect peaks that are greater than minimum peak height.
        
    mpd : int, default=1
        detect peaks that are at least separated by minimum peak distance (in number of data).
        
    threshold : int, default=0
        detect peaks (valleys) that are greater (smaller) than `threshold in relation to their immediate neighbors.
        
    edge : str, default=rising
        for a flat peak, keep only the rising edge ('rising'), only the falling edge ('falling'), both edges ('both'), or don't detect a flat peak (None).
        
    kpsh : bool, default=False
        keep peaks with same height even if they are closer than `mpd`.
        
    valley : bool, default=False
        if True (1), detect valleys (local minima) instead of peaks.
    Returns
    ---------
    ind : 1D array_like
        indeces of the peaks in `x`.
    Modified from 
   ----------------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb
    
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    return ind

def postprocesser(yh1, yh2, yh3, det_th=0.3, p_th=0.3, p_mpd=300, s_th=0.3, s_mpd=300):

    """ 
    modified from https://github.com/smousavi05/EQTransformer
    Postprocessing to detection and phase picking
    """         
             
    detection = trigger_onset(yh1, det_th, det_th)
    pp_arr = _detect_peaks(yh2, mph=p_th, mpd=p_mpd)
    ss_arr = _detect_peaks(yh3, mph=s_th, mpd=s_mpd)
          
    P_PICKS = {}
    S_PICKS = {}
    EVENTS = {}
    matches = {}

    # P
    if len(pp_arr) > 0:
        for pick in range(len(pp_arr)): 
            pauto = pp_arr[pick]
            if pauto: 
                P_prob = np.round(yh2[int(pauto)], 3) 
                P_PICKS.update({pauto : [P_prob]})                 
    # S         
    if len(ss_arr) > 0:            
        for pick in range(len(ss_arr)):        
            sauto = ss_arr[pick]
            if sauto: 
                S_prob = np.round(yh3[int(sauto)], 3) 
                S_PICKS.update({sauto : [S_prob]})             
            
    if len(detection) > 0:
        # merge close detections
        for ev in range(1,len(detection)):
            if detection[ev][0] - detection[ev-1][1] < 50:
                detection[ev-1][1] = detection[ev][1]
        for ev in range(len(detection)):                                 
            D_prob = np.mean(yh1[detection[ev][0]:detection[ev][1]])
            D_prob = np.round(D_prob, 3)
            EVENTS.update({ detection[ev][0] : [D_prob, detection[ev][1]]})            
    
    # matching the detection and picks
    for ev in EVENTS:
        bg = ev
        ed = EVENTS[ev][1]
        """
        if ed > 0.9*len(yh3):
            continue
        if bg < 0.1*len(yh3):
            continue
        """
        S_error = None
        P_error = None        
        if int(ed-bg) >= 100:               
            candidate_Ss = {}
            min_dev = 9999
            peak_prob = -10.0
            closest_Ss = None
            closest_Ss_prob = None
            for Ss, S_val in S_PICKS.items():
                if Ss > ed - 200 and Ss <= ed + 300:
                    temp_dev = np.abs(Ss-ed)
                    if temp_dev < min_dev:
                        min_dev = temp_dev
                        peak_prob = S_val[0]

                        closest_Ss = Ss
                        closest_Ss_prob = S_val
            if min_dev < 9999:
                candidate_Ss = ({closest_Ss : closest_Ss_prob})
            else:
                candidate_Ss = {None:[None]}

            candidate_Ps = {}
            min_dev = 9999
            peak_prob = -10.0
            closest_Ps = None
            closest_Ps_prob = None
            for Ps, P_val in P_PICKS.items():
                if Ps > bg-200 and Ps < ed:
                    temp_dev = np.abs(Ps-bg)
                    if temp_dev < min_dev:
                        min_dev = temp_dev
                        closest_Ps = Ps
                        closest_Ps_prob = P_val
                        peak_prob = P_val[0]
            if min_dev < 9999:
                candidate_Ps = ({closest_Ps : closest_Ps_prob})
            else:
                candidate_Ps = {None:[None]}
                    
            if list(candidate_Ss)[0] or list(candidate_Ps)[0]:

                matches.update({
                                bg:[ed, 
                                    EVENTS[ev][0], 
                                
                                    list(candidate_Ps)[0],  
                                    candidate_Ps[list(candidate_Ps)[0]][0],   
                                                
                                    list(candidate_Ss)[0],  
                                    candidate_Ss[list(candidate_Ss)[0]][0],
                                                ] })
    
    return matches

