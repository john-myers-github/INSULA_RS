"""
A set of functions that you
"""
import numpy as np
import statsmodels.api as sm
import pdb
import numexpr
from scipy.signal import argrelmax
from scipy.stats import ttest_ind
from xarray import concat
from ptsa.data.timeseries import TimeSeries
from ptsa.data.filters import MorletWaveletFilter
from tqdm import tqdm
import pycircstat
from numpy import exp
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from math import cos,sin



def par_find_peaks_by_chan(p_spect_array, frequencies, std_thresh=1.):
    """
    Parameters
    ----------
    p_spect_array: numpy.ndarray
        An array with dimensions frequencies x channels
    frequencies: numpy.ndarray
        An array of the frequencies used
    std_thresh: float
        Threshold in number of standard deviations above the corrected power spectra to be counted as a peak

    Returns
    -------
    peaks_all_chans: numpy.ndarray with type bool
        An array of booleans the same shape as p_spect_array, specifying if there is a peak at a given frequency
        and electrode
    """

    peaks_all_chans = np.zeros(p_spect_array.shape).astype(bool)
    for i, chan_data in enumerate(p_spect_array.T):
        x = sm.tools.tools.add_constant(np.log10(frequencies))
        model_res = sm.RLM(chan_data, x).fit()
        peak_inds = argrelmax(model_res.resid)
        peaks = np.zeros(x.shape[0], dtype=bool)
        peaks[peak_inds] = True
        above_thresh = model_res.resid > (np.std(model_res.resid) * std_thresh)
        peaks_all_chans[:,i] = peaks & above_thresh
    return peaks_all_chans


def par_robust_reg(info):
    """
    Parallelizable robust regression function

    info: two element list. first element, power spectra: # freqs x # elecs. Second element: log transformed freqs

    returns intercepts, slopes, resids
    """

    p_spects = info[0]
    x = sm.tools.tools.add_constant(info[1])

    # holds slope of fit line
    slopes = np.empty((p_spects.shape[1]))
    slopes[:] = np.nan

    # holds residuals
    resids = np.empty((p_spects.shape[0], p_spects.shape[1]))
    resids[:] = np.nan

    # holds intercepts
    intercepts = np.empty((p_spects.shape[1]))
    intercepts[:] = np.nan

    # holds mean height of fit line
    bband_power = np.empty((p_spects.shape[1]))
    bband_power[:] = np.nan

    # loop over every electrode
    for i, y in enumerate(p_spects.T):
        model_res = sm.RLM(y, x).fit()
        intercepts[i] = model_res.params[0]
        slopes[i] = model_res.params[1]
        bband_power[i] = model_res.fittedvalues.mean()
        resids[:, i] = model_res.resid

    return intercepts, slopes, resids, bband_power


def par_robust_reg_no_low_freqs(info):
    """
    Parallelizable robust regression function

    info: two element list. first element, power spectra: # freqs x # elecs. Second element: log transformed freqs

    returns intercepts, slopes, resids
    """

    p_spects = info[0]
    x = sm.tools.tools.add_constant(info[1])
    freq_inds = info[2]

    # holds slope of fit line
    slopes = np.empty((p_spects.shape[1]))
    slopes[:] = np.nan

    # holds residuals
    resids = np.empty((p_spects.shape[0], p_spects.shape[1]))
    resids[:] = np.nan

    # holds intercepts
    intercepts = np.empty((p_spects.shape[1]))
    intercepts[:] = np.nan

    # holds mean height of fit line
    bband_power = np.empty((p_spects.shape[1]))
    bband_power[:] = np.nan

    # loop over every electrode
    for i, y in enumerate(p_spects.T):
        model_res = sm.RLM(y[freq_inds], x[freq_inds]).fit()
        intercepts[i] = model_res.params[0]
        slopes[i] = model_res.params[1]
        bband_power[i] = model_res.fittedvalues.mean()
        resids[:, i] = y - ((x[:, 1]*model_res.params[1]) + model_res.params[0])

    return intercepts, slopes, resids, bband_power


def my_local_max(arr):
    """
    Returns indices of local maxima in a 1D array. Unlike scipy.signal.argrelmax, this does not ignore consecutive
    values that are peaks. It finds the last repetition.

    """
    b1 = arr[:-1] <= arr[1:]
    b2 = arr[:-1] > arr[1:]
    k = np.where(b1[:-1] & b2[1:])[0] + 1
    if arr[0] > arr[1]:
        k = np.append(k, 0)
    if arr[-1] > arr[-2]:
        k = np.append(k, len(arr) - 1)
    return k


def par_find_peaks(info):
    """
    Parallelizable peak picking function, uses robust reg but returns

    """

    # def moving_average(a, n=3):
    #     ret = np.cumsum(a, dtype=float)
    #     ret[n:] = ret[n:] - ret[:-n]
    #     return ret[n - 1:] / n

    p_spect = info[0]
    x = sm.tools.tools.add_constant(info[1])
    model_res = sm.RLM(p_spect, x).fit()
    peak_inds = my_local_max(model_res.resid)
    peaks = np.zeros(x.shape[0], dtype=bool)
    peaks[peak_inds] = True
    above_thresh = model_res.resid > np.std(model_res.resid)
    peaks = peaks & above_thresh
    return peaks

def circ_lin_regress(phases, coords, theta_r, params):
    """
    Performs 2D circular linear regression.

    This is ported from Honghui's matlab code.

    :param phases:
    :param coords:
    :return:
    """

    n = phases.shape[1]
    pos_x = np.expand_dims(coords[:, 0], 1)
    pos_y = np.expand_dims(coords[:, 1], 1)

    # compute predicted phases for angle and phase offset
    x = np.expand_dims(phases, 2) - params[:, 0] * pos_x - params[:, 1] * pos_y

    # Compute resultant vector length. This is faster than calling pycircstat.resultant_vector_length
    # now = time.time()
    x1 = numexpr.evaluate('sum(cos(x) / n, axis=1)')
    x1 = numexpr.evaluate('x1 ** 2')
    x2 = numexpr.evaluate('sum(sin(x) / n, axis=1)')
    x2 = numexpr.evaluate('x2 ** 2')
    Rs = numexpr.evaluate('-sqrt(x1 + x2)')


    # for each time and event, find the parameters with the smallest -R
    min_vals = theta_r[np.argmin(Rs, axis=1)]

    sl = min_vals[:, 1] * np.array([np.cos(min_vals[:, 0]), np.sin((min_vals[:, 0]))])
    offs = np.arctan2(np.sum(np.sin(phases.T - sl[0, :] * pos_x - sl[1, :] * pos_y), axis=0),
                      np.sum(np.cos(phases.T - sl[0, :] * pos_x - sl[1, :] * pos_y), axis=0))
    pos_circ = np.mod(sl[0, :] * pos_x + sl[1, :] * pos_y + offs, 2 * np.pi)

    # compute circular correlation coefficient between actual phases and predicited phases
    circ_corr_coef = pycircstat.corrcc(phases.T, pos_circ, axis=0)

    # compute adjusted r square
    # pdb.set_trace()
    r2_adj = circ_corr_coef ** 2
    # r2_adj = 1 - ((1 - circ_corr_coef ** 2) * (n - 1)) / (n - 4)

    wave_ang = min_vals[:, 0]
    wave_freq = min_vals[:, 1]
    return wave_ang, wave_freq, r2_adj,offs


def channel_peak(data):
    params,rs=findpeak(data)
    peaks = np.zeros(len(data)).astype(bool)
    #params=[ int(round(i)) for i in params]
    if len(params)==5:
        if 5<params[3]<195 and params[2]>.3 and params[4]<50:
            peaks[int(round(params[3]))]=True
    elif len(params)==8:
        if 5<params[3]<195 and params[2]>.3 and params[4]<50:
            peaks[int(round(params[3]))]=True
        if 5<params[6]<195 and params[5]>.3 and params[7]<50:
            peaks[int(round(params[6]))]=True
    elif len(params)==11:
        if 5<params[3]<195 and params[2]>.3 and params[4]<50:
            peaks[int(round(params[3]))]=True
        if 5<params[6]<195 and params[5]>.3 and params[7]<50:
            peaks[int(round(params[6]))]=True
        if 5<params[9]<195 and params[8]>.3 and params[10]<50:
            peaks[int(round(params[9]))]=True
    return peaks
def findpeak(data):
    line,peaks,heights=robustfit(data)
    #when no peaks found from robustfit, return no peaks, extremely rare ad extremely annoying!!!
    if len(peaks)==0 or line[1]>0 or line[0]<0:
        return line,0
    x=np.array(list(range(len(data))))
    y=data
    # Just fit 1 peak
    initial=[line[1],line[0],1,peaks[0],5]
    upper=[0,np.inf,np.inf,199,np.inf]
    lower=[-np.inf,0,0,0,0]
    best_vals,cov = curve_fit(oneOverF1, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF1(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,5,200)
    res=best_vals
    RS=Rs
    # if 1 peaks are good enough, return
    if RS>.999:
        return res,RS

    # fit 2 peaks
    upper=[0,np.inf,np.inf,199,np.inf,np.inf,199,np.inf]
    lower=[-np.inf,0,0,0,0,0,0,0]
    if len(peaks)>=2:
        initial=[line[1],line[0],1,peaks[0],5,1,peaks[1],5]
        best_vals,cov = curve_fit(oneOverF2, x, y,p0=initial,bounds=(lower,upper))
        Rs=np.corrcoef(y,oneOverF2(x,*best_vals))[0,1]**2;
        Rs=Radjust(Rs,8,200)
        if RS<Rs:
            res=best_vals
            RS=Rs
        
    # fit 2 peaks with interpolation
    # first interpolation point
    initial=[line[1],line[0],1,peaks[0],5,1,peaks[0]/2,5]
    best_vals,cov = curve_fit(oneOverF2, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF2(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,8,200)
    if RS<Rs:
        res=best_vals
        RS=Rs

    # second interpolation point
    initial=[line[1],line[0],1,peaks[0],5,1,(peaks[0]+199)/2,5]
    best_vals,cov = curve_fit(oneOverF2, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF2(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,8,200)
    if RS<Rs:
        res=best_vals
        RS=Rs
    
    # if 2 peaks are good enough, return
    if RS>.999 or len(peaks)<2:
        return res,RS

    # fit 3 peaks
    
    upper=[0,np.inf,np.inf,199,199,np.inf,199,199,np.inf,199,199]
    lower=[-np.inf,0,0,0,0,0,0,0,0,0,0]
    if len(peaks)>2:
        initial=[line[1],line[0],1,peaks[0],5,1,peaks[1],5,1,peaks[2],5]
        best_vals,cov = curve_fit(oneOverF3, x, y,p0=initial,bounds=(lower,upper))
        Rs=np.corrcoef(y,oneOverF3(x,*best_vals))[0,1]**2;
        Rs=Radjust(Rs,11,200)
        if RS<Rs:
            res=best_vals
            RS=Rs

    # fit 2 peaks plus different interpolation points cause the 3rd peak is not that reliable.
    Peak=sorted(peaks[:2])
    # first interpolation point
    initial=[line[1],line[0],1,Peak[0]/2,5,1,Peak[0],5,1,Peak[1],5]
    best_vals,cov = curve_fit(oneOverF3, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF3(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,11,len(data))
    if RS<Rs:
        res=best_vals
        RS=Rs
    
    # second interpolation point
    initial=[line[1],line[0],1,Peak[0],5,1,sum(Peak)/2,5,1,Peak[1],5]
    best_vals,cov = curve_fit(oneOverF3, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF3(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,11,200)
    if RS<Rs:
        res=best_vals
        RS=Rs

    # third interpolation point
    initial=[line[1],line[0],1,Peak[0],5,1,Peak[1],5,1,Peak[1]/2+99.5,5]
    best_vals,cov = curve_fit(oneOverF3, x, y,p0=initial,bounds=(lower,upper))
    Rs=np.corrcoef(y,oneOverF3(x,*best_vals))[0,1]**2;
    Rs=Radjust(Rs,11,200)
    if RS<Rs:
        res=best_vals
        RS=Rs
    return res,RS
def oneOverF1(x,slope,intercept, amp1, cen1, wid1):
     return x*slope+intercept+amp1 * exp(-.5*(x-cen1)**2 / wid1**2)

def oneOverF2(x,slope,intercept, amp1, cen1, wid1,amp2, cen2, wid2):
     return x*slope+intercept+amp1 * exp(-.5*(x-cen1)**2 / wid1**2)+amp2 * exp(-.5*(x-cen2)**2 / wid2**2)
    
def oneOverF3(x,slope,intercept, amp1, cen1, wid1,amp2, cen2, wid2,amp3, cen3, wid3):
     return x*slope+intercept+amp1 * exp(-.5*(x-cen1)**2 / wid1**2)+amp2 * exp(-.5*(x-cen2)**2 / wid2**2)+amp3 * exp(-.5*(x-cen3)**2 / wid3**2)
    
def Radjust(Rs,k,n):
    return 1-(1-Rs)*(n-1)/(n-k-1)

def robustfit(data):
    x = sm.tools.tools.add_constant(np.array(range(200)))
    model_res = sm.RLM(data, x).fit()
    peak_inds = np.squeeze(argrelmax(model_res.resid))
    if peak_inds.size==1:
        peak_inds=np.expand_dims(peak_inds,0)
    peak_inds=sorted(peak_inds,key=lambda i:model_res.resid[i])[::-1]
    return model_res.params,peak_inds[:3],model_res.resid[peak_inds[:3]]

def par_find_peaks_by_chan2(p_spect_array, frequencies, std_thresh=1.):
    """
    Parameters
    ----------
    p_spect_array: numpy.ndarray
        An array with dimensions frequencies x channels
    frequencies: numpy.ndarray
        An array of the frequencies used
    std_thresh: float
        Threshold in number of standard deviations above the corrected power spectra to be counted as a peak

    Returns
    -------
    peaks_all_chans: numpy.ndarray with type bool
        An array of booleans the same shape as p_spect_array, specifying if there is a peak at a given frequency
        and electrode
    """

    peaks_all_chans = np.zeros(p_spect_array.shape).astype(bool)
    for i, chan_data in enumerate(p_spect_array.T):
        peaks_all_chans[:,i] = channel_peak(chan_data)
    return peaks_all_chans


def rbar(x):
    n=len(x)
    x1 = numexpr.evaluate('sum(cos(x) / n, axis=0)')
    x1 = numexpr.evaluate('x1 ** 2')
    x2 = numexpr.evaluate('sum(sin(x) / n, axis=0)')
    x2 = numexpr.evaluate('x2 ** 2')
    Rs = numexpr.evaluate('sqrt(x1 + x2)')
    return Rs
def circ_hist2(data,n=20,x=0,y=0,radius=1,color='k'):
    bins=np.linspace(-np.pi,np.pi,n+1)
    counts,bins=np.histogram(data,bins)
    counts=radius*counts/np.max(counts)
    for i in range(len(counts)):
        plt.plot([x,x+counts[i]*cos(bins[i])],[y,y+counts[i]*sin(bins[i])],color=color,linewidth=.3)
        plt.plot([x,x+counts[i]*cos(bins[i+1])],[y,y+counts[i]*sin(bins[i+1])],color=color,linewidth=.3)
        arch=np.linspace(bins[i],bins[i+1],100)
        plt.plot(x+counts[i]*np.cos(arch),y+counts[i]*np.sin(arch),color=color,linewidth=.5)
    plt.axis('equal')


def circ_hist(data,n=20,x=0,y=0,radius=1,color='k',alpha=.3):
    # circular histogram, data need to be within -pi to pi
    bins=np.linspace(-np.pi,np.pi,n+1)
    counts,bins=np.histogram(data,bins)
    counts=radius*counts/np.max(counts)
    for i in range(len(counts)):
        X=[x,x+counts[i]*cos(bins[i]),x+counts[i]*cos(bins[i+1]),x]
        Y=[y,y+counts[i]*sin(bins[i]),y+counts[i]*sin(bins[i+1]),y]
        plt.fill(X,Y,color=color,alpha=alpha,linewidth=0)
    plt.axis('equal')
def adjust(rs,n,k):
    return  1 - ((1 - rs) * (n - 1)) / (n - 1-k)
