import numexpr
import os
import warnings
import numpy as np
import pandas as pd
import xarray as xr
import bottleneck as bn
import h5py
import math
import sys

from ptsa.data.filters import ButterworthFilter
from ptsa.data.filters import MorletWaveletFilter
from ptsa.data.filters import ResampleFilter
from ptsa.data.timeseries import TimeSeries
from ptsa.data.readers.index import JsonIndexReader
from ptsa.data.readers import BaseEventReader

from cmlreaders import CMLReader, get_data_index
from scipy.stats.mstats import zscore
from scipy.io import loadmat
from scipy.signal import butter, sosfilt
from tqdm import tqdm
from glob import glob

from scipy.stats import norm   
    
def truncate(number, digits) -> float:
    stepper = 10.0 ** digits
    return math.trunc(stepper * number) / stepper





