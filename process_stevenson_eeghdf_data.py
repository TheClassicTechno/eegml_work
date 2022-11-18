#%% [markdown]
# ### Example of preprocessing code for stevenson neonatal EEG dataset
# author: Lucas Orts, 2019
# from [source repo](https://gitlab.lee-messer.net/eegml/leorts-eeg-clustering-2019/-/blob/master/stevenson_processing.py)
#
# Lucas's strategy I believe was to download and preprocess the EEG data into
# whitened data and save it as .npy files for easy loading

#%%
import os
import os.path as op
from IPython.display import clear_output

import numpy as np
import pandas as pd
import h5py
import eeghdf

def load_stevenson_eeg(path, selected_eegs):
    """Extracts the stevenson data from the hdf folder and csv files.
    Parameters
    ----------
    path : path to directory that contains the hdf and csv files
    selected_eegs : the range of eegs to be extracted
    Returns
    -------
    hdf_files : a list of the h5py files
    acc_annot : a pandas dataframe containing the accumulative
                agreement of the experts
    """
    
    hdf_files = []
    for i in selected_eegs:
        hdf_files.append(h5py.File('%s/hdf/eeg%s.eeg.h5' % (path, i), 'r'))
    
    annots = [('annotations_2017_A.csv'), ('annotations_2017_B.csv'), ('annotations_2017_C.csv')]
    acc_annot = 0
    for annot in annots:
        with open(op.join(path,annot)) as csvfile:
            acc_annot += pd.read_csv(csvfile) # Accumalitive Dataframe of all 3 expert annotations
    return hdf_files, acc_annot


def load_stevenson_eeghdf(path, selected_eegs):
    """Extracts the stevenson data from the hdf folder and csv files.
    Parameters
    ----------
    path : path to directory that contains the hdf and csv files
    selected_eegs : the range of eegs to be extracted
    Returns
    -------
    eeghdf_files : a list of the h5py files
    acc_annot : a pandas dataframe containing the accumulative
                agreement of the experts
    """
    eeghdf_files = []
    for ii in selected_eegs:
        eeghdf_files.append(eeghdf.Eeghdf('%s/hdf_w_annotations/eeg%s.annot.eeg.h5' % (path, ii)))
    
    annots = [('annotations_2017_A.csv'), ('annotations_2017_B.csv'), ('annotations_2017_C.csv')]
    acc_annot = 0
    for annot in annots:
        with open(op.join(path,annot)) as csvfile:
            acc_annot += pd.read_csv(csvfile) # Accumalitive Dataframe of all 3 expert annotations
    return eeghdf_files, acc_annot
    

def whiten_X(X):
    X = X - np.mean(X, axis = 0) # zero-center the data
    cov = np.dot(X.T, X) / X.shape[0]
    U,S,V = np.linalg.svd(cov)
    X_rot = np.dot(X, U) # decorrelate the data
    X_white = X_rot / np.sqrt(S + 1)
    return X_white

def process_stevenson_data(hdf_files, s_freq, seg_time, n_channels, whiten):
    data = []
    mapping = []
    seg_len = seg_time * s_freq
    for i in range(len(hdf_files)):
        #clear_output(wait=True)
        #print('Loading EEG #%s' % (i))
        rec = (hdf_files[i])['record-0']
        n_segments = (rec['signals'][0].size - 1) // seg_len # Number of data points divided by segment length
        eeg_data = np.empty([n_segments,0,seg_len])
        for ch in range(n_channels): 
            sig = rec['signals'][ch]
            ch_data = [sig[n:n+seg_len] for n in range(0,sig.size,seg_len) # Breaks up eeg into 5 second segments
                       if n < sig.size - seg_len] # Conditional removes remainder not at an even 5 seconds
            ch_data = np.concatenate(ch_data).reshape(n_segments, seg_len)
            if whiten is True: # Whiten channel data
                ch_data = whiten_X(ch_data)
            else:
                ch_data = ch_data - np.mean(ch_data, axis = 0) + 1e-7 # zero-centering channel data
            ch_data = ch_data.reshape(n_segments, 1, seg_len)
            eeg_data = np.concatenate((eeg_data,ch_data),axis = 1)
        data.append(eeg_data)
        filename = hdf_files[i].filename[-12:]
        eeg_mapping = [(filename, n,n+seg_time) for n in range(0, sig.size//s_freq, seg_time)
                       if n < sig.size/s_freq - seg_time]
        mapping.append(eeg_mapping)
    data = np.concatenate(data)
    mapping = np.concatenate(mapping)
    return data, mapping

def process_stevenson_eeghdf_data(eeghdf_files, s_freq, seg_time, n_channels, whiten):
    """based on Lucas's process_stevenson_data but use eeghdf info. He wrote it in C sort of style
    adapting to eeghdf files so for example s_freq is contained
    slices up the files in to equal length segments of length seg_time
    """
    data = []
    mapping = []
    seg_len = seg_time * s_freq
    for ii, hf in enumerate(eeghdf_files):
        #clear_output(wait=True)
        #print('Loading EEG #%s' % (i))
        n_segments = (hf.phys_signals[0].size - 1) // seg_len # Number of data points divided by segment length
        eeg_data = np.empty([n_segments,0,seg_len])
        for ch in range(n_channels): 
            sig = hf.phys_signals[ch]
            ch_data = [sig[n:n+seg_len] for n in range(0,sig.size,seg_len) # Breaks up eeg into 5 second segments
                       if n < sig.size - seg_len] # Conditional removes remainder not at an even 5 seconds
            ch_data = np.concatenate(ch_data).reshape(n_segments, seg_len)
            if whiten is True: # Whiten channel data
                ch_data = whiten_X(ch_data)
            else:
                ch_data = ch_data - np.mean(ch_data, axis = 0) + 1e-7 # zero-centering channel data
            ch_data = ch_data.reshape(n_segments, 1, seg_len)
            eeg_data = np.concatenate((eeg_data,ch_data),axis = 1)
        data.append(eeg_data)
        filename = hf.file_name[-12:]
        eeg_mapping = [(filename, n,n+seg_time) for n in range(0, sig.size//s_freq, seg_time)
                       if n < sig.size/s_freq - seg_time]
        mapping.append(eeg_mapping)
    data = np.concatenate(data)
    mapping = np.concatenate(mapping)
    return data, mapping


def process_stevenson_labels(annot, selected_eegs, sz_agree_crit, seg_time):
    print('Loading Labels')
    y = np.empty([0])
    for i in selected_eegs:
        col = (annot.iloc[:, i-1:i]).dropna() # Grabbing a column and dropping NaN values
        tf_col = (col >= sz_agree_crit).values # Makes T/F seizure vector based on the agreement criteria 
        s_maj = seg_time//2 + 1 # segment parts necessary for majority
        y_seg = [1 if tf_col[t:t+seg_time].sum() >= s_maj else 0 for t in range(0,tf_col.size,seg_time)
                if t < tf_col.size - seg_time]
        y = np.concatenate((y,y_seg))
    return y

def stevenson_eeg_preprocess(path, s_freq = 256, seg_time = 5, selected_eegs = np.arange(1,80),
                             n_channels = 19, sz_agree_criteria = 3, whiten = True):
    """Extracts and preprocesses the stevenon neonatal data.
    Parameters
    ----------
    path : path to directory that contains the hdf and csv files
    s_freq : sampling frequency (default: 256)
    seg_time : duration in seconds of the segments that data, y, and mapping represent (default: 5)
    selected_eegs : the range of eegs to be extracted (default: np.arange(1,80))
    n_channels : the number of channels (default: 19)
    sz_agree_criteria : the number of experts that must agree for a data segment to be marked
                        as a sz (default: 3)
    whiten : If true data is whitened. If false, data is only zero-centered (default: True) 

    Returns
    -------
    data : ndarray, shape (n_segments, n_channels, n_times)
    y : ndarray, shape (n_segments,), 1 for seizure 0 for non seizure
    mapping : ndarray, shape (n_segments, 3), contains filename and segment start and end times
    """
    hdf, annot = load_stevenson_eeg(path, selected_eegs)
    data, mapping = process_stevenson_data(hdf, s_freq, seg_time, n_channels, whiten)
    y = process_stevenson_labels(annot, selected_eegs, sz_agree_criteria, seg_time)
    #clear_output(wait=True)
    print('Loading Complete')
    return data, y, mapping

