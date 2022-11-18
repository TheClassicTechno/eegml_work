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
