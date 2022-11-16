ii = 30


from eeghdf.reader import Eeghdf
#if use 19=ii, may result in length error  bc run into beginning of file- high pass

#need window size 1 second
win_size_sec =1
win_size_samples = int(win_size_sec*eeg.sample_frequency) #will be same number samples in a window
#from 10th to 11th second


def show_ith_window(i, eeg:Eeghdf):
  ax= stacklineplot.stackplot(eeg.rawsignals[:,int((ii*win_size_sec)*eeg.sample_frequency):int(((ii+1)*win_size_sec)*eeg.sample_frequency)],seconds=win_size_sec)
  ax.set_title(f"{ii}th second of eeg file ith window length of "+ f"{win_size_sec}")

show_ith_window(ii, eeg)

def get_window(eeg):
  return eeg.rawsignals[:,int((ii*win_size_sec)*eeg.sample_frequency):int((ii+1)*win_size_sec*eeg.sample_frequency)]
  
index=0
#testing functions in google colab cell.
def common_avg_reference(eeg):
  while ii <= 65:
    eeg_signal_labels = eeg.shortcut_signal_labels[:-2]
    thiswindow = get_window(eeg)

    thiswindow_avg = np.mean(thiswindow[0:-2, :], axis=0) #create new electrode with this avg
    crthiswindow=thiswindow - thiswindow_avg


    ax = stacklineplot.stackplot(crthiswindow)



    #EEG76_1sec.insert(index, crthiswindow)

    #EEG76_marks.insert(index,ii)
    EEG_1sec.insert(index, crthiswindow)

    EEG_marks.insert(index,ii)
    index=index+1
    ii += 5



