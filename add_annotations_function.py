

# %%

def add_annotations(eeg_num):
    """this function was designed to load hdf5 files from the stevenson EEG datatset
    and add annotations
    """
     
      # %%
      # run through all [1...79] inclusive eegs
      #for eeg_num in range(start_eeg_num,end_eeg_num+1):
          #print(eeg_num)

      # %%
      # for eeg_num in eeg_columns_num:
      #     print(f"{eeg_num=}")
      #     dfa75, dfb75, dfc75 = dfa[str(eeg_num)].dropna(), dfb[str(eeg_num)].dropna(), dfc[str(eeg_num)].dropna()


      dfaeeg, dfbeeg, dfceeg = dfa[str(eeg_num)].dropna(), dfb[str(eeg_num)].dropna(), dfc[str(eeg_num)].dropna()    

     

      
      # %%
      eegannot_df = pd.concat([dfaeeg, dfbeeg, dfceeg], axis=1)

      # %%
      eegannot_df.columns = ['A','B','C']

      # %%


      # %%
      eegannotmean = eegannot_df.mean(axis=1)



      # %%
      eegannot_df['consensus'] = eegannotmean > 0.5 

      # %%
      
      # %%
      eegannot_df['consensus'] = eegannot_df.consensus.apply(lambda tf: 1 if tf else 0)

      # %%




      # %%
      # test out what iloc does - it gets a row of values as a series Dataframe.index for our simple index does nothing
      # df.index[ii] -> ii

      s = eegannot_df.iloc[20]

    

      # %%
      bigarr = eegannot_df.to_numpy()

      # %%
     

      # %%
      eeg_consensus_label = bigarr[:,3]

      # %%
      # checking things out
      # # where does the label say there are seizures
      np.where(eeg_consensus_label == 1)

      # %%
      # dict of all labels
      eeg_labels = {eeg_num: eeg_consensus_label}

 
      # %%



      # %%



      shutil.copy(eeg_file_path+f"eeg{eeg_num}.eeg.h5", eeg_annot_path+f"eeg{eeg_num}.annot.eeg.h5") # !cp eeg{eeg_num}.eeg.h5 eeg{eeg_num}.annot.eeg.h5

      # %%
      hf = h5py.File(eeg_annot_path+f'eeg{eeg_num}.annot.eeg.h5', 'r+')
      

      # %%
      

      # %%
      rec = hf['record-0']

     

      # %%
      # make sure record datetimes are in isoforamt
      rec.attrs['end_isodatetime'] = dateutil.parser.parse(rec.attrs['end_isodatetime']).isoformat()
      rec.attrs['end_isodatetime']

      # %%
      rec.attrs['start_isodatetime'] = dateutil.parser.parse(rec.attrs['start_isodatetime']).isoformat()
      rec.attrs['start_isodatetime']

      # %%


      # %%
      rec.create_group("extensions")
      recext = rec['extensions']

      AN_GROUP_NAME= 'dense_seizure_annotations'
      recext.create_group(AN_GROUP_NAME)

     

      # %%

      annots = recext[AN_GROUP_NAME]

      # %%
      AN_GROUP_DESCRIPTION = \
      """
      annotations from Stevenson et al, neonatal dataset for annotators A, B, C labeling each second with a 1 or 0 to signal
      a seizure or non-seizure during that period of time. annotarr.shape = (n_annotators,n_seconds) where 3 is the number of
      annotators. thus annotarr[2,200] is an integer (a zero) which indicates what annotator C marked at the 201st second of 
      the study as not having a seizure.
      """

      annots.attrs['description'] = AN_GROUP_DESCRIPTION

      # %%
      abc_arr = bigarr[:,:3].T

      # %%
      

      # %%
      annots.create_dataset("annotarr", data=abc_arr, dtype=np.int8)

      # %%
      dset = annots['annotarr']

      # %%
      dset[2,200:210]

      # %%
      avg_abc = np.mean(dset[:,:], axis=0)

      # %%
      avg_abc.dtype, avg_abc.shape

      # %%
      

      # %%
      cons = np.zeros_like(avg_abc, dtype='int8')
      print("cons.dtype=", cons.dtype, "cons.shape=", cons.shape)

      cons[np.where(avg_abc > 0.5)] = 1

      # %%
      

      # %%
      

      # %%
      annots.create_dataset("consensus_sz_labels", data=cons, dtype=np.int8)

      # %%
      list(hf.attrs.keys())

      # %%
      hf.attrs['EEGHDFversion']

      # %%
      rec = hf['record-0']

      # %%
      list(rec.attrs.items())

      # %%
      rec.attrs['studyadmincode']

      # %%
      rec.attrs['studyadmincode'] = str(eeg_num)

      # %%
      rec.attrs['studyadmincode']

      # %%
      pt = hf['patient']
      list(pt.attrs.items())

      # %%
      patient_clin = clinical.iloc[eeg_num-1] # remember subject and eeg numbering starts at 1
      patient_clin

      # %%
      pt.attrs['patientcode'] = str(patient_clin.ID)
      pt.attrs['patient_name'] = f"{patient_clin.ID}, Subject"
      pt.attrs['patientcode'],  pt.attrs['patient_name']

      # %%
      pt.attrs['gender'] = patient_clin.Gender.upper() if type(patient_clin.Gender) == str else 'U' 
      pt.attrs['gender']

      # %%
      # pt.attrs['EEG file'] = patient_clin['EEG file']
      # pt.attrs['EEG file']
      ga_str = patient_clin['GA (weeks)'] # note that '' are translated into float nan
      if type(ga_str) == str:

          wks_list = [float(xx) for xx in ga_str.split() if xx.isnumeric()]
          wks = np.mean(wks_list)
          wks*7 # convert to days
          pt.attrs['gestatational_age_at_birth_days'] = wks*7
          if wks < 37:
              pt.attrs['born_premature'] = "true"
          if wks >= 37:
              pt.attrs['born_premature'] = "false"
      else: # not available information
          pt.attrs['gestatational_age_at_birth_days'] = -1.0
          pt.attrs['born_premature'] = "unknown"
      

      # %%
      pt.create_group("extensions")
      ptext = pt['extensions']

      # %%
      ptext

      # %%

      auto_add = [ 'EEG file', 'Gender', 'BW (g)', 'GA (weeks)',
          'EEG to PMA (weeks)', 'Diagnosis', 'Neuroimaging Findings',
          'PNA of Imaging (days)', 'Number of Reviewers Annotating Seizure',
          'Primary Localisation', 'Other']

      auto_add = patient_clin.index.to_list()

      for ind in auto_add:
          val = str(patient_clin[ind])
          if val == 'nan':
              ptext.attrs[ind] = ''
          else:
              ptext.attrs[ind] = val





      # %%
      patient_clin.index.to_list()

      # %%
      list(ptext.attrs.items())

      # %%
      list(rec.items())

      # %%
      list(rec.attrs.items())

      # %%

      # %%
      # fix datetimes

      # %%


      # %%
      # finish up
      hf.close()

      # %%



