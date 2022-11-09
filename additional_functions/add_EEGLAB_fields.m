function[EEG] = add_EEGLAB_fields(EEG)
EEG.trials = 1; 
EEG.setname = 'EEG_444702';
EEG.icawinv = [];
EEG.icaweights = []; 
EEG.icasphere =  []; 
EEG.nbchan = size(EEG.data, 1);
EEG.xmax = length(EEG.data); 
EEG.xmin= 0; 
EEG.icaact = []; 
EEG.comments = '';
EEG.pnts = size(EEG.data, 2);

end