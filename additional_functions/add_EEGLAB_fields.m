function[EEG] = add_EEGLAB_fields(EEG)

    EEG.trialinfo = ''; 
    EEG.trials = size(EEG.data, 3); 
    EEG.epoch = EEG.trials; 
    EEG.setname = '*********';
    EEG.icawinv = [];
    EEG.icaweights = []; 
    EEG.icasphere =  []; 
    EEG.nbchan = size(EEG.data, 1);
    EEG.xmax = length(EEG.data)/EEG.srate; 
    EEG.xmin= 0; 
    EEG.icaact = []; 
    EEG.comments = '';
    EEG.pnts = size(EEG.data, 2);
    EEG.subject = [];
    EEG.condition = []; 
    EEG.group = [];
    EEG.session = [];
    

end