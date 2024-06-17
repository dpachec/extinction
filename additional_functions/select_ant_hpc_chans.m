function [EEG] = select_ant_hpc_chans(EEG)

    chansCoord = cat(1, EEG.chanlocs.mniCoordR); 
    ids = chansCoord(:, 2) > -22.5; 

    EEG.chanlocs = EEG.chanlocs(ids); 
    EEG.data = EEG.data(ids, :, :); 
    EEG.nbchan = size(EEG.data, 1); 


end