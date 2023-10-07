function [phaTS] = extract_pha_EXT(EEG, cfg)
nChans = size(EEG.data, 1); nTrials = size(EEG.data, 3); 
for chani = 1:nChans
    for triali = 1:nTrials
        data        = squeeze(EEG.data(chani, :, triali)); % until 1s
        if isempty(find(isnan(data)))
            dataF       = eegfilt(data,1000,cfg.freqs(1),cfg.freqs(end)); 
            dataHA      = angle(hilbert(dataF));
            dataHADS      = downsample(dataHA, 10);
            phaTS(triali, chani, :) = dataHADS; 
        end
    end
end


end