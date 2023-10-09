function [phaTS] = extract_pha_EXT(EEG, cfg)
nChans = size(EEG.data, 1); nTrials = size(EEG.data, 3); nTimepoints = size(EEG.data, 2)/10;
for chani = 1:nChans
    for triali = 1:nTrials
        data        = squeeze(EEG.data(chani, :, triali)); % until 1s
        if isempty(find(isnan(data)))
            dataF       = eegfilt(data,1000,cfg.freqs(1),cfg.freqs(end)); 
            dataHA      = angle(hilbert(dataF));
            dataHADS      = downsample(dataHA, 10);
            phaTS(triali, chani, :) = dataHADS; 
        else
            phaTS(triali, chani, :) = nan(1, nTimepoints);
        end
    end
end


end