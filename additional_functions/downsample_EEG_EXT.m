function [EEG] = downsample_EEG_EXT(EEG)


    for chani = 1:size(EEG.data, 1)
        data_prev = squeeze(EEG.data(chani, :, :));
        data(chani, :,:) = downsample(data_prev, 10);
    end
    
    EEG.data = data; 


end