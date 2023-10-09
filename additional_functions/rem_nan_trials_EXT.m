function[EEG] = rem_nan_trials_EXT (EEG)
    for triali = 1:size(EEG.data, 3)
        count = 1; 
        for triali = 1:size(EEG.data, 3)
            dtr = squeeze(EEG.data(:, :,triali)); 
            if find(isnan(dtr))
                id2rem(count,:) = triali; 
                count = count+1;
            end
        end
    end
    if exist('id2rem')
        disp (['removing ' num2str(length(id2rem)) ' trials'])
        EEG.data(:, :, id2rem) = []; 
        EEG.event(id2rem) = []; 
        EEG.nan = length(id2rem); 
    end

end