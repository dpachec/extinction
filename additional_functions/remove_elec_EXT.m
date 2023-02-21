
function [EEG] = remove_elec_EXT(EEG, thresChan)

for chani = 1:size(EEG.data, 1)
    for triali = 1:size(EEG.data, 3)
        d2t = EEG.data(chani, 3001:5000, triali); 
        idN = find(isnan(d2t)); 
        if length(idN)>10
            d2check(triali, chani, :) = 1; 
        else
            d2check(triali, chani, :) = 0; 
        end
    end
end

totalNofNanTrls = sum(d2check);
chans2rem = totalNofNanTrls > thresChan; %% if the total number of trials with nans is higher than this, remove channel

EEG.chanlocs(chans2rem) = []; 
EEG.data(chans2rem, :, :) = []; 



