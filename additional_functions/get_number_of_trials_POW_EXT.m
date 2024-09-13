
function [allTrials2Check] = get_number_of_trials_POW_EXT(EEG, ids1, ids2)

    totalNIds1 = sum(ids1); 
    totalNIds2 = sum(ids2); 

    d2check = EEG.power(ids1, :, :, :); 
    nanTrials = zeros(size(d2check, 1), 1); 
    for triali = 1:size(d2check, 1)
        tP = d2check(triali, :, :, 250:350); 
        x = find(isnan(tP)); 
        if ~isempty(x)
            nanTrials(triali, 1) = 1; 
        end
    end
    if exist('nanTrials')
        allTrials2Check(:, 1) = totalNIds1-sum(nanTrials); 
        %allTrials2Check{subji, 1} = nanTrials; 
    end
    
    d2check = EEG.power(ids2, :, :, :); 
    nanTrials = zeros(size(d2check, 1), 1); 
    for triali = 1:size(d2check, 1)
        tP = d2check(triali, :, :, 250:350); 
        x = find(isnan(tP)); 
        if ~isempty(x)
            nanTrials(triali, 1) = 1; 
        end
    end
    if exist('nanTrials')
        allTrials2Check(:, 2) = totalNIds2 - sum(nanTrials); 
        %allTrials2Check{subji, 2} = nanTrials; 
    end

end