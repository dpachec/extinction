function [nTrials] = compute_trial_number_EXT(ids)

    for i = 1:length(ids)
        c1 = ids{i};
        if ~isempty(c1)
            nTrials(i,1) = size(c1{:, 1}, 1);
            nTrials(i,2) = size(c1{:, 2}, 1);
        else
            nTrials(i, 1) = nan; 
            nTrials(i, 2) = nan; 
        end
    
    end

end