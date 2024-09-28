function [out_rsa, ids] = rem_nan_subj_EXT(out_rsa ,sub2exc, nTrials, minTr)

% % % % %  remove hack 
% % first zero pad for regions with less than 50 subjects

if ~iscell(out_rsa)
    out_rsa(size(out_rsa,1) + 1 : 50, :, :, :) = zeros(length(size(out_rsa,1) + 1 : 50), size(out_rsa, 2), size(out_rsa, 3), size(out_rsa, 4));
    
    out_rsa(sub2exc, :, :, :) = zeros(length(sub2exc),size(out_rsa, 2), size(out_rsa, 3), size(out_rsa, 4));
    ids = []; 
    for subji = 1:size(out_rsa, 1)
        cond1 = squeeze(out_rsa(subji, 1, :, :)); 
        if size(out_rsa, 2) == 2
            cond2 = squeeze(out_rsa(subji, 2, :, :)); 
        end
        if cond1(1) == 0
            ids = [ids subji];
        end
    end
else

    out_rsa(size(out_rsa,1) + 1 : 50) = cell(1,50- size(out_rsa,1) ); 

    %sub2exc2 = find(nTrials(:, 1) <minTr); 
    sub2exc2 = find(nTrials(:, 1) <minTr | nTrials(:, 2) <minTr ); 
    sub2exc3 = [sub2exc, sub2exc2']; 
    sub2exc = unique(sub2exc3)'; 

    out_rsa(sub2exc) =  cell(1,length(sub2exc)); 
    ids = find(cellfun('isempty', out_rsa)); 
    out_rsa(ids) = []; 
    
end