function [ids] = rem_nan_subj_EXT2(out_rsa)

% % % % %  remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)
    cond1 = squeeze(out_rsa(subji, :,1, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end
end
