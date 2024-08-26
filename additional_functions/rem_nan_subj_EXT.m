function [ids] = rem_nan_subj_EXT(out_rsa ,sub2exc)

% % % % %  remove hack 
%out_rsa(sub2exc, :, :, :) = zeros(1,size(out_rsa, 2), size(out_rsa, 3), size(out_rsa, 4));
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
