
function [cond1 cond2] = remove_half_matrix(cond1, cond2)
% % remove half of the matrix
for subji = 1:size(cond1, 1)
    rdm2Tril = squeeze(cond1(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond1TR(subji, :, :) = rdm2Tril;

    rdm2Tril = squeeze(cond2(subji, :, :)); 
    rdm2Tril = tril(rdm2Tril);
    rdm2Tril(rdm2Tril==0) = nan; 
    cond2TR(subji, :, :) = rdm2Tril;

end

cond1 = cond1TR; 
cond2 = cond2TR;