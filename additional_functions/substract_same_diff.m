function [rsaDiff]  = substract_same_diff(out_rsa, t2AV); 


    cond1 = cellfun(@mean, cellfun(@(x) x{1}, out_rsa, 'un', 0), 'un', 0); 
    cond2 = cellfun(@mean, cellfun(@(x) x{2}, out_rsa, 'un', 0), 'un', 0); 
    
    % check if nans
    clength1 = double(string(cellfun(@length, cond1, 'un', 0)));
    clength2 = double(string(cellfun(@length, cond2, 'un', 0)));
    nanSubj1 = clength1 == 1; 
    nanSubj2 = clength2 == 1; 

    nanSubjPrev = nanSubj1 | nanSubj2; 
    nanSubj = find(nanSubjPrev); 


    cond1(nanSubj) = []; 
    cond2(nanSubj) = []; 

    cond1 = cat(1, cond1{:}); 
    cond2 = cat(1, cond2{:}); 
    
    
    for subji = 1:size(cond1, 1)
       cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
       cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
    end       
    cond1 = cond1B(:, 1:40); cond2 = cond2B(:, 1:40); 
    
    

    



    if isempty(t2AV)
        
      
    else

      rsaDiff = mean(cond1(:, t2AV), 2) - mean(cond2(:, t2AV), 2); 

    end
    


end