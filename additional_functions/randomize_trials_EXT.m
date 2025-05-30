function [out_rsa_perm] = randomize_trials_EXT(out_rsa, cov)



    for subji = 1:size(out_rsa, 1)
        orSc1 = out_rsa{subji}{1}; 
        orSc2 = out_rsa{subji}{2}; 
        
        junts = cat(1, orSc1, orSc2);
        realCondMapping = [zeros(1,size(orSc1, 1)) ones(1, size(orSc2, 1))]';
        fakeCondMapping = realCondMapping(randperm(length(realCondMapping)));

        if strcmp(cov, 'C')
            out_rsa_perm{subji,:}{1} = junts(fakeCondMapping ==0,6:40,6:40);
            out_rsa_perm{subji,:}{2} = junts(fakeCondMapping ==1,6:40,6:40); 
        elseif strcmp(cov, 'V')
            out_rsa_perm{subji,:}{1} = junts(fakeCondMapping ==0,6:45,6:45);
            out_rsa_perm{subji,:}{2} = junts(fakeCondMapping ==1,6:45,6:45); 
        end

    end


end