function [h tObs d2pm1 d2pm2 se1 se2]  = compute_real_differences_EXT(out_rsa, t2AV); 


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
        
        d2pm1	= squeeze(mean(cond1,'omitnan'));
        d2pm2	= squeeze(mean(cond2,'omitnan'));
        d2pstd1	= std(cond1, 'omitnan');
        d2pstd2	= std(cond2, 'omitnan');
        disp(['Number of subjects: ' num2str(size(cond1, 1))]); 
        se1 = d2pstd1/sqrt(size(cond1, 1));
        se2 = d2pstd2/sqrt(size(cond1, 1));
        
        [h p ci ts] = ttest(cond1, cond2); 
        h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
            allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        if exist('allSTs')
            [max2u id] = max(abs(allSTs));
            tObs = allSTs(id); 
        else
            tObs = 0; 
        end

    else

        d2pm1	= squeeze(mean(cond1(:, t2AV), 2,'omitnan'));
        d2pm2	= squeeze(mean(cond2(:, t2AV), 2,'omitnan'));
        d2pstd1	= std(mean(cond1(:, t2AV), 2),[], 'omitnan');
        d2pstd2	= std(mean(cond2(:, t2AV), 2),[], 'omitnan');
        disp(['Number of subjects: ' num2str(size(cond1, 1))]); 
        se1 = d2pstd1/sqrt(size(cond1, 1));
        se2 = d2pstd2/sqrt(size(cond1, 1));
        
        [h p ci ts] = ttest(mean(cond1(:, t2AV), 2), mean(cond2(:, t2AV), 2)); 
        h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
            allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        if exist('allSTs')
            [max2u id] = max(abs(allSTs));
            tObs = allSTs(id); 
        else
            tObs = 0; 
        end



    end
    


end