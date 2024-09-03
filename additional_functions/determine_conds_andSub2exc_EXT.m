function [cond1 cond2 sub2exc] = determine_conds_andSub2exc_EXT(c1, c2, paths)
    
    cond1 = load ([paths.results.trial_based c1]);
    cond2 = load ([paths.results.trial_based c2]);

    if isfield(cond1, 'ctxTRALL')
        cond1 = cond1.ctxTRALL; 
    else
        cond1 = cond1.itstaTRALL; 
    end
    
    if isfield(cond2, 'ctxTRALL')
        cond2 = cond2.ctxTRALL; 
    else
        cond2 = cond2.itstaTRALL; 
    end
    
    if length(cond1) < 50 
        cond1{50,2}= []; 
    end
    if length(cond2) < 50 
        cond2{50,2}= []; 
    end


    c1s = strsplit(c1, '_'); c2s = strsplit(c2, '_'); 
    if strcmp(c1s{1}, 'trlSTA') & strcmp(c2s{1}, 'trlSTA')
        if (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'TMP')
            sub2exc = [19 21 24];
        elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'TMP')
            sub2exc = [6 19 21 24];
        elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'TMP')
            sub2exc = [19 24 25];
        elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'OFC')) | strcmp(c1s{2}, 'OFC') & strcmp(c2s{2}, 'TMP')
            sub2exc = [];
        elseif (strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')
            sub2exc = [19];        
        elseif (strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'HPC')
            sub2exc = [19];        
        else
            sub2exc = []; 
        end
    elseif strcmp(c1s{1}, 'trlCTX') & strcmp(c2s{1}, 'trlCTX')
        if (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')
            sub2exc = [19];
        elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'PFC')
            sub2exc = []; % no subjects with few trials in this contrast 
        elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'TMP')
            sub2exc = [19 24 25];
        elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'OFC')) | strcmp(c1s{2}, 'OFC') & strcmp(c2s{2}, 'PFC')
            sub2exc = [19];
        else
            sub2exc = []; 
        end
    elseif strcmp(c1s{1}, 'trlCTX') & strcmp(c2s{1}, 'trlSTA')
        if (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')
            sub2exc = [];
        else
            sub2exc = []; 
        end
    end

end