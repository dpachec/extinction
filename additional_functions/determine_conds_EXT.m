function [cond1 cond2] = determine_conds_EXT(c1, c2, paths)
    

    sub2exc = [11 37]';


    cond1 = load ([paths.results.trial_based c1]);
    if ~strcmp(c2, 'none')
        cond2 = load ([paths.results.trial_based c2]);
    else 
        cond2 = struct; 
    end
    

    if isfield(cond1, 'ctxTRALL')
        cond1 = cond1.ctxTRALL; 
    else
        if isfield(cond1, 'itstaTRALL')
            cond1 = cond1.itstaTRALL; 
        elseif isfield(cond1, 'trlSTA_ET')
                cond1 = cond1.trlSTA_ET; 
                disp('cond1 trlSTA_ET')
        elseif isfield(cond1, 'trlSTA_TE')
                cond1 = cond1.trlSTA_TE; 
                disp('cond1 trlSTA_TE')
        elseif isfield(cond1, 'trlSTA_TA')
                cond1 = cond1.trlSTA_TA; 
                disp('cond1 trlSTA_TA')
        end
    end 
    
    
    if isfield(cond2, 'ctxTRALL')
        cond2 = cond2.ctxTRALL; 
    else
        if isfield(cond2, 'itstaTRALL')
            cond2 = cond2.itstaTRALL; 
        elseif isfield(cond2, 'trlSTA_ET')
            cond2 = cond2.trlSTA_ET; 
        elseif isfield(cond2, 'trlSTA_AT')
            cond2 = cond2.trlSTA_AT; 
        else
            cond2 = cell(1); 

        end
    end
    
    if length(cond1) < 50 
        cond1{50,2}= []; 
    end
    if length(cond2) < 50 
        cond2{50,2}= []; 
    end


    c1s = strsplit(c1, '_');
    c2s = strsplit(c2, '_'); 


    % % % 
    % % % if strcmp(c1s{1}, 'trlSTA') & strcmp(c2s{1}, 'trlSTA')
    % % %     if (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [19 21 24 27 37];
    % % %     elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [6 19 21 24 27 37];
    % % %     elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [19 24 25 27 37];
    % % %     elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'OFC')) | strcmp(c1s{2}, 'OFC') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [27 37];
    % % %     elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'OCC')) | strcmp(c1s{2}, 'OCC') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [25 27 37];            
    % % %     elseif (strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')
    % % %         sub2exc = [];
    % % %     else
    % % %         sub2exc = [];
    % % %     end
    % % % elseif strcmp(c1s{1}, 'trlCTX') & strcmp(c2s{1}, 'trlCTX')
    % % %     if (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = [19];
    % % %     elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = []; % no subjects with few trials in this contrast 
    % % %     elseif (strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'PFC')) | strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'TMP')
    % % %         sub2exc = [19 24 25];
    % % %     elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'OFC')) | strcmp(c1s{2}, 'OFC') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = [19];
    % % %     else
    % % %         sub2exc = []; 
    % % %     end
    % % % elseif (strcmp(c1s{1}, 'trlCTX') & strcmp(c2s{1}, 'trlSTA')) |  (strcmp(c2s{1}, 'trlCTX') & strcmp(c1s{1}, 'trlSTA'))
    % % %     if (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'AMY')) | strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = [19 27 37];
    % % %     elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'TMP')) | strcmp(c1s{2}, 'TMP') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = [24 25 29 27 37]; 
    % % %     elseif (strcmp(c1s{2}, 'PFC') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'PFC')
    % % %         sub2exc = [27 37]; 
    % % %     elseif (strcmp(c1s{2}, 'AMY') & strcmp(c2s{2}, 'HPC')) | strcmp(c1s{2}, 'HPC') & strcmp(c2s{2}, 'AMY')
    % % %         %sub2exc = [11 14 21 27 29]; %less than 10 trials
    % % %         sub2exc = [27 37]; 
    % % %     else
    % % %         sub2exc = []; 
    % % %     end
    % % % elseif strcmp(c1s{1}, 'trlSTA') 
    % % %     if strcmp(c1s{2}, 'PFC')
    % % %         sub2exc = [19 27];
    % % %     elseif strcmp(c1s{2}, 'TMP')
    % % %         sub2exc = [27 37]; 
    % % %     elseif strcmp(c1s{2}, 'AMY')
    % % %         sub2exc = [27 37];             
    % % %     elseif strcmp(c1s{2}, 'HPC')
    % % %         sub2exc = [27 37];  
    % % %     elseif strcmp(c1s{2}, 'OFC')
    % % %         sub2exc = [27 ];              
    % % %     end
    % % % 
    % % % end
    % % % 
    % % % 


    

end