
function [ALLEEG] = loadTracesEXT(roi, LT, paths)

    if strcmp(roi, 'AMY')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'AMY' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'AMY' '_C']; 
        end
    elseif strcmp(roi, 'HPC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'HPC' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'HPC' '_C']; 
        end
        
    elseif strcmp(roi, 'OFC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'OFC' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'OFC' '_C']; 
        end
    
    elseif strcmp(roi, 'TMP')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'TMP' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'TMP' '_C']; 
        end
    elseif strcmp(roi, 'OCC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'OCC' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'OCC' '_C']; 
        end
    elseif strcmp(roi, 'FRO')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'FRO' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'FRO' '_C']; 
        end    
    elseif strcmp(roi, 'PFC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'PFC' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'PFC' '_C']; 
        end  
    
    end 
        load ([paths.results.traces file2load]); 

end