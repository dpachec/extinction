
function [ALLEEG] = loadTracesEXT(roi, LT, paths)

    if strcmp(roi, 'AMY')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'AMY' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'AMY' '_C_6_4']; 
        end
    elseif strcmp(roi, 'HPC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'HPC' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'HPC' '_C_6_4']; 
        end
        
    elseif strcmp(roi, 'OFC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'OFC' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'OFC' '_C_6_4']; 
        end
    
    elseif strcmp(roi, 'TMP')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'TMP' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'TMP' '_C_6_4']; 
        end
    elseif strcmp(roi, 'OCC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'OCC' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'OCC' '_C_6_4']; 
        end
    elseif strcmp(roi, 'FRO')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'FRO' '_V_6_4']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'FRO' '_C_6_4']; 
        end        
    
    end 
        load ([paths.results.traces file2load]); 

end