
function [ALLEEG] = loadTracesEXT(roi, LT, paths)

    if strcmp(roi, 'AMY')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'Amygdala' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'Amygdala' '_C']; 
        end
    elseif strcmp(roi, 'HPC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'Hippocampus' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'Hippocampus' '_C']; 
        end
        
    elseif strcmp(roi, 'OFC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'orbitofrontal' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'orbitofrontal' '_C']; 
        end
    
    elseif strcmp(roi, 'VVS')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'VVS' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'VVS' '_C']; 
        end
    elseif strcmp(roi, 'OCC')
        if strcmp(LT, 'V')
            file2load = ['TR_' 'occipital' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['TR_' 'occipital' '_C']; 
        end
    
    end 
        load ([paths.results.traces file2load]); 

end