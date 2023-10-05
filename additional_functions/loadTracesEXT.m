
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
        
    elseif strcmp(roi, 'ORB')
        if strcmp(LT, 'V')
            file2load = ['allS_' 'orbitofrontal' '_V']; 
        elseif strcmp(LT, 'C')
            file2load = ['allS_' 'orbitofrontal' '_C']; 
        end
        
    end
        %file2load = ['allS_' 'superiorfrontal' '_C']; 
        
        load ([paths.results.traces file2load]); 

end