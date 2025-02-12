function [ids2comp] = extract_ids_EXT(EEG, contrast)

chansLab = {EEG.chanlocs.fsLabel}';

if strcmp(contrast, 'AMY-PFC')
    allPFCLabs = {'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole' };
    chAMYidLeft = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Left')); 
    chHPCidLeft = find(contains(chansLab, allPFCLabs) & contains(chansLab, 'lh')); 
    chAMYidRight = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Right')); 
    chHPCidRight = find(contains(chansLab, allPFCLabs) & contains(chansLab, 'rh')); 
end
if strcmp(contrast, 'AMY-TMP')
    allTMPLabs = {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };
    chAMYidLeft = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Left')); 
    chHPCidLeft = find(contains(chansLab, allTMPLabs) & contains(chansLab, 'lh')); 
    chAMYidRight = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Right')); 
    chHPCidRight = find(contains(chansLab, allTMPLabs) & contains(chansLab, 'rh')); 
end
if strcmp(contrast, 'AMY-HPC')
    chAMYidLeft = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Left')); 
    chHPCidLeft = find(contains(chansLab, 'Hippocampus') & contains(chansLab, 'Left')); 
    chAMYidRight = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Right')); 
    chHPCidRight = find(contains(chansLab, 'Hippocampus') & contains(chansLab, 'Right')); 
end
if strcmp(contrast, 'AMY-OFC')
    chAMYidLeft = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Left')); 
    chHPCidLeft = find(contains(chansLab, 'orbitofrontal') & contains(chansLab, 'Left')); 
    chAMYidRight = find(contains(chansLab, 'Amygdala') & contains(chansLab, 'Right')); 
    chHPCidRight = find(contains(chansLab, 'orbitofrontal') & contains(chansLab, 'Right')); 
end

% % % % only anterior channels 
% mniCoord = cell2mat({EEG.chanlocs.mniCoord}');
% yCoord = mniCoord(:, 2);
% yCond = yCoord > -21; 
% chAMYidLeft = intersect(chAMYidLeft, find(yCond));
% chHPCidLeft = intersect(chHPCidLeft, find(yCond));
% chAMYidRight = intersect(chAMYidRight, find(yCond));
% chHPCidRight = intersect(chHPCidRight, find(yCond));

%if (~isempty(chAMYidLeft) & ~isempty(chHPCidLeft) ) | (~isempty(chAMYidRight) & ~isempty(chHPCidRight) )
if (~isempty(chAMYidLeft) & ~isempty(chHPCidLeft)) | (~isempty(chAMYidRight) & ~isempty(chHPCidRight))
    [A,B] = meshgrid(chAMYidLeft,chHPCidLeft);c=cat(2,A',B');dLeft=reshape(c,[],2);
    [A,B] = meshgrid(chAMYidRight,chHPCidRight);c=cat(2,A',B');dRight=reshape(c,[],2);
    ids2comp = [dLeft; dRight];
else
    ids2comp = []; 
end

end