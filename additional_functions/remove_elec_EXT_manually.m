function [EEG] = remove_elec_EXT_manually (EEG, subji)


        if subji == 6
            chans2exc = {'POL A''2'; 'POL A''3'}; 
        end
%         if subji == 7
%             chans2exc = {'POL A''5'}; 
%         end
%         if subji == 11
%             chans2exc = {'POL A5'}; 
%         end
%         if subji == 12
%             chans2exc = {'POL A4'; 'POL A5'}; 
%         end
        if subji == 17
            chans2exc = {'POL A''3'; 'POL A''4' }; 
        end
        if subji == 23
            chans2exc = {'POL A''1'}; 
        end
%         if subji == 27
%             chans2exc = {'POL A3'; 'POL A4'}; 
%         end
%         if subji == 28
%             chans2exc = {'POL A4'; 'POL A5'; 'POL A6'; 'POL A7'}; 
%         end
%         if subji == 39
%             chans2exc = {'AmT2_3'}; 
%         end
%         if subji == 41
%             chans2exc = {'AmT2_4'}; 
%         end   
%         if subji == 44
%             chans2exc = {'Am2g_3'}; 
%         end    
%         if subji == 47
%             chans2exc = {'AmT2_3'}; 
%         end    
        
        if exist('chans2exc')
            [chansI ids] = intersect({EEG.chanlocs.labels}', chans2exc); 
            EEG.data(ids, :) = []; 
            EEG.chanlocs(ids) = []; 
            EEG.nbchans = size(EEG.data,1); 
        end










end