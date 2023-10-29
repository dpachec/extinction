function [EEG] = remove_elec_EXT_manually (EEG, subji)

        if subji ==1
          chans2exc = {'POL F''4'; 'POL F''5'}; 
        end
        if subji ==2
          chans2exc = {'POL H''7'}; 
        end
        if subji == 6
            chans2exc = {'POL A''1'; 'POL A''2'; 'POL A''3'}; 
        end
%         if subji == 7
%             chans2exc = {'POL A''5'}; 
%         end
        if subji == 9
            chans2exc = {'POL F''6'; 'POL F''7'};
        end
         if subji == 10
             chans2exc = {'POL J''1'}; 
         end
%         if subji == 12
%             chans2exc = {'POL A4'; 'POL A5'}; 
%         end
        if subji == 14
            chans2exc = {'POL H''6'; 'POL T''5'; 'POL V''5'}; 
        end
        if subji == 15
            chans2exc = {'POL A9'}; 
        end
        if subji == 16
            chans2exc = {'POL F''6'; 'POL F''3';'POL J''1'; 'POL J''2'}; 
        end
        if subji == 17
            chans2exc = {'POL A''3'; 'POL A''4' }; 
        end
        if subji == 19
            chans2exc = {'POL H''9'};
        end
        if subji == 21
            chans2exc = {'POL H''3'};
        end
        if subji == 23
            chans2exc = {'POL A''1'}; 
        end
        if subji == 24
            chans2exc = {'POL Z''1'; 'POL F''13'; 'POL T''11'}; 
        end
        if subji == 25
            chans2exc = {'POL F''10'; 'POL L''12'}; 
        end
        if subji == 26
            chans2exc = {'POL B''2'}; 
        end
        if subji == 27
            chans2exc = {'POL V4'; 'POL W11'}; 
        end
%         if subji == 28
%             chans2exc = {'POL A4'; 'POL A5'; 'POL A6'; 'POL A7'}; 
%         end
        if subji == 31
            chans2exc = {'AmT2_7'; 'HmT2_6'; 'HpT2_8'; 'TObp_4'; 'TBmd_1'}; 
        end
        if subji == 32
            chans2exc = {'GAd_7'; 'T2md_3'}; 
        end
        if subji == 33
            chans2exc = {'TPod_1'};
        end
        if subji == 38
            chans2exc = {'T2md_4'; 'T2md_5'};
        end
%         if subji == 39
%             chans2exc = {'AmT2_3'}; 
%         end
        if subji == 40
            chans2exc = {'TBp_3'};
        end

%         if subji == 41
%             chans2exc = {'AmT2_4'}; 
%         end   
        if subji == 43
            chans2exc = {'HaT2_6'; 'HaT2_8'; 'TBa_4'; 'TBm_4'; 'TBm_5'}; 
        end   
        if subji == 44
            chans2exc = {'Am2g_8'; 'Ha2d_8'; 'TPod_7'; 'TBad_7'}; 
        end
        if subji == 45
            chans2exc = {'TBm_5'};
        end    
        if subji == 46
            chans2exc = {'TBad_2'; 'TBag_4'; 'TPog_4'};
        end    
        if subji == 47
            chans2exc = {'AmT2_8'}; 
        end    
        
        if exist('chans2exc')
            [chansI ids] = intersect({EEG.chanlocs.labels}', chans2exc); 
            EEG.data(ids, :) = []; 
            EEG.chanlocs(ids) = []; 
            EEG.nbchans = size(EEG.data,1); 
        end










end