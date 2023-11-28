function [EEG] = remove_elec_EXT_manually (EEG, subji)


        if subji ==1
          chans2exc = {'POL F''4'; 'POL F''5'; ; 'POL A''1';'POL V''5';'POL T''10'; 'POL L''14';'POL F''2';'POL F''3'}; 
        end
        if subji ==2
          chans2exc = {'POL H''7'; 'POL X''5'; 'POL B''9'; 'POL X''6'; 'POL X''7'}; 
        end
        if subji ==3
          chans2exc = {'POL L2'; 'POL L3'; 'POL L4'; 'POL L8'; 'POL F11'; 'POL R12'; 'POL R13';'POL H3';'POL C10';'POL H12'}; 
        end
        if subji ==4
          chans2exc = {'POL B1'; 'POL L16'; 'POL S9'; 'POL S7'; 'POL S1'; 'POL B5'; 'POL S5'; 'POL S2'}; 
        end
        if subji ==5
          chans2exc = {'POL J1' ; 'POL J2' ; 'POL J3' ; 'POL A5';'POL V5';'POL V6';'POL L4';'POL V8';'POL L10';'POL L11'};
        end
        if subji == 6
            chans2exc = {'POL A''1'; 'POL A''2'; 'POL A''3'; 'POL H3'; 'POL J1'; 'POL J4';'POL L5';'POL W2';'POL W3';'POL W4';'POL W5';'POL W6';...
                         'POL W7';'POL V5'; 'POL V6'; 'POL H''12';'POL A''8';'POL F''7';'POL A''9'; 'POL A13';}; 
        end
         if subji == 7
             chans2exc = {'POL J''5'; 'POL V''5';'POL V''6'; 'POL X''10'}; 
         end
         if subji == 8
             chans2exc = {'POL H''1'; 'POL Y''1'; 'POL C''2'; 'POL X15'; 'POL X14';'POL R13';'POL R14';'POL Q14';'POL Q12'};
         end
        if subji == 9
            chans2exc = {'POL F''6'; 'POL F''7'; 'POL X''9';'POL J''14';'POL J''13'; 'POL J''8';'POL L''14'; 'POL X''1'; 'POL X''5'}; 
        end
         if subji == 10
             chans2exc = {'POL J''1';'POL WJ2';'POL J1';'POL J2';'POL J3';'POL H11';'POL A''11';'POL H''10'}; 
         end
         
         if subji == 11
             chans2exc = {'POL H1' ; 'POL H6';'POL H9';'POL H10'; 'POL V''3'; 'POL E''2'; 'POL E''3';'POL E''4';'POL Z''6';'POL A10'}; 
         end
         if subji == 12
             chans2exc = {'POL H1'; 'POL H''1'; 'POL H''2';  'POL H''3'; 'POL A5'; 'POL A''1'; 'POL A''2'; 'POL R1'; 'POL V''7';'POL R12'; 'POL E''8'; ...
                          'POL R''14'; 'POL R''13'; 'POL V''10'}; 
         end
        if subji == 13
             chans2exc = {'POL X''6'; 'POL Y''5';'POL L''7' ;'POL T''9';'POL Y''9';'POL T''8'; 'POL A''9'; 'POL H''3'; ...
                            'POL F''2'; 'POL F''3'; 'POL F''4'; 'POL F''5'; 'POL Y''7'; 'POL Y''8'};
        end
        if subji == 14
            chans2exc = {'POL H''6'; 'POL T''5'; 'POL V''5'; 'POL J''6'; 'POL F''4';'POL X''6'}; 
        end
        if subji == 15
            chans2exc = {'POL A9'; 'POL A3'; 'POL A4'; 'POL H1';'POL X1';'POL X2';'POL X3';'POL X4';'POL X5';'POL X6'; 'POL Y4'; 'POL U10'; 'POL O8'; ...
                            'POL S12'; 'POL V4'; 'POL V5'; 'POL V6'; 'POL Z10'}; 
        end
        if subji == 16
            chans2exc = {'POL F''6'; 'POL F''3';'POL J''1'; 'POL J''2'; 'POL A''4';  'POL Z''2'; 'POL Z''3'; 'POL Z''4'; 'POL Z''5' ; ...
                         'POL K''1'; 'POL K''2'; 'POL K''3'; 'POL O''7'; 'POL K''7';'POL J''3';'POL J''4';'POL A''11';'POL F''10'; ...
                         'POL Y''8';'POL J''5';'POL J''9'; 'POL H''5'}; 
        end
        if subji == 17
            chans2exc = {'POL A''3'; 'POL A''4';'POL R''13';'POL R''14';'POL R''15'; 'POL V''9'}; 
        end
        if subji == 18
            chans2exc = {'POL J5'; 'POL H1'; 'POL N5'; 'POL N6';'POL J11';'POL J12';'POL H10';'POL V5';'POL V7';'POL R13';'POL F12'}; 
        end
        if subji == 19
            chans2exc = {'POL H''9'; 'POL Z''6'; 'POL Z''7';'POL F''11'}; 
        end
        if subji == 20
            chans2exc = {'POL B''1'; 'POL F''1'; 'POL Y''17'; 'POL B''13'; 'POL P''1'; ...
                         'POL X''11'; 'POL X''16'; 'POL X''17'; 'POL L''2'; 'POL F''2';'POL B''11'; 'POL F''7'}; 
        end
        if subji == 21
            chans2exc = {'POL H''3' ; 'POL J''3'; 'POL A''5'; 'POL P''10'; 'POL V''4';'POL Y''8';'POL J''11'}; 
        end
        if subji == 22
            chans2exc = {'POL H1'; 'POL V6';'POL V''5'}; 
        end
        if subji == 23
            chans2exc = {'POL A''1'; 'POL J''1'; 'POL J''2'; 'POL J''3'; 'POL Z''1'; 'POL V''8'; 'POL V''6'; 'POL T''4'}; 
        end
        if subji == 24
            chans2exc = {'POL Z''1'; 'POL F''13'; 'POL T''11'; 'POL A''1'; 'POL R''5'; 'POL X''4'; 'POL V''7'; 'POL Y''4';'POL F''1';'POL H''3'}; 
        end
        if subji == 25
            chans2exc = {'POL F''10'; 'POL L''12'; 'POL X''5'; 'POL X''6'; 'POL V''5'; 'POL T''4';'POL W''6';'POL V''9'}; 
        end
        if subji == 26
            chans2exc = {'POL B''2'; 'POL R''11'; 'POL X''4'; 'POL O''4'; 'POL B''8';'POL R''9';'POL R''10';'POL W''7';'POL T''7'; 'POL J''8'; 'POL Z''9'; 'POL Z''10'}; 
        end

        if subji == 27
            chans2exc = {'POL V4'; 'POL W11'; 'POL A4'; 'POL Z1'; 'POL Z2'; 'POL Z7'; 'POL Z10'; 'POL Z11'; 'POL H4'; 'POL V''6'; 'POL W5'; ...
                         'POL J11'; 'POL J12'; 'POL H''4';'POL H''8';'POL H9'; 'POL H6'; 'POL W6'; 'POL W8'; 'POL H''5'; 'POL H''7';'POL W9'};
        end
         if subji == 28
             chans2exc = {'POL A7'; 'POL N2' ; 'POL N3' ; 'POL N4' ; 'POL N5' ; 'POL N6' ; 'POL N7'; 'POL Z15'; 'POL B10'; 'POL F14'; ...
                            'POL K15'; 'POL X11';'POL Y13'; 'POL B9'; 'POL R9';'POL A11'}; 
         end
        if subji == 29
            chans2exc = {'POL H3' ; 'POL H''4'; 'POL A4'; 'POL X6';'POL W14';'POL W15';'POL J10'; 'POL V''6'; ; 'POL A11';'POL H2';'POL J3'};
        end
        if subji == 30
            chans2exc = {'POL J4' ; 'POL J5'; 'POL J11';'POL V5';'POL F11';'POL A11'};
        end
        if subji == 31
            chans2exc = {'AmT2_7'; 'HmT2_6'; 'HpT2_8'; 'TObp_4'; 'TBmd_1'; 'OLa_1'; 'HpT2_4';'HmT2_7';'HmT2_4'}; 
        end
        if subji == 32
            chans2exc = {'GAd_7'; 'T2md_3';'TBmd_1';'Hm2g_5';'Hm2g_4';'Hm2g_5';'T2md_2';'TBpd_5'}; 
        end
        if subji == 33
            chans2exc = {'TPod_1'; 'HaT2_3'; 'AmT2_1'; 'AmT2_3'; 'HmT2_2';'HaT2_2'};
        end
        if subji == 34
            chans2exc = {'TBa_1';'LMI2_8';'TBa_7';'TBp_2'};
        end
        if subji == 35
            chans2exc = {'HaT2_3'; 'AmT2_2';'F2a_1'; 'F2a_3'; 'F3a_6';'PMs_5';'F2a_2'; 'F2p_2';'F2p_1'; 'F1a_2'; 'F1p_2';'HaT2_7'};

        end
        if subji == 36
            chans2exc = {'HaT1_1'; 'HaT1_2'; 'InAm_2'; 'InaF_4'; 'TPol_2';'HmT2_3';'HmT2_4';'TBm_1';'TBa_2'};
        end
        if subji == 37
            chans2exc = {'F1m_4'; 'F1p_1'; 'F1a_2'; 'F1p_2'; 'F1m_2'; 'F1m_3'; };
        end
        if subji == 38
            chans2exc = {'T2md_4'; 'T2md_5'; 'AmT2_2'; 'AmT2_4'; 'T2mg_1';'TPod_2';'TPod_6';'T2Bg_2';'T2pg_8';'TPod_5';'TPog_3'};
        end
        if subji == 39
            chans2exc = {'HaT1_3';'TBa_1'; 'TBa_3'; 'TBp_2'}; 
        end
        if subji == 40
            chans2exc = {'TBp_3';'TBp_6'; 'TBa_6'; 'TBm_5'};
        end

         if subji == 41
             chans2exc = {'AmT2_4'; 'TBm_3';'InT1_1';'InT1_2';'TBp_5';'TBm_6';'TBm_4'}; 
         end   
        if subji == 43
            chans2exc = {'HaT2_6'; 'HaT2_8'; 'TBa_4'; 'TBm_4'; 'TBm_5'; 'HaT2_2';'HaT2_3';'Am2d_4' ;'HaT2_1'; 'Ha1d_3';'T1a_1';'TBa_2'; 'HmT1_8'; 'TBa_1'; 'TBa_3' }; 
        end   
        if subji == 44
            chans2exc = {'Am2g_8'; 'Ha2d_8'; 'TPod_7'; 'TBad_7'; 'Ha2g_3'; 'Ha2d_2'; 'Ha2d_3'; 'Ha2d_4'; 'Am2g_3';'TPog_1';'T1pd_2';'T1pd_1';'Ha2d_7';'TBag_6'; ...
                         'TBad_6'; 'TBag_3'; 'TBag_4'}; 
        end
        if subji == 45
            chans2exc = {'TBm_5'; 'AmT2_3';'HaT2_5';'TBm_1';'HmT2_3';'AmT2_6'; 'TBa_2'; 'TBa_4'; 'TBa_5'; 'TBa_6'};
        end    
        if subji == 46
            chans2exc = {'TBad_2'; 'TBag_4'; 'TPog_4'; 'TBmg_2'; 'Ha2d_1'; 'Ha2d_2';'Ha2g_2';'Am2g_2';'HLe2_3';'TBag_3';'HLe2_6'; ...
                        'TBad_4'; 'TBad_5'; 'TPog_5'; 'TPog_6'};
        end    
        if subji == 47
            chans2exc = {'AmT2_8'; 'HaT1_4';'TBa_2';'TBa_5'}; 
        end   
        if subji == 48
            chans2exc = {'TBa_3'; 'HmT2_4';'TBa_7'; 'TBm_1'; 'TBm_2'}; 
        end   
        if subji == 49
            chans2exc = {'HmT2_4';'TBp_7';'TBp_8'}; 
        end   
        if subji == 50
            chans2exc = {'TBa_1';'O3a_1';'TBp_3';'O3p_3'; 'TBp_1'; 'TBp_2'}; 
        end  
        
        
        if exist('chans2exc')
            [chansI ids] = intersect({EEG.chanlocs.labels}', chans2exc); 
            EEG.data(ids, :) = []; 
            EEG.chanlocs(ids) = []; 
            EEG.nbchans = size(EEG.data,1); 
        end










end