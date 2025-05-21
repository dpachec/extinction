%% EXPORT TRACES IN LOOP
%% Loads amplitude time series and export epoched data for specific ROIs
clear , clc

c2u = 'C';

listF2sav = {   
                %{'Amygdala' }; 
                %{'Hippocampus'}; 
                %{'orbitofrontal'}; 
                %{'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };
                %{'occipital' '-cuneus' 'lingual' 'pericalcarine'}; % '-' needed for cuneus, to not be confounded with precuneus
                %{'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
                
                    
                % % % COMBINATIONS
                {'Amygdala' 'Hippocampus'}; 
                {'Amygdala' 'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };
                {'Amygdala' 'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
                {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' ...
                 'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
                {'Hippocampus' 'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };
                {'orbitofrontal' 'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' };

            };   
n2SAV = {
            %'AMY'; 
            %'HPC'; 
            %'OFC'; 
            %'TMP'; 
            %'OCC'; 
            %'PFC'; 
            %'PFCO';
            'AMY-HPC';
            'AMY-TMP'; 
            'AMY-PFC';
            'PFC-TMP';
            'HPC-TMP'
            'OFC-TMP'

            
            
            };



allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18', ...
            'p_sub19', 'p_sub20', 'p_sub21'}'; %subject 8 has differnet format (see below)}';


paths = load_paths_EXT; 
for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi paths n2SAV c2u allsubs

    sROI       = listF2sav{listi}; 

    
    
    for subji = 1:length(allsubs)
        subji
        sub = allsubs{subji}; 
        cd([ paths.iEEG])
        load ([sub '_iEEG.mat']);
    
        EEG = remove_elec_EXT_manually(EEG, subji); 
    
        %chansLab = {EEG.chanlocs.fsLabelsR}';
        chansLab = {EEG.chanlocs.fsLabel}';
        selChans = contains(chansLab, sROI);
    
        
        if find(selChans)
            
            EEG.chanlocs = EEG.chanlocs(selChans);
            EEG.data = EEG.data(selChans, :); %contains nans
            
            %epoch data
            Ev = [{EEG.event.type}]'; 
            Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            clen = cellfun(@length, Ev1); 
            EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
            Ev2 = cat(1, Ev1{:});
           
            Ev2(:, 10) = erase(Ev2(:, 10), ' '); %paris subjects have a space in the last character of the event WHY??
            ids = strcmp(Ev2(:, 10), c2u); 
            EEG.event = EEG.event(ids);
    
            %epoch data and markers
            [EEG id2check] = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
            
            EEG = remove_elec_EXT(EEG, 50); %thres channels is 1/5 of 192 = 38
    
    
            if ~isempty(EEG.data)
                nChans(subji, :) = size(EEG.data, 2);
                ALLEEG{subji,:} = EEG; 
            end
        
        end
    
    
    end
    
    mkdir(paths.results.traces)
    filename = [paths.results.traces 'TR_' n2SAV{listi} '_' c2u];
    nSub = sum(cell2mat(cellfun(@(x) ~isempty(x), ALLEEG, 'un', 0)));
    totalChans = sum(nChans);
    save(filename, 'ALLEEG', 'nSub', 'nChans', 'totalChans', '-v7.3');
    
    
    
end
disp('done all files');
cd (paths.github)
