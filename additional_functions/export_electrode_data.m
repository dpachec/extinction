%% 
%%
clear, clc

paths = load_paths_EXT; 


%files2load = {['TR_' 'OFC' '_C']; ['TR_' 'AMY' '_C']; ['TR_' 'HPC' '_C']; ['TR_' 'TMP' '_C']; ['TR_' 'OCC' '_C']; ['TR_' 'FRO' '_C']; }
files2load = {['TR_' 'HPC' '_C']; };

allElec = []; 
for filei = 1:length(files2load)

    load ([paths.results.traces files2load{filei}]); 

    ROIp = strsplit(files2load{filei}, '_');
    ROI = ROIp{2}; 

    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji}; 
        if ~isempty(EEG)
            chSub = EEG.chanlocs; 
            chSub2 = squeeze(struct2cell(chSub))';
            for i = 1:size(chSub2, 1)
                chSub2{i, 7} = subji; 
                chSub2{i, 8} = ROI; 
            end
            allElec = [allElec; chSub2];
        end
    end

end

%%create one CSV file with all electrodes
writecell(allElec, 'allElec.csv')






%% % % % % % % % FIRST LOAD FILE - > START HERE
clear, clc

paths = load_paths_EXT; 
file2load = ['TR_' 'OCC' '_C']; 
load ([paths.results.traces file2load]); 

%% create one CSV file with all electrodes in one region
clearvars -except ALLEEG
clc

allElec = []; 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG)
        chSub = EEG.chanlocs; 
        %chSub2 = struct2cell(chSub)';
        chSub2 = squeeze(struct2cell(chSub))';
        for i = 1:size(chSub2, 1)
            chSub2{i, 7} = subji; 
        end
        allElec = [allElec; chSub2];
    end
end

writecell(allElec, 'allElec.csv')

%% create one CSV file with all electrodes in ALL REGIONS


%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast
clear , clc

listF2sav = {
                ['TR_' 'OCC' '_C']; 
                ['TR_' 'PFC' '_C']; 
                ['TR_' 'OFC' '_C']; 
                ['TR_' 'TMP' '_C']; 
                ['TR_' 'AMY' '_C']; 
                ['TR_' 'HPC' '_C']; 
};   


t1 = datetime; 
allElec = []; 
paths = load_paths_EXT; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 ALLEEG allElec paths
    f2sav       = listF2sav{listi}; 
    load ([paths.results.traces f2sav])
    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji}; 
        if ~isempty(EEG)
            chSub = EEG.chanlocs; 
            %chSub2 = struct2cell(chSub)';
            chSub2 = squeeze(struct2cell(chSub))';
            for i = 1:size(chSub2, 1)
                chSub2{i, 7} = subji; 
            end
            allElec = [allElec; chSub2];
        end
    end
end


writecell(allElec, 'allElec.csv')

%% create one CSV file with all electrodes for AMY-HPC
clear, clc
paths = load_paths_EXT; 
file2load = ['TR_' 'AMY-HPC' '_C']; 
load ([paths.results.traces file2load]); 
clearvars -except ALLEEG
clc

allElec = []; 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji}; 
    
    if ~isempty(EEG) 
        chansLab = {EEG.chanlocs.fsLabel}';
        chAMYid = contains(chansLab, 'Amygdala'); 
        chHPCid = contains(chansLab, 'Hippocampus'); 
    
        if ~isempty(chAMYid) & ~isempty(chHPCid) 
            chSub = EEG.chanlocs; 
            %chSub2 = struct2cell(chSub)';
            chSub2 = squeeze(struct2cell(chSub))';
            for i = 1:size(chSub2, 1)
                chSub2{i, 7} = subji; 
            end
            allElec = [allElec; chSub2];
        end
    end
end

writecell(allElec, 'allElecAMY-HPC.csv')


%%
count = 0; 
for subji = 1:47
    if ~isempty(ALLEEG_AMY{subji}) & ~isempty(ALLEEG_HPC{subji})
        count = count+1; 
    end

end