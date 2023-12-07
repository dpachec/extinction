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