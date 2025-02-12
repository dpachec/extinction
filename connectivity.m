%%
%% load traces and filter

clear

listF2sav = {

'TR_AMY-PFC_C_4-8'; 
'TR_AMY-HPC_C_4-8'; 
'TR_AMY-TMP_C_4-8'; 
'TR_PFC-TMP_C_4-8'; 
'TR_HPC-TMP_C_4-8'; 
'TR_OFC-TMP_C_4-8'; 
};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
        
    f2sav       = listF2sav{listi}; 
    f2uBef = strsplit(f2sav, '_'); 
    f2uBef2 = f2uBef{end};  
    f2u = [double(string(f2uBef2(1))) double(string(f2uBef2(3)))];

    f2load = char(string(join(f2uBef(1:3), '_'))); 
    load ([paths.results.traces f2load]); 
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);


        
        if ~isempty(ALLEEG{subji})
            EEG = ALLEEG{subji}; 
            Ev = extract_event_EXT(EEG);
            data= EEG.data;
            % nans should be cleaned as these are not included 
            allIDs2exc = []; 
            for triali = 1:size(data, 3)
                dataTR = squeeze(data(:, :, triali)); 
                if ~isempty(find(isnan(dataTR))) 
                    allIDs2exc = [allIDs2exc triali]; 
                
                end

            end
            
            data(:, :,allIDs2exc ) = []; 
            Ev(allIDs2exc, :) = []; 
            

            % % % Filter data
             for chani = 1:size(data, 1)
                 for triali = 1:size(data, 3)
                    amyD = squeeze(data(chani, :, triali));
                    amyDF = eegfilt (amyD,1000, f2u(1), f2u(end));
                    dataF(chani, :, triali) = amyDF; 
                 end
             end

             ALLFILT{subji, 1} = dataF; 
             ALLFILT{subji, 2} = Ev; 

        end
    end

    save([ paths.results.rsa f2sav 'Hz.mat'], 'ALLFILT');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%%

clear 
f2load = 'TR_AMY-HPC_C_4-8Hz'; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = 'TR_AMY-HPC_C'; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 

%%Connectivity across time for every trial
clc
clearvars -except ALLFILT ALLEEG file2load paths f2load

foi = [4:8];
toi = [3001:4750]; % full trial 


for subji = 1:length(ALLFILT)
    disp(['Subji: ' num2str(subji)])

    data    = ALLFILT{subji, 1};
    Ev      = ALLFILT{subji, 2};
    EEG     = ALLEEG{subji}; 

    clear ids2comp
    if ~isempty(data)
        cbH = strsplit(f2load, '_'); 
        contrast = cbH{2}; 
        
        ids2u = Ev(:, 2) == 2; 
        data = data(:, :, ids2u);% only during extinction
        
        ids2comp = extract_ids_EXT(EEG, contrast);
        
        if ~isempty(ids2comp) 
            CON_B = compute_connectivity_EXT(data, toi, [], foi, ids2comp, 'Trials'); % only time in this band
            %size(CON_B, 1)
            CON(subji, :) =  mean(CON_B, 'all');
            %CON{subji, 1} = mean(CON_B, 'all');
            %CON{subji, 2} = Ev(ids2u, :); 
        end
   end
        
end

CON(CON==0) = []; 
%save(['CON_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end))], 'CON');



%%permutations

nPerm = 1000; 
clear CON_P; 
clc 

for permi = 1:nPerm
    disp(['Permi: ' num2str(permi)])
    clear CON_BB
    for subji = 1:length(ALLFILT)
        clear dataPerm; 
        
        data    = ALLFILT{subji, 1};
        Ev      = ALLFILT{subji, 2};
        EEG     = ALLEEG{subji}; 
        
        clear ids2comp
        if ~isempty(data)
            for t = 1:size(data, 3)
                % Apply circshift along the time dimension (the 2nd dimension).
                shiftBy =round(size(data, 2)*rand());
                dataPerm(:, :, t) = circshift(data(:, :, t), [0, shiftBy, 0]);
            end
            
            cbH = strsplit(f2load, '_'); 
            contrast = cbH{2}; 
            
            ids2u = Ev(:, 2) == 2; 
            dataPerm = dataPerm(:, :, ids2u);% only during extinction
            
            ids2comp = extract_ids_EXT(EEG, contrast);
            
            if ~isempty(ids2comp) 
                CON_B = compute_connectivity_EXT(dataPerm, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                CON_BB(subji, :) =  mean(CON_B, 'all');
            end
        end
    end
    %size(CON_BB)
    CON_BB(CON_BB==0) = []; 
    %size(CON_BB)
    CON_P(:, permi, :) = CON_BB; 
end



save(['CON_P_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) '_' num2str(nPerm), '_p' ], 'CON_P');









clear 
f2load = 'TR_AMY-TMP_C_4-8Hz'; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = 'TR_AMY-TMP_C'; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 

%%Connectivity across time for every trial
clc
clearvars -except ALLFILT ALLEEG file2load paths f2load

foi = [4:8];
toi = [3001:4750]; % full trial 


for subji = 1:length(ALLFILT)
    disp(['Subji: ' num2str(subji)])

    data    = ALLFILT{subji, 1};
    Ev      = ALLFILT{subji, 2};
    EEG     = ALLEEG{subji}; 

    clear ids2comp
    if ~isempty(data)
        cbH = strsplit(f2load, '_'); 
        contrast = cbH{2}; 
        
        ids2u = Ev(:, 2) == 2; 
        data = data(:, :, ids2u);% only during extinction
        
        ids2comp = extract_ids_EXT(EEG, contrast);
        
        if ~isempty(ids2comp) 
            CON_B = compute_connectivity_EXT(data, toi, [], foi, ids2comp, 'Trials'); % only time in this band
            %size(CON_B, 1)
            CON(subji, :) =  mean(CON_B, 'all');
            %CON{subji, 1} = mean(CON_B, 'all');
            %CON{subji, 2} = Ev(ids2u, :); 
        end
   end
        
end

CON(CON==0) = []; 
%save(['CON_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end))], 'CON');



%%permutations

nPerm = 1000; 
clear CON_P; clear dataPerm; 
clc 

for permi = 1:nPerm
    disp(['Permi: ' num2str(permi)])
    clear CON_BB
    for subji = 1:length(ALLFILT)
        clear dataPerm; 
        
        data    = ALLFILT{subji, 1};
        Ev      = ALLFILT{subji, 2};
        EEG     = ALLEEG{subji}; 
        
        clear ids2comp
        if ~isempty(data)
            for t = 1:size(data, 3)
                % Apply circshift along the time dimension (the 2nd dimension).
                shiftBy =round(size(data, 2)*rand());
                dataPerm(:, :, t) = circshift(data(:, :, t), [0, shiftBy, 0]);
            end
            
            cbH = strsplit(f2load, '_'); 
            contrast = cbH{2}; 
            
            ids2u = Ev(:, 2) == 2; 
            dataPerm = dataPerm(:, :, ids2u);% only during extinction
            
            ids2comp = extract_ids_EXT(EEG, contrast);
            
            if ~isempty(ids2comp) 
                CON_B = compute_connectivity_EXT(dataPerm, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                CON_BB(subji, :) =  mean(CON_B, 'all');
            end
        end
    end
    %size(CON_BB)
    CON_BB(CON_BB==0) = []; 
    %size(CON_BB)
    CON_P(:, permi, :) = CON_BB; 
end



save(['CON_P_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) '_' num2str(nPerm), '_p' ], 'CON_P');












clear 
f2load = 'TR_HPC-TMP_C_4-8Hz'; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = 'TR_HPC-TMP_C'; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 

%%Connectivity across time for every trial
clc
clearvars -except ALLFILT ALLEEG file2load paths f2load

foi = [4:8];
toi = [3001:4750]; % full trial 


for subji = 1:length(ALLFILT)
    disp(['Subji: ' num2str(subji)])

    data    = ALLFILT{subji, 1};
    Ev      = ALLFILT{subji, 2};
    EEG     = ALLEEG{subji}; 

    clear ids2comp
    if ~isempty(data)
        cbH = strsplit(f2load, '_'); 
        contrast = cbH{2}; 
        
        ids2u = Ev(:, 2) == 2; 
        data = data(:, :, ids2u);% only during extinction
        
        ids2comp = extract_ids_EXT(EEG, contrast);
        
        if ~isempty(ids2comp) 
            CON_B = compute_connectivity_EXT(data, toi, [], foi, ids2comp, 'Trials'); % only time in this band
            %size(CON_B, 1)
            CON(subji, :) =  mean(CON_B, 'all');
            %CON{subji, 1} = mean(CON_B, 'all');
            %CON{subji, 2} = Ev(ids2u, :); 
        end
   end
        
end

CON(CON==0) = []; 
%save(['CON_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end))], 'CON');



%%permutations

nPerm = 1000; 
clear CON_P; clear dataPerm; 
clc 

for permi = 1:nPerm
    disp(['Permi: ' num2str(permi)])
    clear CON_BB
    for subji = 1:length(ALLFILT)
        clear dataPerm; 
        
        data    = ALLFILT{subji, 1};
        Ev      = ALLFILT{subji, 2};
        EEG     = ALLEEG{subji}; 
        
        clear ids2comp
        if ~isempty(data)
            for t = 1:size(data, 3)
                % Apply circshift along the time dimension (the 2nd dimension).
                shiftBy =round(size(data, 2)*rand());
                dataPerm(:, :, t) = circshift(data(:, :, t), [0, shiftBy, 0]);
            end
            
            cbH = strsplit(f2load, '_'); 
            contrast = cbH{2}; 
            
            ids2u = Ev(:, 2) == 2; 
            dataPerm = dataPerm(:, :, ids2u);% only during extinction
            
            ids2comp = extract_ids_EXT(EEG, contrast);
            
            if ~isempty(ids2comp) 
                CON_B = compute_connectivity_EXT(dataPerm, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                CON_BB(subji, :) =  mean(CON_B, 'all');
            end
        end
    end
    %size(CON_BB)
    CON_BB(CON_BB==0) = []; 
    %size(CON_BB)
    CON_P(:, permi, :) = CON_BB; 
end



save(['CON_P_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) '_' num2str(nPerm), '_p' ], 'CON_P');











clear 
f2load = 'TR_OFC-TMP_C_4-8Hz'; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = 'TR_OFC-TMP_C'; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 

%%Connectivity across time for every trial
clc
clearvars -except ALLFILT ALLEEG file2load paths f2load

foi = [4:8];
toi = [3001:4750]; % full trial 


for subji = 1:length(ALLFILT)
    disp(['Subji: ' num2str(subji)])

    data    = ALLFILT{subji, 1};
    Ev      = ALLFILT{subji, 2};
    EEG     = ALLEEG{subji}; 

    clear ids2comp
    if ~isempty(data)
        cbH = strsplit(f2load, '_'); 
        contrast = cbH{2}; 
        
        ids2u = Ev(:, 2) == 2; 
        data = data(:, :, ids2u);% only during extinction
        
        ids2comp = extract_ids_EXT(EEG, contrast);
        
        if ~isempty(ids2comp) 
            CON_B = compute_connectivity_EXT(data, toi, [], foi, ids2comp, 'Trials'); % only time in this band
            %size(CON_B, 1)
            CON(subji, :) =  mean(CON_B, 'all');
            %CON{subji, 1} = mean(CON_B, 'all');
            %CON{subji, 2} = Ev(ids2u, :); 
        end
   end
        
end

CON(CON==0) = []; 
%save(['CON_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end))], 'CON');



%%permutations

nPerm = 1000; 
clear CON_P; clear dataPerm; 
clc 

for permi = 1:nPerm
    disp(['Permi: ' num2str(permi)])
    clear CON_BB
    for subji = 1:length(ALLFILT)
        clear dataPerm; 
        
        data    = ALLFILT{subji, 1};
        Ev      = ALLFILT{subji, 2};
        EEG     = ALLEEG{subji}; 
        
        clear ids2comp
        if ~isempty(data)
            for t = 1:size(data, 3)
                % Apply circshift along the time dimension (the 2nd dimension).
                shiftBy =round(size(data, 2)*rand());
                dataPerm(:, :, t) = circshift(data(:, :, t), [0, shiftBy, 0]);
            end
            
            cbH = strsplit(f2load, '_'); 
            contrast = cbH{2}; 
            
            ids2u = Ev(:, 2) == 2; 
            dataPerm = dataPerm(:, :, ids2u);% only during extinction
            
            ids2comp = extract_ids_EXT(EEG, contrast);
            
            if ~isempty(ids2comp) 
                CON_B = compute_connectivity_EXT(dataPerm, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                CON_BB(subji, :) =  mean(CON_B, 'all');
            end
        end
    end
    %size(CON_BB)
    CON_BB(CON_BB==0) = []; 
    %size(CON_BB)
    CON_P(:, permi, :) = CON_BB; 
end



save(['CON_P_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) '_' num2str(nPerm), '_p' ], 'CON_P');










%%




 