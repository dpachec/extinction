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
                    amyDF = eegfilt (amyD,1;;, f2u(1), f2u(end));
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
f2load = 'TR_AMY-TMP_C_4-8Hz'; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = 'TR_HPC-TMP_C'; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 

%% Connectivity across time for every trial
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




%%
clear 

ctr2u = ['AMY-TMP']; 

f2load = ['TR_' ctr2u '_C_4-8Hz']; 
paths = load_paths_EXT; 
load ([paths.results.rsa f2load]); 

f2load1 = ['TR_' ctr2u '_C']; % hack because electrode data is missing in ALLFIT
load ([paths.results.traces f2load1]); 


%% Connectivity across time for every trial and CONDITION
clc
clearvars -except ALLFILT ALLEEG file2load paths f2load ctr2u

foi = [4:8];
toi = [3001:4750]; % full trial 


for condi = 1:3
    clear CON
    for subji = 1:length(ALLFILT)
        disp(['Subji: ' num2str(subji)])
    
        data    = ALLFILT{subji, 1};
        Ev      = ALLFILT{subji, 2};
        EEG     = ALLEEG{subji}; 
    
        clear ids2comp 
        if ~isempty(data)
            cbH = strsplit(f2load, '_'); 
            contrast = cbH{2}; 
            
            
            ids2u = Ev(:, 2) == 2 &  Ev(:, 6) ==  condi; 
            dataC = data(:, :, ids2u);% only during extinction
            
            ids2comp = extract_ids_EXT(EEG, contrast);
            
            if ~isempty(ids2comp) & ~isempty(dataC) 
                CON_B = compute_connectivity_EXT(dataC, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                %size(CON_B, 1)
                CON(subji, :) =  mean(CON_B, 'all');
                %CON{subji, 1} = mean(CON_B, 'all');
                %CON{subji, 2} = Ev(ids2u, :); 
            end
        end
            
   end
   CON(CON==0) = nan; 
   CON3{condi} = CON; 
end


save(['CON3_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) ], 'CON3');

%% addmissinsubjectjustincase
%CON3{2}(end+1) = nan; 
%CON3{3}(end+1) = nan; 

% allL = cellfun(@length, CON3)
% [minL id] = min(allL); maxL = max(allL); 
% if minL < maxL
%     CON3{id}(end+ (maxL-minL)) = nan; 
% end

%%remove subjects with nans
sub2exc = []; 
for coni = 1:3
    CON = CON3{coni}; 
    ids = find(isnan(CON)); 
    sub2exc = [sub2exc ; ids]; 
end

%%ANOVA REPEATED MEASURES
clc 
data = [CON3{1}, CON3{2}, CON3{3}];
data(sub2exc, :) = []; 
nSubj = size(data, 1); 
d4anova = data(:);
d4anova = d4anova(:); 
d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];

[p f] = RMAOV1(d4anova);

boxplot(data)
set(gca, 'FontSize', 16)
%title(['F=' num2str(f) '        ' 'p=' num2str(p)])
exportgraphics(gcf, 'myP.png', 'Resolution',150)

%disp(['F(2, ' num2str ])






%% IN LOOP FOR ALL PAIR OF ROIS
clear 

%roiPAIRS = {['AMY-TMP'] ['AMY-HPC'] ['AMY-PFC'] ['TMP-HPC'] ['TMP-OFC'] ['TMP-PFC'] }; 
%toiALLRS = {[3001:4750] [3501:4750] [4001:4650] [3001:4450] [3150:4000] [3001:3800] }

roiPAIRS = { ['HPC-TMP'] ['OFC-TMP'] ['PFC-TMP'] }; 
toiALLRS = {  [3001:4450] [3150:4000] [3001:3800] }


for contri = 1:6

    clearvars -except roiPAIRS toiALLRS contri

    ctr2u = roiPAIRS{contri}; 
    toi   = toiALLRS{contri};
    
    f2load = ['TR_' ctr2u '_C_4-8Hz']; 
    paths = load_paths_EXT; 
    load ([paths.results.rsa f2load]); 

    f2load1 = ['TR_' ctr2u '_C']; % hack because electrode data is missing in ALLFIT
    load ([paths.results.traces f2load1]); 
    
    %%Connectivity across time for every trial and CONDITION
    clc
    clearvars -except ALLFILT ALLEEG file2load paths f2load ctr2u toi roiPAIRS toiALLRS
    
    foi = [4:8];
    
    
    for condi = 1:3
        clear CON
        for subji = 1:length(ALLFILT)
            disp(['Subji: ' num2str(subji)])
        
            data    = ALLFILT{subji, 1};
            Ev      = ALLFILT{subji, 2};
            EEG     = ALLEEG{subji}; 
        
            clear ids2comp 
            if ~isempty(data)
                cbH = strsplit(f2load, '_'); 
                contrast = cbH{2}; 
                
                
                ids2u = Ev(:, 2) == 2 &  Ev(:, 6) ==  condi; 
                dataC = data(:, :, ids2u);% only during extinction
                
                ids2comp = extract_ids_EXT(EEG, contrast);
                
                if ~isempty(ids2comp) & ~isempty(dataC) 
                    CON_B = compute_connectivity_EXT(dataC, toi, [], foi, ids2comp, 'Trials'); % only time in this band
                    %size(CON_B, 1)
                    CON(subji, :) =  mean(CON_B, 'all');
                    %CON{subji, 1} = mean(CON_B, 'all');
                    %CON{subji, 2} = Ev(ids2u, :); 
                end
            end
                
       end
       CON(CON==0) = nan; 
       CON3{condi} = CON; 
    end
    
    
    save(['CON3_' contrast '_' num2str(foi(1)) '-' num2str(foi(end)) '_' num2str(toi(1)) '-' num2str(toi(end)) ], 'CON3');
    
    
    

end







%% addmissinsubjectjustincase
clearvars -except CON3

% allL = cellfun(@length, CON3)
% [minL id] = min(allL); maxL = max(allL); 
% if minL < maxL
%     CON3{id}(end+ (maxL-minL)) = nan; 
% end

%%remove subjects with nans
sub2exc = []; 
for coni = 1:3
    CON = CON3{coni}; 
    ids = find(isnan(CON)); 
    sub2exc = [sub2exc ; ids]; 
end

%%ANOVA REPEATED MEASURES

data = [CON3{1}, CON3{2}, CON3{3}];
data(sub2exc, :) = []; 
nSubj = size(data, 1); 
d4anova = data(:);
d4anova = d4anova(:); 
d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];

[p f] = RMAOV1(d4anova);

boxplot(data)
set(gca, 'FontSize', 16)
%title(['F=' num2str(f) '        ' 'p=' num2str(p)])
exportgraphics(gcf, 'myP.png', 'Resolution',150)

%disp(['F(2, ' num2str ])


%%