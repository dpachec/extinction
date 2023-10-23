%% 
%% EXPORT TRACES IN LOOP

clear , clc

c2u = 'C';

listF2sav = {   %{'Amygdala'}; 
%                 {'Hippocampus'}; 
%                 {'orbitofrontal'}; 
%                 {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'parahippocampal' 'entorhinal' };
%                 {'occipital' 'cuneus' 'lingual' 'pericalcarine' 'bankssts'}; 
                {'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
            };   
%n2SAV = {'AMY'; 'HPC'; 'OFC'; 'TMP'; 'OCC'};
n2SAV = {'FRO'};


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';

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

%%EXPORT TRACES IN LOOP

clear , clc

c2u = 'V';

listF2sav = {   %{'Amygdala'}; 
%                 {'Hippocampus'}; 
%                 {'orbitofrontal'}; 
%                 {'inferiortemporal' 'middletemporal' 'superiortemporal' 'transversetemporal' 'fusiform' 'temporalpole' 'parahippocampal' 'entorhinal' };
%                 {'occipital' 'cuneus' 'lingual' 'pericalcarine' 'bankssts'}; 
                {'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' 'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
            };   
%n2SAV = {'AMY'; 'HPC'; 'OFC'; 'TMP'; 'OCC'};
n2SAV = {'FRO'};


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';

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


%%
clear, close all
paths = load_paths_EXT; 

c2u = 'V';

sROI = {'Amygdala'}; 

%sROI = {'Hippocampus'}; 

%sROI = {'orbitofrontal'}; 

%sROI = {'occipital'}; 

%sROI = { 'inferiortemporal' 'middletemporal' 'superiortemporal' 'bankssts' 'fusiform' 'temporalpole' 'occipital' 'lingual' 'parahippocampal' 'cuneus' 'pericalcarine' };


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


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

sROI = char(join(sROI, '_'));
mkdir(paths.results.traces)
filename = [paths.results.traces 'TR_' sROI '_' c2u];
nSub = sum(cell2mat(cellfun(@(x) ~isempty(x), ALLEEG, 'un', 0)));
totalChans = sum(nChans);
save(filename, 'ALLEEG', 'nSub', 'nChans', 'totalChans', '-v7.3');

disp('done');

cd (paths.github)

%% % % % % % % % FIRST LOAD FILE - > START HERE
clear, clc

paths = load_paths_EXT; 
file2load = ['TR_' 'AMY' '_C']; 
%file2load = ['TR_' 'HPC' '_C']; 
%file2load = ['TR_' 'OFC' '_C']; 
%file2load = ['TR_' 'OCC' '_C']; 

load ([paths.results.traces file2load]); 

%%
count = 0; 
for subji = 1:47
    if ~isempty(ALLEEG_AMY{subji}) & ~isempty(ALLEEG_HPC{subji})
        count = count+1; 
    end

end

%% Plot all ERPs

for subji = 1:length(ALLEEG)

    EEG = ALLEEG{subji}; 

    if ~isempty(EEG)
        
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        
        % % %   % % Acquisition
        ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        ids3 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; % Cs+Cs- & Cs-Cs-
        ids4 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; % Cs+Cs-
        ids5 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-

        tfDCH1 = mean(EEG.data(: ,:, ids1), 3, 'omitnan'); 
        tfDCH2 = mean(EEG.data(: ,:, ids2), 3, 'omitnan'); 
        tfDCH3 = mean(EEG.data(: ,:, ids3), 3, 'omitnan'); 
        tfDCH4 = mean(EEG.data(: ,:, ids4), 3, 'omitnan'); 
        tfDCH5 = mean(EEG.data(: ,:, ids5), 3, 'omitnan'); 

        allERPs{subji} = [tfDCH1; tfDCH2; tfDCH3; tfDCH4; tfDCH5];

        nChans = size(EEG.data, 1);
        for chani = 1:nChans

            t = tiledlayout(5, 1); set(gcf, 'Position', [100 100 500 1000])
            nexttile
            plot(tfDCH1(chani, :)); 
            nexttile
            plot(tfDCH2(chani, :)); 
            nexttile
            plot(tfDCH3(chani, :)); 
            nexttile
            plot(tfDCH4(chani, :)); 
            nexttile
            plot(tfDCH5(chani, :)); 
            
            figName = [num2str(subji) '_' num2str(chani)]; 
            title(t, figName, 'Interpreter','none'); 
            mkdir(paths.results.tracesPlots)
            exportgraphics(gcf, [paths.results.tracesPlots figName '.png'], 'Resolution', 150); 
            close all; 

        end
        
        %c1{subji,:} = tfDCH1; 
        %2{subji,:} = tfDCH2; 

    end



end

disp ('done plotting traces')

%% Plot all SPECTROGRAMS

for subji = 1:length(ALLEEG)

    EEG = ALLEEG{subji}; 

    if ~isempty(EEG)
        
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        
        % % %   % % Acquisition
        ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        ids3 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; % Cs+Cs- & Cs-Cs-
        ids4 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; % Cs+Cs-
        ids5 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-

        EEG = extract_power_EXT(EEG, 0.01); 
        EEG = normalize_EXT(EEG);  %across trials

        tfDCH1 = mean(EEG.power(ids1, :, :, :), 1, 'omitnan'); 
        tfDCH2 = mean(EEG.power(ids2, :, :, :), 1, 'omitnan'); 
        tfDCH3 = mean(EEG.power(ids3, :, :, :), 1, 'omitnan'); 
        tfDCH4 = mean(EEG.power(ids4, :, :, :), 1, 'omitnan'); 
        tfDCH5 = mean(EEG.power(ids5, :, :, :), 1, 'omitnan'); 

        nChans = size(EEG.data, 1);
        for chani = 1:nChans

            t = tiledlayout(5, 1); set(gcf, 'Position', [100 100 500 1000])
            nexttile
            imagesc(squeeze(tfDCH1(:, chani, :,:)));
            nexttile
            imagesc(squeeze(tfDCH2(:, chani, :,:))); 
            nexttile
            imagesc(squeeze(tfDCH3(:, chani, :,:))); 
            nexttile
            imagesc(squeeze(tfDCH4(:, chani, :,:))); 
            nexttile
            imagesc(squeeze(tfDCH5(:, chani, :,:))); 
            
            figName = [num2str(subji) '_' num2str(chani)]; 
            title(t, figName, 'Interpreter','none'); 
            exportgraphics(gcf, [paths.results.tracesPlots figName '.png'], 'Resolution', 150); 
            close all; 

        end
        
        %c1{subji,:} = tfDCH1; 
        %2{subji,:} = tfDCH2; 

    end



end

disp ('done plotting spectrograms')

%% Save without plotting

for subji = 1:length(ALLEEG)

    EEG = ALLEEG{subji}; 

    if ~isempty(EEG)
        
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        
        % % %   % % Acquisition
        ids1 = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        ids2 = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');

        % % % % % % Extinction
        ids3 = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids2 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; % Cs+Cs- & Cs-Cs-
        ids4 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; % Cs+Cs-
        ids5 = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; % Cs-Cs-

        tfDCH1 = mean(EEG.data(: ,:, ids1), 3, 'omitnan'); 
        tfDCH2 = mean(EEG.data(: ,:, ids2), 3, 'omitnan'); 
        tfDCH3 = mean(EEG.data(: ,:, ids3), 3, 'omitnan'); 
        tfDCH4 = mean(EEG.data(: ,:, ids4), 3, 'omitnan'); 
        tfDCH5 = mean(EEG.data(: ,:, ids5), 3, 'omitnan'); 

        allERPs{subji, 1} = tfDCH1; 
        allERPs{subji, 2} = tfDCH2; 
        allERPs{subji, 3} = tfDCH3; 
        allERPs{subji, 4} = tfDCH4; 
        allERPs{subji, 5} = tfDCH5; 

        
    end



end


%% plot grand average for CS+ and CS- during Acquisition

c2p1 = allERPs(:,3);
c2pm = cellfun(@(x) mean(x, 1), c2p1, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm1 = cat(1, c2pm{:});
d2p1 = mean(c2pm1); 


c2p2 = allERPs(:, 5);
c2pm = cellfun(@(x) mean(x, 1), c2p2, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm2 = cat(1, c2pm{:});
d2p2 = mean(c2pm2); 



[h p ci ts] = ttest(c2pm1, c2pm2)
h = squeeze(h); t = squeeze(ts.tstat); 
hb = h; hb(h==0) = nan; hb(hb==1) = 0; 

clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
    allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_sum_obs= allSTs(id); 

figure()
plot(d2p1); hold on; 
plot(d2p2); 
plot(hb, LineWidth=4)
set(gca, xlim=[2500 4500])


%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    
    c1B = c2pm1; 
    c2B = c2pm2; 
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :, :);
           c1B(subji, :, :) = c2B(subji, :, :);
           c2B(subji, :, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1B, c2B); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = max(abs(allSTs));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

disp('done')

%%

clear p ratings2u mcsP

ratings2u = max_clust_sum_obs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm





%% 
figure
histogram(max_clust_sum_perm); hold on; 
scatter(max_clust_obs,0, 200, 'filled','r');
set(gca, 'FontSize', 14);
xlabel('T')
exportgraphics(gcf, [paths.results.power 'myP.png'], 'Resolution',150)













































%%


































