

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

%% PLOT ALL ERP FOR ALL FILES IN LOOP 

clear , clc


listF2sav = {   'TR_OFC_C_6_4'
                'TR_FRO_C_6_4'
                'TR_PFC_C_6_4'
                'TR_TMP_C_6_4'
                'TR_OCC_C_6_4'
                'TR_AMY_C_6_4'
                'TR_HPC_C_6_4'
            
            };   

paths = load_paths_EXT; 

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi paths 
    file2load = listF2sav{listi}; 
    load ([paths.results.traces file2load]); 
    
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
                
                figName = [file2load '_s' num2str(subji) '_ch' num2str(chani)]; 
                title(t, figName, 'Interpreter','none'); 
                mkdir(paths.results.tracesPlots)
                exportgraphics(gcf, [paths.results.tracesPlots figName '.png'], 'Resolution', 150); 
                close all; 
    
            end
            
            %c1{subji,:} = tfDCH1; 
            %2{subji,:} = tfDCH2; 
    
        end
    
    
    
    end

disp ('done plotting traces from all subjects')


end





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

%% PLOT ALL SPECTROGRAMS FOR ALL FILES IN LOOP 

clear , clc


listF2sav = {   'TR_OFC_C_6_4'
                'TR_FRO_C_6_4'
                'TR_PFC_C_6_4'
                'TR_TMP_C_6_4'
                'TR_OCC_C_6_4'
                'TR_AMY_C_6_4'
                'TR_HPC_C_6_4'
            
            };   

paths = load_paths_EXT; 

for listi = 1:length(listF2sav)
    disp(['File > ' num2str(listi) '      ' listF2sav{listi}]);
    clearvars -except listF2sav listi paths 
    file2load = listF2sav{listi}; 
    load ([paths.results.traces file2load]); 
    
    
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

        nChans = size(EEG.power, 2);
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
            
            figName = [file2load '_s' num2str(subji) '_ch' num2str(chani)]; 
            title(t, figName, 'Interpreter','none'); 
            exportgraphics(gcf, [paths.results.powerPlots figName '.png'], 'Resolution', 150); 
            close all; 

        end
        
        %c1{subji,:} = tfDCH1; 
        %2{subji,:} = tfDCH2; 

    end



end
disp ('done plotting traces from all subjects')


end



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


































