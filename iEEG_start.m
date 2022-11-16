%% 
%% 
clear, close all
paths = load_paths; 

c2u = 'C';
sROI = {'Amygdala'}; % case sensitive 

%  sROI = { 'inferiortemporal' 'middletemporal' 'superiortemporal' 'bankssts' 'ctx-lh-fusiform' 'ctx-lh-temporalpole' 
%           'inferiorparietal' 'lateraloccipital' 'lingual' 'parahippocampal' 'cuneus' 'pericalcarine' 'entorhinal'};

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

    chansLab = {EEG.chanlocs.fsLabelsR}';
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
        EEG.event = EEG.event(ids)

        %epoch data and markers
        EEG = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
        EEG = extract_power_EXT(EEG, 0.01); 
        
        %remove amplitude data 
        EEG = rmfield(EEG, 'data');
        EEG = normalize_EXT(EEG);
    
        ALLEEG{subji,:} = EEG; 
    end


end

sROI = char(join(sROI, '_'));
filename = [paths.iEEGRes.power 'allS_' sROI '_' c2u];
save(filename, 'ALLEEG', '-v7.3');

cd (paths.github)

%% delete trigger labels in events

x = [{EEG2.event.type}]; 
ids2rem = strcmp(x, 'trigger');
EEG2.event(ids2rem) = []; 

%% check that markers are ok
%data2check = [EEG.data(1, :); EEG.markers_artifacts(1,:)*1000]; 
data2check = [EEG1.data(1, :)]; 
eegplot(data2check, 'srate', EEG1.srate, 'winlength', 50, 'spacing', 1000, 'events', EEG1.event);




%% PLOT grand average for each condition
paths = load_paths; 
file2load = ['allS_' 'orbitofrontal' '_U']; 
%load ([paths.iEEGRes.power file2load]); 
clearvars -except ALLEEG paths file2load

c2u = file2load(end);


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    
    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        if ndims(EEG.power) == 3
            tmph(:, 1, :, :)  = EEG.power; 
            EEG.power = tmph; 
        end

        ids1 = strcmp(Ev2(:, 10), c2u) & ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) & strcmp(Ev2(:, 2), '1');
        tfDCH1 = mean(EEG.power(ids1, :, : ,:), 'omitnan'); 
        tfDTF1 = squeeze(mean(tfDCH1, 2, 'omitnan'));
        ids2 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '3')  & strcmp(Ev2(:, 2), '1');
        tfDCH2 = mean(EEG.power(ids2, :, : ,:), 'omitnan'); 
        tfDTF2 = squeeze(mean(tfDCH2, 2, 'omitnan'));

        c1(subji, :, :) = tfDTF1; 
        c2(subji, :, :) = tfDTF2; 
        
    end

end


cd (paths.github)

%%

sub2exc = []

c1B = c1; c2B = c2; 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));


[h p ci ts] = ttest(c1, c2); 
h = squeeze(h); t = squeeze(ts.tstat);

%h = zeros(54, 300);
%h(clustinfo.PixelIdxList{45}) = 1; 



%times = -1:.01:1.99; 
times = -3:.01:3.99; 
freqs = 1:54;
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 800])
nexttile
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; colorbar
nexttile
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; colorbar
nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; colorbar
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 1); set(gca, 'clim', [-3 4])
colormap(brewermap([],'Spectral'))
set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'});

%exportgraphics(gcf, [paths.iEEGRes.power file2load '.png'], 'Resolution',150)
exportgraphics(gcf, [paths.iEEGRes.power  'myP.png'], 'Resolution',150)



%% stats 
clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 


%% permutations 

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = c1(:,:,101:300); 
    c2B = c2(:,:,101:300); 
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :, :);
           c1B(subji, :, :) = c2B(subji, :, :);
           c2B(subji, :, :) = tmp; 
        end
    end
    
    CP1	= squeeze(mean(c1B, 'omitnan'));
    CP2	= squeeze(mean(c2B, 'omitnan'));
    [hPerm p ci tsPerm] = ttest(CP1, CP2); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
    else
        allSTs = 0; 
    end
    max_clust_sum_perm(permi,:) = abs(allSTs(id)); 

end

%%
clear p mcsR mcsP

mcsR = max_clust_obs; 
mcsP = max_clust_sum_perm;

allAb = mcsP(abs(mcsP) > abs(mcsR))';
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;




%% 
figure
histogram(max_clust_sum_perm); hold on; 
scatter(abs(max_clust_obs),0, 'filled','r');
set(gca, 'FontSize', 14);
xlabel('T')
exportgraphics(gcf, [paths.iEEGRes.power 'myP.png'], 'Resolution',150)



%% plot 2 bars with mean theta 







%% count channels 

clear nChans
for subji = 1:length(ALLEEG)

    EEG = ALLEEG{subji};
    if~isempty(EEG)
        nChans(subji,:) = length(EEG.chanlocs);
    end


end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])







%% Searchlight analysis 






























%% 
d2p = squeeze(EEG.data(1, :, 1)); 
figure;plot(d2p)



%%
data2check = [EEG.data(1, :); EEG.marker_artifacts(1,:)*1000]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000);










%% plot example trial in one subject (WITH MORE THAN 1 electrode)

%EEG = ALLEEG{1}; 
tr = 20; 
ch = 2; 

% figure
% d2p	= squeeze(EEG.dsPower(tr, ch, : ,:));
% myCmap = colormap(brewermap([],'YlOrRd'));
% colormap(myCmap)
% contourf(1:300, 1:54, d2p, 40, 'linecolor', 'none'); colorbar

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
d2p	= squeeze(EEG.power(tr, ch, : ,:));
contourf(1:70, 1:54, d2p, 40, 'linecolor', 'none'); colorbar


%% plot example trial in one subject (ONLY 1 electrode)

%EEG = ALLEEG{1}; 
tr = 12; 

figure
d2p	= squeeze(EEG.dsPower(tr, : ,:));
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(1:300, 1:54, d2p, 40, 'linecolor', 'none'); colorbar

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
d2p	= squeeze(EEG.power(tr,: ,:));
contourf(1:3000, 1:54, d2p, 40, 'linecolor', 'none'); colorbar


%% plot two different conditions (average in all amygdala channels and across trials)


Ev = [{EEG.event.type}]';
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
Ev2 = cat(1, Ev1{:})
ids = strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1') % CS+CS+ during acquisition
d2p1	= squeeze(mean(mean(EEG.power(ids, :, : ,:))));
ids = strcmp(Ev2(:, 6), '2')  & strcmp(Ev2(:, 2), '1') % CS+CS- during acquisition
d2p2	= squeeze(mean(mean(EEG.power(ids, :, : ,:))));

myCmap = colormap(brewermap([],'YlOrRd')); colormap(myCmap)

tiledlayout(2, 1,'TileSpacing','compact');
nexttile
contourf(1:701, 1:54, d2p1, 40, 'linecolor', 'none'); hold on; %colorbar
nexttile
contourf(1:701, 1:54, d2p2, 40, 'linecolor', 'none'); hold on; %colorbar












%% power analysis all channels / SETTINGS FOR PLOTTING IN BRAIN PLOTS 
clear, close all
paths = load_paths; 

c2u = 'C';
sROI = 'allChannels'; % don't change > for specific ROIs check next block 

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


for subji = 1:length(allsubs)
    subji
    sub = allsubs{subji}; 
    cd([paths.iEEG])
    load ([sub '_iEEG.mat']);

    
    Ev = [{EEG.event.type}]'; 
    Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
    clen = cellfun(@length, Ev1); 
    EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
    Ev2 = cat(1, Ev1{:});
   
    Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    ids = strcmp(Ev2(:, 10), c2u); 
    EEG.event = EEG.event(ids)

    %epoch data and markers
    EEG = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
    EEG = extract_power_EXT(EEG, 0.1); 
    EEG.power = EEG.power(:, :, :, 31:40);
    
    EEG = rmfield(EEG, 'data');
    EEG = normalize_EXT(EEG);

    ALLEEG{subji,:} = EEG; 


end




sROI = char(join(sROI, '_'));
filename = [paths.iEEGRes.power 'allS_' sROI '_' c2u];
save(filename, "ALLEEG");

cd (paths.github)


%% 