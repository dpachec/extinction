%% 
%% 
clear, close all
paths = load_paths; 

c2u = 'U';

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17', 'p_sub18'}';


for subji = 1:length(allsubs)
    subji
    sub = allsubs{subji}; 
    cd([ paths.fiEEG])
    load ([sub '_iEEG.mat']);

    %select amygdala electrodes
    chansLab = {EEG.chanlocs.fsLabel}';
    selChans = contains(chansLab, 'Amygdala');

    EEG = artifact_detection(EEG, 4, 200, 100);

    if find(selChans)
        EEG.chanlocs = EEG.chanlocs(selChans);
        EEG.data = EEG.data(selChans, :);
        EEG.marker_artifacts = EEG.marker_artifacts(selChans,:);
    
        %epoch data
        Ev = [{EEG.event.type}]'; 
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        clen = cellfun(@length, Ev1); 
        EEG.event = EEG.event(clen==10); Ev1 = Ev1(clen==10);
        Ev2 = cat(1, Ev1{:});
       
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %paris subjects havespace in the last character of the event WHY??
        ids = strcmp(Ev2(:, 10), 'U'); 
        EEG.event = EEG.event(ids)
        EEG = pop_epoch( EEG, {}, [-3 4], 'newname', 'verbose', 'epochinfo', 'yes');
        % remove trials with artifacts
        EEG = extract_power_EXT(EEG, 0.01); 
        %remove edge artifacts
        if ndims(EEG.power) == 4
            EEG.power = EEG.power(:, :, :, 201:500);
        else
            EEG.power = EEG.power(:, :, 201:500);
        end
        %remove amplitude data 
        EEG = rmfield(EEG, 'data');
        EEG = rmfield(EEG, 'marker_artifacts');
        EEG = normalize_EXT(EEG);
    
        ALLEEG{subji,:} = EEG; 
    end


end

filename = [paths.iEEGRes.power 'allS_' c2u];
save(filename, "ALLEEG");



%% check that markers are ok
data2check = [EEG.data(1, :); EEG.marker_artifacts(1,:)*1000]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000);







%% plot example trial in one subject

EEG = ALLEEG{1}; 
tr = 1; 
ch = 1; 
d2p	= squeeze(EEG.power(tr, ch, : ,:));

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
contourf(1:701, 1:54, d2p, 40, 'linecolor', 'none'); hold on; %colorbar



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


%% PLOT grand average for each condition
paths = load_paths; 

if ~exist('ALLEEG') load ([paths.iEEGRes.power 'allS_U']); end
clearvars -except ALLEEG


c2u = 'U'

count = 1; 
for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};

    if ~isempty(EEG)
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??

        if ndims(EEG.power) == 4
            ids1 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1'); % CS+CS+ during acquisition
            d2p1	= squeeze(mean(mean(EEG.power(ids1, :, : ,:), 'omitnan'), 'omitnan'));
            ids2 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '3')  & strcmp(Ev2(:, 2), '1');  % CS+CS- during acquisition
            d2p2	= squeeze(mean(mean(EEG.power(ids2, :, : ,:), 'omitnan'), 'omitnan'));
        else
            ids1 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1');  % CS+CS+ during acquisition
            d2p1	= squeeze(mean(EEG.power(ids1, :, : ), 'omitnan'));
            ids2 = strcmp(Ev2(:, 10), c2u) & strcmp(Ev2(:, 6), '3')  & strcmp(Ev2(:, 2), '1');  % CS+CS- during acquisition
            d2p2	= squeeze(mean(EEG.power(ids2, :, : ), 'omitnan'));

        end
    
        c1(count, :, :) = d2p1; 
        c2(count, :, :) = d2p2; 
        count = count+1; 
    end

end

%%

d2p1	= squeeze(mean(c1, 'omitnan'));
d2p2	= squeeze(mean(c2, 'omitnan'));

figure
%myCmap = colormap(brewermap([],'YlOrRd'));
%colormap(myCmap)
contourf(1:300, 1:54, d2p1, 40, 'linecolor', 'none'); hold on; colorbar
figure
contourf(1:300, 1:54, d2p2, 40, 'linecolor', 'none'); hold on; colorbar

figure
contourf(1:300, 1:54, d2p1-d2p2, 40, 'linecolor', 'none'); hold on; colorbar


%% stats 

[h p ci ts] = ttest(c1, c2); 
h = squeeze(h); t = squeeze(ts.tstat);
figure
contourf(1:300, 1:54, h, 40, 'linecolor', 'none'); hold on; colorbar
















































%% 
d2p = squeeze(EEG.data(1, :, 1)); 
figure;plot(d2p)



%%
data2check = [EEG.data(1, :); EEG.marker_artifacts(1,:)*1000]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000);






















%% 