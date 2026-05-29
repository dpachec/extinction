%% % % % % % % % FIRST LOAD FILE - > START HERE
%% count number of trials for each condition > POWER
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'TMP' '_C']; 

load ([paths.results.power file2load]); 

%% count number of subjects with only 1 electode in AMY
clearvars -except ALLEEG paths file2load
clc 
allPats = []; 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        
        chans = [{EEG.chanlocs.fsLabel}]'; 
        if length(chans)== 1
            allPats = [allPats subji]; 
        end
    end
end


%% count number of electrodes in Guangzhou and Paris
clearvars -except ALLEEG paths file2load
clc 
all_chansChina = []; 
all_chansParis = []; 
allPatsChina = []; 
allPatsParis = []; 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        if subji < 31
            chansChina = [{EEG.chanlocs.fsLabel}]'; 
            all_chansChina = [all_chansChina ; chansChina(:)]
            allPatsChina = [allPatsChina subji]; 
        else
            chansParis = [{EEG.chanlocs.fsLabel}]'; 
            all_chansParis = [all_chansParis ; chansParis(:)]
            allPatsParis = [allPatsParis subji]; 
        end
    end
end

%% barplot
% Plot with a pie chart
data = [length(allPatsChina) length(allPatsParis)]; 
bar(1, data, 'stacked');
xlabel('Entries');
ylabel('Counts');
%title('Counts of Each Entry in the Cell Array');

exportgraphics(gca, 'myP.png', 'Resolution',300)


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

nChans2 = nChans(any(nChans,2),:);
disp (['Elec per subject in amy: ' num2str(mean(nChans2)) ' ± ' num2str(std(nChans2))])



%% Compute number of trials for each condition

clearvars -except ALLEEG paths  totalChans nChans nSub nL file2load


for subji = 1:length(ALLEEG)
    
    EEG = ALLEEG{subji};
    
    if ~isempty(EEG)
        
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); %sub33 has some space in the last character of the event WHY??


        % % % % CS++
        ids1 = (strcmp(Ev2(:, 2), '1')  | strcmp(Ev2(:, 2), '2') ) &  strcmp(Ev2(:, 6), '1')  ;

        % % % % CS+-
        ids2 = (strcmp(Ev2(:, 2), '1')  | strcmp(Ev2(:, 2), '2') ) &  strcmp(Ev2(:, 6), '2')  ;

        % % % % CS--
        ids3 = (strcmp(Ev2(:, 2), '1')  | strcmp(Ev2(:, 2), '2') ) &  strcmp(Ev2(:, 6), '3')  ;


        
        
        nTrials(subji, :) = get_number_of_trials_POW_EXT(EEG, ids1, ids2, ids3); 



                
    end

end


nTrials2 = nTrials(any(nTrials,2),:);
disp (['Trials per subject cs++: ' num2str(mean(nTrials2(:, 1))) ' ± ' num2str(std(nTrials2(:, 1)))])
disp (['Tirals per subject cs+-: ' num2str(mean(nTrials2(:, 2))) ' ± ' num2str(std(nTrials2(:, 2)))])
%disp (['Elec per subject cs--: ' num2str(mean(nTrials2(:, 3))) ' ± ' num2str(std(nTrials2(:, 3)))])





%% count number of trials for each condition > Stability

clear, clc

paths = load_paths_EXT; 
f2sav = 'RSA_AMY_C_1-44_1_0_500-50_1_T_SICSPA-SICSMA';
load ([ paths.results.rsa f2sav '.mat']);
cfg = getParams_EXT(f2sav);
sub2exc = [37]; %
nTrials = compute_trial_number_EXT(ids); 

nTrials2 = nTrials(any(nTrials,2),:);
disp (['Trials per subject cs+: ' num2str(mean(nTrials2(:, 1))) ' ± ' num2str(std(nTrials2(:, 1)))])
disp (['Tirals per subject cs-: ' num2str(mean(nTrials2(:, 2))) ' ± ' num2str(std(nTrials2(:, 2)))])





%% count channels 

clear nChans
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if~isempty(EEG)
        length(EEG.chanlocs)
        nChans(subji,:) = length(EEG.chanlocs);
    end
end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])



%% COUNT INTERDISTANCES IN REAL DATA

clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 

load ([paths.results.power file2load]); 


%% count interdistances at each electrode
clearvars -except ALLEEG paths file2load
clc 
tic
all_distsGuangzhou = []; 
for subji = 1:30
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 
        distsLeft = pdist(elecLeft, 'euclidean');
        distsRight = pdist(elecRight, 'euclidean');
        LRD = [distsLeft distsRight]; 
        all_distsGuangzhou = [all_distsGuangzhou LRD];
    end
end
all_distsParis = []; 
for subji = 31:length(ALLEEG)
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 
        distsLeft = pdist(elecLeft, 'euclidean');
        distsRight = pdist(elecRight, 'euclidean');
        LRD = [distsLeft distsRight]; 
        all_distsParis = [all_distsParis LRD];
    end
end

toc

%% plot 2 bar
clc
figure(2); set(gcf,'Position', [0 0 500 650]); 
scatter(1, all_distsGuangzhou, 600, 'k.'); hold on; 
scatter(2, all_distsParis, 600, 'k.'); 
%set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [0 10] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

set(gca, 'LineWidth', 3);

exportgraphics(gca, 'myP.png', 'Resolution',300)



%% COUNT INTERDISTANCES FOR ALL ELECTRODES TOGETHER SEPARATELY FOR PARIS AND GUANGZHOU

clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'TMP' '_C']; 

load ([paths.results.power file2load]); 


%% count interdistances at ALL electrodes together in one hemisphere
clearvars -except ALLEEG paths file2load
clc 
tic
all_EL_Guangzhou = []; 
all_ER_Guangzhou = []; 
for subji = 1:30
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 
        
        all_EL_Guangzhou = cat(1, all_EL_Guangzhou, elecLeft);
        all_ER_Guangzhou = cat(1, all_ER_Guangzhou, elecRight);
    end
end
all_EL_Paris = []; 
all_ER_Paris = [];
for subji = 31:length(ALLEEG)
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 

        all_EL_Paris = cat(1, all_EL_Paris, elecLeft);
        all_ER_Paris = cat(1, all_ER_Paris, elecRight);

    end
end


distsGuangzhouLeft = pdist(all_EL_Guangzhou, 'euclidean');
distsGuangzhouRight = pdist(all_ER_Guangzhou, 'euclidean');
distsParisLeft = pdist(all_EL_Paris, 'euclidean');
distsParisRight = pdist(all_ER_Paris, 'euclidean');

toc

%% 
allData = [distsGuangzhouLeft'; distsGuangzhouRight'; distsParisLeft'; distsParisRight'];

% Create corresponding group labels
groupLabels = [repmat(1, length(distsGuangzhouLeft), 1); 
               repmat(2, length(distsGuangzhouRight), 1); 
               repmat(3, length(distsParisLeft), 1); 
               repmat(4, length(distsParisRight), 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(allData, groupLabels);

% Display the ANOVA table
%disp(tbl);

% Post-hoc multiple comparison test (Tukey's test)
multcompare(stats);

%% Histograms

dG = [distsGuangzhouLeft'; distsGuangzhouRight'];
dP = [distsParisLeft'; distsParisRight'];
[h p ci ts] = ttest2(dG, dP)

mean(dG)
mean(dP)

figure
histogram(dG); hold on; 
histogram(dP);
set(gca, 'FontSize', 15)
myD = cohens_d(dG, dP);

exportgraphics(gca, 'myP.png', 'Resolution',300)


%% Two-Factor ANOVA with Unequal Group Sizes

factorA = [zeros(1, length(distsGuangzhouLeft) + length(distsGuangzhouRight)) ones(1, length(distsParisLeft) + length(distsParisRight))]; 
factorB = [zeros(1, length(distsGuangzhouLeft)) ones(1, length(distsGuangzhouRight)) zeros(1, length(distsParisLeft)) ones(1, length(distsParisRight))]; 
allData = [distsGuangzhouLeft'; distsGuangzhouRight'; distsParisLeft'; distsParisRight'];

% Run Two-Way ANOVA (Type III SS)
[p, tbl, stats] = anovan(allData, {factorA, factorB}, 'model', 'interaction', ...
    'varnames', {'FactorA', 'FactorB'});

% Display results
disp(tbl);

% Post-hoc multiple comparisons
multcompare(stats, 'Dimension', 1);  % Compare Factor A
%multcompare(stats, 'Dimension', 2);  % Compare Factor B


%% violin Plot
clc
data = {distsGuangzhouLeft, distsGuangzhouRight, distsParisLeft, distsParisRight}; 

for violi = 1:4

    violinplot(violi, data{violi}); hold on; 
    plot([violi-0.2 violi+0.2], [median(data{violi}) median(data{violi})], 'k', 'LineWidth', 2);

end
xlim([0.5 numel(data) + 0.5]);
xticks(1:numel(data));
set(gca, 'xticklabels', {'Guangzhou Left', 'Guangzhou Right', 'Paris Left', 'Paris Right'});
ax = gca;
set(gca, 'FontSize', 14)
ylabel('Inter-electrode distances');
title('AMY')
ax = gca;
hold off;

exportgraphics(gca, 'myP.png', 'Resolution',300)




%% COUNT INTERDISTANCES IN EACH SUBJECT and Hemisphere

clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'AMY' '_C']; 

load ([paths.results.power file2load]); 

%% count interdistances IN EACH SUBJECT and Hemisphere
clearvars -except ALLEEG paths file2load
clc 
tic

for subji = 1:30
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 
        
        all_EL_Guangzhou{subji, :} = elecLeft;
        all_ER_Guangzhou{subji, :} = elecRight;
    end
end

for subji = 31:length(ALLEEG)
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 

        all_EL_Paris{subji, :} = elecLeft;
        all_ER_Paris{subji, :} = elecRight;

    end
end


if length(all_EL_Guangzhou) < 30
    all_EL_Guangzhou{30} = []; 
end
if length(all_ER_Guangzhou) < 30
    all_ER_Guangzhou{30} = []; 
end
if length(all_EL_Paris) < 50
    all_EL_Paris{50} = []; 
end
if length(all_ER_Paris) < 50
    all_ER_Paris{50} = []; 
end


%%
clear elecH
for subji = 1:30
    elecH = all_EL_Guangzhou{subji}; 
    if ~isempty(elecH)
        distS_LG(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
distS_LG(distS_LG == 0) = []; 
distS_LG(isnan(distS_LG)) = []; 

clear elecH
for subji = 1:30
    elecH = all_ER_Guangzhou{subji}; 
    if ~isempty(elecH)
        distS_RG(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
distS_RG(distS_RG == 0) = []; 
distS_RG(isnan(distS_RG)) = []; 

clear elecH
for subji = 31:50
    elecH = all_EL_Paris{subji}; 
    if ~isempty(elecH)
        distS_LP(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
distS_LP(distS_LP == 0) = []; 
distS_LP(isnan(distS_LP)) = []; 

clear elecH
for subji = 31:50
    elecH = all_ER_Paris{subji}; 
    if ~isempty(elecH)
        distS_RP(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
distS_RP(distS_RP == 0) = []; 
distS_RP(isnan(distS_RP)) = []; 




%% violin Plot
clc
data = {distS_LG, distS_RG, distS_LP, distS_RP}; 

for violi = 1:4

    violinplot(violi, data{violi}); hold on; 
    plot([violi-0.2 violi+0.2], [median(data{violi}) median(data{violi})], 'k', 'LineWidth', 2);

end
xlim([0.5 numel(data) + 0.5]);
xticks(1:numel(data));
set(gca, 'xticklabels', {'Guangzhou Left', 'Guangzhou Right', 'Paris Left', 'Paris Right'});
ax = gca;
set(gca, 'FontSize', 14)
ylabel('Inter-electrode distances');
title('AMY')
ax = gca;
hold off;

exportgraphics(gca, 'myP.png', 'Resolution',300)


%% 
allData = [distS_LG; distS_RG; distS_LP; distS_RP];

% Create corresponding group labels
groupLabels = [repmat(1, length(distS_LG), 1); 
               repmat(2, length(distS_RG), 1); 
               repmat(3, length(distS_LP), 1); 
               repmat(4, length(distS_RP), 1)];

% Perform one-way ANOVA
[p, tbl, stats] = anova1(allData, groupLabels);

% Display the ANOVA table
%disp(tbl);

% Post-hoc multiple comparison test (Tukey's test)
multcompare(stats);





%% NO HEMISPHERE DISTINCTION
clear distS_LG distS_RG distS_LP distS_RP
clear elecH
for subji = 1:30
    elecH = all_EL_Guangzhou{subji}; 
    if ~isempty(elecH)
        distS_LG(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
%distS_LG(distS_LG == 0) = []; 
%distS_LG(isnan(distS_LG)) = []; 

clear elecH
for subji = 1:30
    elecH = all_ER_Guangzhou{subji}; 
    if ~isempty(elecH)
        distS_RG(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
%distS_RG(distS_RG == 0) = []; 
%distS_RG(isnan(distS_RG)) = []; 

clear elecH
for subji = 31:50
    elecH = all_EL_Paris{subji}; 
    if ~isempty(elecH)
        distS_LP(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
%distS_LP(distS_LP == 0) = []; 
%distS_LP(isnan(distS_LP)) = []; 

clear elecH
for subji = 31:50
    elecH = all_ER_Paris{subji}; 
    if ~isempty(elecH)
        distS_RP(subji, :) = mean(pdist(elecH, 'euclidean'));
    end
end
%distS_RP(distS_RP == 0) = []; 
%distS_RP(isnan(distS_RP)) = []; 

distS_LG(end+1:30) = nan; 
distS_RG(end+1:30) = nan; 
distS_LP(end+1:50) = nan; 
distS_RP(end+1:50) = nan; 

distS_LG(distS_LG==0) = nan; 
distS_RG(distS_RG==0) = nan; 
distS_LP(distS_LP==0) = nan; 
distS_RP(distS_RP==0) = nan; 


% % % average in each subject first
distS_G = mean([distS_LG distS_RG], 2, 'omitnan');
% % % average in each subject first
distS_P = mean([distS_LP distS_RP], 2, 'omitnan');


distS_G(distS_G == 0) = []; 
distS_G(isnan(distS_G)) = []; 

distS_P(distS_P == 0) = []; 
distS_P(isnan(distS_P)) = []; 

%% 

[h p ci ts] = ttest2(distS_G, distS_P)

%% Histograms

dataG = distS_G; 
dataP = distS_P; 
[h p ci ts] = ttest2(dataG, dataP)

figure
histogram(dataG); hold on; 
histogram(dataP);
set(gca, 'FontSize', 15)
myD = cohens_d(dataG, dataP)

exportgraphics(gca, 'myP.png', 'Resolution',300)

%% 
bar (1, mean(dataG)); hold on
bar(2, mean(dataP))

%% Two-Factor ANOVA with Unequal Group Sizes

factorA = [zeros(1, length(distS_LG) + length(distS_RG)) ones(1, length(distS_LP) + length(distS_RP))]; 
factorB = [zeros(1, length(distS_LG)) ones(1, length(distS_RG)) zeros(1, length(distS_LP)) ones(1, length(distS_RP))]; 
allData = [distS_LG; distS_RG; distS_LP; distS_RP];

% Run Two-Way ANOVA (Type III SS)
[p, tbl, stats] = anovan(allData, {factorA, factorB}, 'model', 'interaction', ...
    'varnames', {'FactorA', 'FactorB'});

% Display results
disp(tbl);

% Post-hoc multiple comparisons
multcompare(stats, 'Dimension', 1);  % Compare Factor A
multcompare(stats, 'Dimension', 2);  % Compare Factor B


%% COMPUTE DELAUNAY DISTANCE FOR EVERY SUBJECT
clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'TMP' '_C']; 

load ([paths.results.power file2load]); 

%% count interdistances IN EACH SUBJECT and Hemisphere
clearvars -except ALLEEG paths file2load
clc 
tic

for subji = 1:30
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 
        
        all_EL_Guangzhou{subji, :} = elecLeft;
        all_ER_Guangzhou{subji, :} = elecRight;
    end
end

for subji = 31:length(ALLEEG)
    EEG = ALLEEG{subji};
    distsLeft = []; 
    distsRight = []; 
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        idsLeft =  find(contains(chans, 'Left') | contains(chans, 'lh')); 
        idsRight =  find(contains(chans, 'Right') | contains(chans, 'rh')); 

        elecLeft = cat(1, EEG.chanlocs(idsLeft).mniCoordR); 
        elecRight= cat(1, EEG.chanlocs(idsRight).mniCoordR); 

        all_EL_Paris{subji, :} = elecLeft;
        all_ER_Paris{subji, :} = elecRight;

    end
end


if length(all_EL_Guangzhou) < 30
    all_EL_Guangzhou{30} = []; 
end
if length(all_ER_Guangzhou) < 30
    all_ER_Guangzhou{30} = []; 
end
if length(all_EL_Paris) < 50
    all_EL_Paris{50} = []; 
end
if length(all_ER_Paris) < 50
    all_ER_Paris{50} = []; 
end


%%
clc
for subji = 1:30
    elech2 = all_EL_Guangzhou{subji}; 
    if ~isempty(elech2) & size(elech2, 1) > 3
        electrode_coords = elech2;
        [K, volume] = convhull(electrode_coords(:,1), electrode_coords(:,2), electrode_coords(:,3));
        num_triangles = size(K, 1);
        triangle_areas = zeros(num_triangles, 1);
        
        for i = 1:num_triangles
            % Extract three vertices of the triangle
            p1 = electrode_coords(K(i,1), :);
            p2 = electrode_coords(K(i,2), :);
            p3 = electrode_coords(K(i,3), :);
        
            % Compute the area using the cross product formula in 3D
            v1 = p2 - p1;
            v2 = p3 - p1;
            triangle_areas(i) = 0.5 * norm(cross(v1, v2));
        end
        
        % Compute the mean triangle area
        mean_triangle_area = mean(triangle_areas);
        
        mTA_LG(subji, :) = mean_triangle_area;
        % Display results
        fprintf('Mean Triangle Area: %.4f\n', mean_triangle_area);
    end

end

for subji = 1:30
    elech2 = all_ER_Guangzhou{subji}; 
    if ~isempty(elech2) & size(elech2, 1) > 3
        electrode_coords = elech2;
        [K, volume] = convhull(electrode_coords(:,1), electrode_coords(:,2), electrode_coords(:,3));
        num_triangles = size(K, 1);
        triangle_areas = zeros(num_triangles, 1);
        
        for i = 1:num_triangles
            % Extract three vertices of the triangle
            p1 = electrode_coords(K(i,1), :);
            p2 = electrode_coords(K(i,2), :);
            p3 = electrode_coords(K(i,3), :);
        
            % Compute the area using the cross product formula in 3D
            v1 = p2 - p1;
            v2 = p3 - p1;
            triangle_areas(i) = 0.5 * norm(cross(v1, v2));
        end
        
        % Compute the mean triangle area
        mean_triangle_area = mean(triangle_areas);
        
        mTA_RG(subji, :) = mean_triangle_area;
        % Display results
        fprintf('Mean Triangle Area: %.4f\n', mean_triangle_area);
    end

end

for subji = 31:50
    elech2 = all_EL_Paris{subji}; 
    if ~isempty(elech2) & size(elech2, 1) > 3
        electrode_coords = elech2;
        [K, volume] = convhull(electrode_coords(:,1), electrode_coords(:,2), electrode_coords(:,3));
        num_triangles = size(K, 1);
        triangle_areas = zeros(num_triangles, 1);
        
        for i = 1:num_triangles
            % Extract three vertices of the triangle
            p1 = electrode_coords(K(i,1), :);
            p2 = electrode_coords(K(i,2), :);
            p3 = electrode_coords(K(i,3), :);
        
            % Compute the area using the cross product formula in 3D
            v1 = p2 - p1;
            v2 = p3 - p1;
            triangle_areas(i) = 0.5 * norm(cross(v1, v2));
        end
        
        % Compute the mean triangle area
        mean_triangle_area = mean(triangle_areas);
        
        mTA_LP(subji, :) = mean_triangle_area;
        % Display results
        fprintf('Mean Triangle Area: %.4f\n', mean_triangle_area);

    end
end

for subji = 31:50
    elech2 = all_ER_Paris{subji}; 
    if ~isempty(elech2) & size(elech2, 1) > 3
        electrode_coords = elech2;
        [K, volume] = convhull(electrode_coords(:,1), electrode_coords(:,2), electrode_coords(:,3));
        num_triangles = size(K, 1);
        triangle_areas = zeros(num_triangles, 1);
        
        for i = 1:num_triangles
            % Extract three vertices of the triangle
            p1 = electrode_coords(K(i,1), :);
            p2 = electrode_coords(K(i,2), :);
            p3 = electrode_coords(K(i,3), :);
        
            % Compute the area using the cross product formula in 3D
            v1 = p2 - p1;
            v2 = p3 - p1;
            triangle_areas(i) = 0.5 * norm(cross(v1, v2));
        end
        
        % Compute the mean triangle area
        mean_triangle_area = mean(triangle_areas);
        
        mTA_RP(subji, :) = mean_triangle_area;
        % Display results
        fprintf('Mean Triangle Area: %.4f\n', mean_triangle_area);
    end
end

mTA_LG(mTA_LG==0) = nan; 
mTA_RG(mTA_RG==0) = nan; 
mTA_LP(mTA_LP==0) = nan; 
mTA_RP(mTA_RP==0) = nan; 

mTA_RG(end+1:30) = nan; 
mTA_LP(end+1:50) = nan; 


%%
mTA_LG(mTA_LG==0) = []; 
mTA_RG(mTA_RG==0) = []; 
mTA_LP(mTA_LP==0) = []; 
mTA_RP(mTA_RP==0) = []; 



%% plot 4 bar
clc
data = {mTA_LG, mTA_RG, mTA_LP, mTA_RP}; 
sizeData = cellfun(@length, data)

for condi = 1:4
    
    mean_S = mean(data{condi}, 1, 'omitnan');
    hb = plot ([condi], data{condi}, 'k'); hold on;
    set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',15);hold on;
    h = bar (condi, mean_S);hold on;
    set(h,'FaceColor', 'none', 'lineWidth', 1);
end
set(gca, 'xtick', [1:4], 'xticklabels', {'Guangzhou Left', 'Guangzhou Right', 'Paris Left', 'Paris Right'})
ylabel('Mean triangle area')
set(gca,'FontSize', 14, 'linew',1, 'xlim', [0 5]);
exportgraphics(gca, 'myP.png', 'Resolution',300)

%% combine 
mTA_G = mean([mTA_LG mTA_RG], 2, 'omitnan'); 
mTA_P = mean([mTA_LP mTA_RP], 2, 'omitnan'); 

%% plot 4 bar
clc
data = {mTA_G, mTA_P}; 
sizeData = cellfun(@length, data)

for condi = 1:2
    
    mean_S = mean(data{condi}, 1, 'omitnan');
    hb = plot ([condi], data{condi}, 'k'); hold on;
    set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',15);hold on;
    h = bar (condi, mean_S);hold on;
    set(h,'FaceColor', 'none', 'lineWidth', 1);
end
set(gca, 'xtick', [1:2], 'xticklabels', {'Guangzhou Left', 'Guangzhou Right', 'Paris Left', 'Paris Right'})
ylabel('Mean triangle area')
set(gca,'FontSize', 14, 'linew',1, 'xlim', [0 3]);
exportgraphics(gca, 'myP.png', 'Resolution',300)

[h p ci ts] = ttest2(mTA_G, mTA_P)

%% Two-Factor ANOVA with Unequal Group Sizes

factorA = [zeros(1, length(mTA_LG) + length(mTA_RG)) ones(1, length(mTA_LP) + length(mTA_RP))]; 
factorB = [zeros(1, length(mTA_LG)) ones(1, length(mTA_RG)) zeros(1, length(mTA_LP)) ones(1, length(mTA_RP))]; 
allData = [mTA_LG; mTA_RG; mTA_LP; mTA_RP];

% Run Two-Way ANOVA (Type III SS)
[p, tbl, stats] = anovan(allData, {factorA, factorB}, 'model', 'interaction', ...
    'varnames', {'FactorA', 'FactorB'});

% Display results
disp(tbl);

% Post-hoc multiple comparisons
multcompare(stats, 'Dimension', 1);  % Compare Factor A
multcompare(stats, 'Dimension', 2);  % Compare Factor B









%% EXAMPLE DIFFERENCE BETWEEN INTERELECTRODE DIFFERENCES AND DELAUNAY Define Two Electrode Configurations
% Configuration 1: Uniform Grid (Evenly spaced)
uniform_coords = [0 0; 0 1; 1 0; 1 1; 0.5 0.5; 0.5 1.5; 1.5 0.5; 1.5 1.5];

% Configuration 2: Clustered (Some regions denser than others)
clustered_coords = [0 0; 0.2 0.1; 0.1 0.2; 0.3 0.1; 1 1; 1.2 1; 1 1.2; 1.1 1.1];

%%Compute Mean Interelectrode Distance for Both Configurations
mean_dist_uniform = mean(pdist(uniform_coords));
mean_dist_clustered = mean(pdist(clustered_coords));

%%Compute Delaunay Triangulation and Mean Triangle Area
tri_uniform = delaunay(uniform_coords(:,1), uniform_coords(:,2));
tri_clustered = delaunay(clustered_coords(:,1), clustered_coords(:,2));

% Function to compute triangle areas
triangle_area = @(p1, p2, p3) 0.5 * abs(det([p1 - p3; p2 - p3]));

% Compute mean triangle area for each configuration
areas_uniform = arrayfun(@(i) ...
    triangle_area(uniform_coords(tri_uniform(i,1),:), ...
                  uniform_coords(tri_uniform(i,2),:), ...
                  uniform_coords(tri_uniform(i,3),:)), ...
    1:size(tri_uniform,1));

areas_clustered = arrayfun(@(i) ...
    triangle_area(clustered_coords(tri_clustered(i,1),:), ...
                  clustered_coords(tri_clustered(i,2),:), ...
                  clustered_coords(tri_clustered(i,3),:)), ...
    1:size(tri_clustered,1));

mean_area_uniform = mean(areas_uniform);
mean_area_clustered = mean(areas_clustered);

%%Display Results
fprintf('Mean Interelectrode Distance:\n');
fprintf('  Uniform: %.4f\n', mean_dist_uniform);
fprintf('  Clustered: %.4f\n\n', mean_dist_clustered);

fprintf('Mean Delaunay Triangle Area:\n');
fprintf('  Uniform: %.4f\n', mean_area_uniform);
fprintf('  Clustered: %.4f\n\n', mean_area_clustered);

%%Visualization
figure;
subplot(1,2,1);
triplot(tri_uniform, uniform_coords(:,1), uniform_coords(:,2), 'b');
hold on; scatter(uniform_coords(:,1), uniform_coords(:,2), 50, 'r', 'filled');
title('Uniform Electrode Placement');
axis equal; grid on;

subplot(1,2,2);
triplot(tri_clustered, clustered_coords(:,1), clustered_coords(:,2), 'b');
hold on; scatter(clustered_coords(:,1), clustered_coords(:,2), 50, 'r', 'filled');
title('Clustered Electrode Placement');
axis equal; grid on;


clc; clear; close all;

%% Define 3D electrode coordinates (Replace with your actual data)
electrode_coords = [
    65.3000  -27.3500  -14.1500;
    34.0500   -2.6500  -43.9500;
    39.5500   -2.5000  -43.8000;
    45.0000   -2.8500  -43.3500;
    30.4000  -32.1000  -20.7500;
    35.6000  -29.5000  -22.2500;
    45.4000  -25.0000  -25.5500;
    50.5500  -22.4500  -27.0500;
    55.7000  -19.9000  -28.3500;
    60.8500  -17.9000  -29.5500;
    46.4000  -49.9000   -7.4500;
    62.7500  -47.9000   -1.9500
];

%%Compute 3D Delaunay Triangulation
DT = delaunayTriangulation(electrode_coords);

%%Extract Triangles from the Convex Hull
K = convexHull(DT);  % Get the surface triangles

%%Compute Triangle Areas
numTriangles = size(K, 1);
triangleAreas = zeros(numTriangles, 1);

for i = 1:numTriangles
    % Get the 3 points of the triangle
    p1 = electrode_coords(K(i,1), :);
    p2 = electrode_coords(K(i,2), :);
    p3 = electrode_coords(K(i,3), :);
    
    % Compute triangle area using the cross product
    v1 = p2 - p1;
    v2 = p3 - p1;
    triangleAreas(i) = 0.5 * norm(cross(v1, v2));  % Area formula
end

%%Compute Mean Triangle Area
meanTriangleArea = mean(triangleAreas);

%%Display Results
fprintf('Mean Triangle Area: %.4f\n', meanTriangleArea);

%%Plot the 3D Delaunay Triangulation
figure;
trisurf(K, electrode_coords(:,1), electrode_coords(:,2), electrode_coords(:,3), ...
    'FaceColor', 'cyan', 'FaceAlpha', 0.7, 'EdgeColor', 'k');
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Delaunay Triangulation');
axis equal;
grid on;


%% COUNT SUBREGIONS in each ROI

clear, clc

paths = load_paths_EXT; 
file2load = ['allS_' 'PFC' '_C']; 

load ([paths.results.power file2load]); 


%%
clearvars -except ALLEEG paths file2load
clc 
all_chans = []; 
all_sub = []; 
for subji = 1:length(ALLEEG)


    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        all_chans = [all_chans ; chans(:)]; 
        all_sub = [all_sub ; repelem(subji, length(chans))']; 
    end
end

%%
all_chans2 = [all_chans num2cell(all_sub)]; 


%%
all_chans2(:, 1) = strrep(all_chans2(:, 1), 'ctx-lh-', '')
all_chans2(:, 1) = strrep(all_chans2(:, 1), 'ctx-rh-', '')

%%
clc

% Find unique entries and index mapping
[uniqueEntries, ~, idx] = unique(all_chans2(:, 1));

%uniqueEntries = strrep(uniqueEntries, 'ctx-', '');
for stri = 1:length(uniqueEntries)
    uniqueEntries{stri} = [ upper(uniqueEntries{stri}(1)) uniqueEntries{stri}(2:end)];
end

% Count the occurrences
counts = accumarray(idx, 1);

figure; set(gcf, 'position', [100 100 1300 700])
p = piechart(counts, uniqueEntries,LabelStyle="name")
%p = piechart(counts, uniqueEntries,LabelStyle="namedata")
%p.Labels = uniqueEntries + "\newline"  + string(counts) ;
p.FontSize = 30; 
p.NumDisplayWedges = 6; % group in "others" those that are too small
p.DisplayOrder = "descend";
myCArray = brewermap(9, 'Set3'); 
colororder(myCArray); 

exportgraphics(gca, 'myP.png', 'Resolution',300)

%% 
clearvars -except ALLEEG all_chans2
clc
all_chans3 = string(all_chans2);
% Get Unique Brain Region Names
uniqueRegions = unique(all_chans3(:, 1));

% Initialize Output Arrays
regionCounts = zeros(length(uniqueRegions), 1); % Number of occurrences
uniqueSubjectsCount = zeros(length(uniqueRegions), 1); % Unique subjects count

% Loop Through Each Unique Brain Region
for i = 1:length(uniqueRegions)
    region = uniqueRegions{i};
    
    % Find Rows Matching the Region
    idx = strcmp(all_chans3(:, 1), region);
    
    % Count Occurrences
    regionCounts(i) = sum(idx);
    
    % Count Unique Subject IDs
    uniqueSubjectsCount(i) = numel(unique(all_chans3(idx, 2)));
end

% Convert Results to Table for Readability
uniqueRegions = char(uniqueRegions); 
uniqueRegions(:, 1) = upper(uniqueRegions(:, 1)); 
resultTable = table(char(uniqueRegions), regionCounts, uniqueSubjectsCount, ...
    'VariableNames', {'BrainRegion', 'Electrodes', 'Subjects'});

% Display the Result
disp(resultTable);



%% SELECT CHANNELS TO PLOT
clearvars -except ALLEEG paths file2load
clc 
for subji = 1:length(ALLEEG)
    EEG = ALLEEG{subji};
    if ~isempty(EEG)
        chans = [{EEG.chanlocs.fsLabel}]'; 
        ids2rem1 = []; 
        %ids2rem1 = contains(chans, 'Left')
        ids2rem1 = contains(chans, 'Right')
        %ids2rem1 = contains(chans, 'Hippocampus')
        %ids2rem = logical(ids2rem1+ids2rem2)
        ids2rem = ids2rem1; 
        if ndims(EEG.power) == 3
            tmph(:, 1, :, :)  = EEG.power; 
            EEG.power = tmph; 
        end
        EEG.chanlocs(ids2rem) = []; 
        if size(EEG.chanlocs, 2) >0 & size(EEG.chanlocs, 1) >0 
            EEG.power(:, ids2rem, :, :) = []; 
            ALLEEG1{subji,:} = EEG; 
        end

    end
end

clear ALLLEEG
ALLEEG = ALLEEG1; 

%% Plot time frequency power for all channels in one subject

EEG = ALLEEG{19}; 
size(EEG.power)
nChans = size(EEG.power, 2); 
nTrials = size(EEG.power, 1); 
for chani = 1:nChans
    for triali = 1:5 %:nTrials
        
        figure()
        imagesc(squeeze(EEG.power(triali, chani, :, :))); 

    end
end




%% plot cluster depic schematic 

eE1 = zeros(44,200); 
eE1(clustinfo.PixelIdxList{32}) = 1; 


figure(); set(gcf, 'Position', [100 100 500 200])
contourf(1:2000, 1:440, myresizem(eE1, 10), 40, 'linecolor', 'none'); hold on; 
contour(1:2000, 1:440,myresizem(eE1, 10), 1, 'Color', [0, 0, 0], 'LineWidth', 6); 
plot([250 250 ],get(gca,'ylim'), 'k:','lineWidth', 8);set(gca, 'clim', [-.125 .125]); hold on; 

set(gca, 'xTick', [], 'xTicklabel', [], 'yTick', [], 'yTicklabel', []);
colormap autumn

exportgraphics(gcf, ['myP.png'], 'Resolution',300)







%% compare during acquisition

[h1 p1 ci ts] = ttest(CSPPE, CSPME); 
t1 = ts.tstat; 

[h2 p2 ci ts] = ttest(CSPPE, CSMME); 
t2 = ts.tstat; 

disp (['t = ' num2str(t1) '  ' ' p = ' num2str(p1)]);
disp (['t = ' num2str(t2) '  ' ' p = ' num2str(p2)]);



%% EXTRACT SINGLE TRIAL POWER LEVELS

clearvars -except ALLEEG paths clustinfo nL 
close all, clc
sub2exc = []; %subj 23-35-38-39-44-45

load clustinfoE_px4

for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2');
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        %length(find(ids))
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:500), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            %cTR = squeeze(mean(powH(triali, :, 3:8, 441:460), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
            nNan{subji,:} = sum(isnan(thPow));

           allPOWAMY{subji, :} = thPow; 


    end

    

end


save ('allPOWAMY', 'allPOWAMY', 'nNan'); 
 


%% Correlate responses and power in cluster

clearvars -except ALLEEG paths clustinfo nL 
close all, clc
sub2exc = []; %subj 23-35-38-39-44-45

load clustinfoE_px4


for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        % % %   % % Acquisition
        %ids = strcmp(Ev2(:, 2), '1'); 
        %ids = strcmp(Ev2(:, 2), '1') &  ( strcmp(Ev2(:, 6), '1')  | strcmp(Ev2(:, 6), '2') ) ;
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3');
    
        % % % % % % Extinction
        ids = strcmp(Ev2(:, 2), '2'); 
        %ids = strcmp(Ev2(:, 2), '2') & strcmp(Ev2(:, 6), '1') ;
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2') ) ; 
        %ids = strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '3') ) ; 
        

        % % % Acquisition AND Extinction
        %ids = strcmp(Ev2(:, 2), '1') | strcmp(Ev2(:, 2), '2');
        %ids = strcmp(Ev2(:, 2), '1') & strcmp(Ev2(:, 6), '3') | ( strcmp(Ev2(:, 2), '2') & ( strcmp(Ev2(:, 6), '2')  | strcmp(Ev2(:, 6), '3') ) );

        %length(find(ids))
        powH = EEG.power(ids, :, :, :);

        
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:500), 2));
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{5}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');

            
            %cTR = squeeze(mean(powH(triali, :, 3:8, 441:460), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        nNan(subji,:) = sum(isnan(thPow));

        ratings2u = double(string(Ev2(ids, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        allNTR(subji,:) = length(thPow);



        allRho(subji, :) = corr(thPow, ratings2u, 'type', 'k');
        
% % % %         
% % % %         figure()
% % % %         scatter(thPow, ratings2u, 150, 'filled');
% % % %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% % % %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% % % %         allSlopes(subji, :) = C(2);
% % % %         allIntercepts(subji, :) = C(1);
% % % %         set(gca, 'ylim', [1 5], 'xlim', [-4 4], 'Fontsize', 24)
% % % %         


    end

    

end
 
idF = allRho==0 | isnan(allRho); 
allRho(idF) = []; 
%allSlopes(idF) = []; 

[h p ci t] = ttest (allRho);
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);



%% plot one bar
data.data = [allRho]; 


figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 2], 'ylim', [-.3 .5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);
set(gca, 'LineWidth', 3);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)










%% ALL FREQUENCIES + ERPs

load allERPs
sub2exc = [];
c1B = c1(:, 1:54, 251:480); c2B = c2(:, 1:54, 251:480); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
d2p1	= squeeze(mean(c1B, 'omitnan'));
d2p2	= squeeze(mean(c2B, 'omitnan'));


[h p ci ts] = ttest(c1B, c2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 

% 

%h = zeros(54, 230);
%h(clustinfo.PixelIdxList{id}) = 1; 

 

c2p1 = allERPs(:,1);
c2pm = cellfun(@(x) mean(x, 1), c2p1, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm1 = cat(1, c2pm{:});
cp1m = downsample(c2pm1', 10)'; 
cp1m = cp1m(:,  251:480); 
d2p1ERP = mean(cp1m); 

c2p2 = allERPs(:, 2);
c2pm = cellfun(@(x) mean(x, 1), c2p2, 'un', 0);
c2pm(cellfun('isempty', c2pm))=[]; 
c2pm2 = cat(1, c2pm{:});
cp2m = downsample(c2pm2', 10)'; 
cp2m = cp2m(:,  251:480); 
d2p2ERP = mean(cp2m); 


times = -.5:.01:1.79; 
%times = -3:.01:3.99
freqs = 1:size(c1B, 2);
figure()
tiledlayout(3, 1,'TileSpacing','loose'); set(gcf, 'Position', [100 100 600 900])
nexttile
contourf(times, freqs, d2p1, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);set(gca, 'clim', [-.125 .125]); hold on; 
yyaxis right
plot(times, d2p1ERP, 'k', LineWidth=2);
set(gca, ylim=[-2 2])
title('CS+')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


nexttile
contourf(times, freqs, d2p2, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'Z-Power';
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3); set(gca, 'clim', [-.125 .125])
yyaxis right
plot(times, d2p2ERP, 'k',LineWidth=2); 
title('CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')


nexttile
contourf(times, freqs, t, 40, 'linecolor', 'none'); hold on; c = colorbar; c.Label.String = 'T';
contour(times, freqs,h, 1, 'Color', [0, 0, 0], 'LineWidth', 2); set(gca, 'clim', [-4 4])
plot([0 0 ],get(gca,'ylim'), 'k:','lineWidth', 3);
yyaxis right
plot(times, d2p1ERP-d2p2ERP,'k', LineWidth=2); 
title('CS+ vs CS-')
xlabel('Time (s)')
ylabel('Frequency (Hz)')
%plot([1.77 1.77 ],get(gca,'ylim'), 'k:','lineWidth', 3);
colormap(brewermap([],'*Spectral'))

set(findobj(gcf,'type','axes'),'FontSize',16, 'ytick', [1 30 54], 'yticklabels', {'1', '30', '150'}, 'xlim', [-.5 1.75]);
%set(findobj(gcf,'type','axes'),'FontSize',18, 'ytick', [3 8], 'yticklabels', {'3', '8'}, 'xlim', [-.5 1.75]);


%exportgraphics(gcf, [paths.results.power file2load '.png'], 'Resolution',150)
exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',300)














%% THETA BAND - LINE PLOTS

sub2exc = [];

c1B = c1(:, 3:8, 201:500); c2B = c2(:, 3:8, 201:500); 
c1B(sub2exc,:,:) = []; c2B(sub2exc,:,:) = []; 

c1B(c1B == 0) = nan; 
c2B(c2B == 0) = nan; 
mc1B = squeeze(mean(c1B, 2, 'omitnan'));
mc2B = squeeze(mean(c2B, 2, 'omitnan'));
mc1B(any(isnan(mc1B), 2), :) = [];
mc2B(any(isnan(mc2B), 2), :) = [];

d2pm1	= squeeze(mean(mc1B,'omitnan'));
d2pm2	= squeeze(mean(mc2B,'omitnan'));
d2pstd1	= std(mc1B);
d2pstd2	= std(mc2B);
se1 = d2pstd1/sqrt(size(mc1B, 1))
se2 = d2pstd2/sqrt(size(mc2B, 1))

[h p ci ts] = ttest(mc1B, mc2B); 
h = squeeze(h); t = squeeze(ts.tstat);

clear allSTs  
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id); 

h(1:150) = 0;

hb = h; hb(h==0) = nan; hb(hb==1) = -.01; 
times = (-1:.01:1.99) + .25;

colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 

%plot(times, d2p1); hold on; 
%plot(times, d2p2); hold on; 

xlabel('Time (s)')
ylabel('Theta Power')
plot (times, hb, 'Linewidth', 7)
set(gca, 'xlim', [-.25 1.75])
set(gca, 'FontSize', 24);

exportgraphics(gcf, [paths.results.power  'myP.png'], 'Resolution',150)




%% permutations 2D (line plot)

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    c1B = squeeze(mc1B(:, 51:200)); 
    c2B = squeeze(mc2B(:, 51:200));
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :);
           c1B(subji, :) = c2B(subji, :);
           c2B(subji, :) = tmp; 
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
        [max2u id] = min(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

%%

clear p ratings2u mcsP

ratings2u = max_clust_obs; 
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





%% Extract data 4 LME


clearvars -except ALLEEG paths  
close all
load clustinfoE_px5
%pi2u = unique([clustinfo_A.PixelIdxList{3} ; clustinfo_E.PixelIdxList{3} ]);
sub2exc = [];

for subji = 1:size(ALLEEG, 1)
    EEG = ALLEEG{subji}; 
    if ~isempty(EEG) & isempty(intersect(sub2exc, subji))
        Ev = [{EEG.event.type}]';
        Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
        Ev2 = cat(1, Ev1{:});
        Ev2(:, 10) = erase(Ev2(:, 10), ' '); 
    
        powH = EEG.power;
        
        for triali = 1:size(powH, 1)
            cTR = squeeze(mean(powH(triali, :, 3:8, 251:480), 2));
            thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{5}), 'all');
            %thPow(triali, :) = mean(cTR(clustinfo.PixelIdxList{4}), 'all');
            %thPow(triali, :) = mean(cTR(pi2u), 'all');
            
            %cTR = squeeze(mean(powH(triali, :, 3:8, 301:400), 2));
            %thPow(triali, :) = mean(cTR, 'all');
    
        end
        
        ratings2u = double(string(Ev2(:, 7)));

        ids2rem = isnan(thPow) | isnan(ratings2u); 
        trialN = double(string(Ev2(:, 1)));
        exph = double(string(Ev2(:, 2)));
        currCS = double(string(Ev2(:, 8)));
        trial_type = double(string(Ev2(:, 6)));

        trialN(ids2rem) = []; 
        trial_type(ids2rem) = []; 
        currCS(ids2rem) = []; 
        exph(ids2rem) = []; 
        thPow(ids2rem) = []; 
        ratings2u(ids2rem) = []; 
        subID = repelem(subji, length(ratings2u))';
        
        
        d4LME = [thPow ratings2u subID trialN exph currCS trial_type];
        % % remove the testing phase for now
        %d4LME(d4LME(:, 5) == 3,:) = []; 

        allData4LME{subji,:} = d4LME;
        
    end
end


%% 

%%fit model
clc

d4LME = cat(1, allData4LME{:});
d4LME(d4LME(:, 2) == 5,:) = []; % % remove ratings of 5 !
d4LME(d4LME(:, 5) == 3,:) = []; 

tbl = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), d4LME(:,5), d4LME(:,6), d4LME(:,7),...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'Phase', 'currCS', 'trial_type'});

%tbl.Ratings = ordinal(tbl.Ratings);
%tbl2.currCS= categorical(tbl2.currCS);
%tbl.trial_type= categorical(tbl.trial_type);
%tbl2.Phase= categorical(tbl2.Phase);


%VarDecompTbl = colldiag(d4LME)
% That can be passed along for visualization
%colldiag_tableplot(VarDecompTbl);



%lme = fitlme(tbl2,'Ratings ~ theta_AMY + trial_type + Phase + trialN + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'Ratings ~ theta_AMY + currCS+ Phase+ trialN + currCS*theta_AMY + Phase*theta_AMY + (1|subID)'); % random intercept model
%lme = fitlme(tbl2,'theta_AMY ~ Ratings + currCS + Phase+ trialN + (1|subID)'); % random intercept model

lme = fitlme(tbl,'theta_AMY ~ Ratings  + trial_type + Phase*currCS + (1|subID)'); % random intercept model

lme
%lme.Coefficients


%%
writetable(tbl, 'd4LME.csv')
writetable(d4LME2, 'd4LME2.csv')


%% 
plotResiduals(lme,'fitted', 'ResidualType', 'Raw')

%plotResiduals(lme,'lagged')



%% FIT MODEL ONLY during extinction
clc

d4LME2 = cat(1, allData4LME{:});

% % % EXTINCTION
d4LME2 = d4LME2(d4LME2(:, 5) == 2, :); 
tbl2 = table(d4LME2(:,1), d4LME2(:,2), d4LME2(:,3), d4LME2(:,4),d4LME2(:,6), d4LME2(:,7), ...
    'VariableNames',{'theta_AMY','Ratings','subID', 'trialN', 'currCS', 'trial_type'});
lme = fitlme(tbl2,'theta_AMY ~ Ratings + trialN  + currCS + currCS*Ratings + (1|subID)'); % random intercept model

% % % % ACQUISITION
%d4LME2 = d4LME2(d4LME2(:, 5) == 1 & d4LME2(:, 6) == 0, :); % only cs-
%tbl2 = table(d4LME2(:,1), d4LME2(:,2), d4LME2(:,3), d4LME2(:,4), 'VariableNames',{'theta_AMY','Ratings','subID', 'trialN'});
%lme = fitlme(tbl2,'theta_AMY ~ Ratings + trialN  + (1|subID)'); % random intercept model


lme







%% 
writetable(tbl2, 'extinction_data', 'FileType','spreadsheet')





%% plot 2 bar
mcsP(mcsP==0) = nan; mcsM(mcsM==0) = nan;
data.data = [mcsP mcsM ];

figure(2); set(gcf,'Position', [0 0 500 650]); 
mean_S = mean(data.data, 1, 'omitnan');
hb = plot ([1 2], data.data); hold on;
set(hb, 'lineWidth', 3, 'Marker', '.', 'MarkerSize',35);hold on;
h = bar (mean_S);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 3);
set(gca,'XTick',[1 2],'XTickLabel',{'', ''}, 'FontSize', 30, 'linew',2, 'xlim', [0 3], 'ylim', [0 5] );
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 3);

[h p ci t] = ttest (data.data(:,1), data.data(:,2));
disp (['t = ' num2str(t.tstat) '  ' ' p = ' num2str(p)]);

set(gca, 'LineWidth', 3);

%export_fig(2, '_2.png','-transparent', '-r80');
%close all;   


%non-learners
nL = mcsP > mcsM; 







%% count channels 

clear nChans
for subji = 1:length(ALLEEG1)

    EEG = ALLEEG1{subji};
    if~isempty(EEG)
        nChans(subji,:) = length(EEG.chanlocs);
    end


end

disp (['Subjects with elec: ' num2str(length(find(nChans)))])
disp (['Total elec num: ' num2str(sum(nChans))])







%% 
d2p = squeeze(EEG.data(1, :, 1)); 
figure;plot(d2p)



%%
data2check = [EEG.data(1, :); EEG.marker_artifacts(1,:)*1000]; 
eegplot(data2check, 'srate', EEG.srate, 'winlength', 50, 'spacing', 1000);










%% plot example trial in one subject (WITH MORE THAN 1 electrode)

EEG = ALLEEG{37}; 
tr = 20; 
ch = 1; 

% figure
% d2p	= squeeze(EEG.dsPower(tr, ch, : ,:));
% myCmap = colormap(brewermap([],'YlOrRd'));
% colormap(myCmap)
% contourf(1:300, 1:54, d2p, 40, 'linecolor', 'none'); colorbar

figure
myCmap = colormap(brewermap([],'YlOrRd'));
colormap(myCmap)
d2p	= squeeze(EEG.power(tr, ch, : ,:));
contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar



%% plot two different conditions (average in all amygdala channels and across trials)


Ev = [{EEG.event.type}]';
Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
Ev2 = cat(1, Ev1{:})
ids = strcmp(Ev2(:, 6), '1')  & strcmp(Ev2(:, 2), '1') % CS+CS+ during acquisition
d2p1	= squeeze(mean(mean(EEG.power(ids, :, : ,:), 'omitnan'), 'omitnan'));
ids = strcmp(Ev2(:, 6), '2')  & strcmp(Ev2(:, 2), '1') % CS+CS- during acquisition
d2p2	= squeeze(mean(mean(EEG.power(ids, :, : ,:), 'omitnan'), 'omitnan'));

myCmap = colormap(brewermap([],'YlOrRd')); colormap(myCmap)

tiledlayout(2, 1,'TileSpacing','compact');
nexttile
contourf(1:700, 1:54, d2p1, 40, 'linecolor', 'none'); hold on; %colorbar
nexttile
contourf(1:700, 1:54, d2p2, 40, 'linecolor', 'none'); hold on; %colorbar



%% plot all electrodes and all trials 

clear, clc

paths = load_paths_EXT; 


%files2load = {['allS_' 'PFC' '_C']; ['allS_' 'AMY' '_C']; ['allS_' 'TMP' '_C']; ['allS_' 'HPC' '_C']; ['allS_' 'OCC' '_C']; ['allS_' 'OFC' '_C']}; 
files2load = {['allS_' 'TMP' '_C']}; 


for listi = 1:length(files2load)

    clearvars -except listi files2load paths

    file2load = files2load{listi}; 

    load ([paths.results.power file2load]); 

    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji}; 
        if ~isempty(EEG)
            nChans = size(EEG.chanlocs,1); 
            nTrials = size(EEG.power, 1); 
        
            for chani = 1:nChans
                for triali = 1%nTrials
                    figure
                    myCmap = colormap(brewermap([],'YlOrRd'));
                    colormap(myCmap)
                    d2p	= squeeze(EEG.power(triali, chani, : ,:));
                    contourf(1:700, 1:54, d2p, 40, 'linecolor', 'none'); colorbar
                    chName = EEG.chanlocs(chani).labels; 
                    fname = [paths.results.allSpectro file2load(6:end) '_' num2str(subji, '%02.f') '_' chName '_' num2str(triali, '%03.f') '.png']; 
                    exportgraphics(gcf, fname, 'Resolution',100)
                    close all; 
                end
            end
        end
    end
end



%%
c = [115, 88, 86, 104, 87, 108, 96, 94, 95, 83, 83, 107, 92, 88, 98, 98, 91, 93, 93, 124, 84, NaN, ...
    82, 75, NaN, NaN, 97, 82, NaN, 109, NaN, NaN, NaN, 75, 118, NaN, NaN, NaN, NaN, 93, 101, 81, NaN, 99, NaN, 88, NaN, NaN, NaN, NaN];

c = [106, 39, 73, 80, 51, 75, 81, 52, 70, 77, 51, 95, 63, 51, 81, 73, 89, 84, 108, 120, 52, NaN, 66, 68, NaN, NaN, 54, 55, NaN, 103];

mean(c, 'omitnan')
std(c, 'omitnan')


%% count number of trials for each condition in LOOP
clear, clc
paths = load_paths_EXT; 
allFiles = { ['allS_' 'AMY' '_C'],  ['allS_' 'TMP' '_C'],  ['allS_' 'PFC' '_C'],  ['allS_' 'HPC' '_C'],  ['allS_' 'OFC' '_C']}
for filei = 1:5
    file2load =allFiles{filei}; 
    
    load ([paths.results.power file2load]); 
    
    %%count number of electrodes in Guangzhou and Paris
    clearvars -except ALLEEG paths file2load filei data allFiles
    clc 
    all_chansChina = []; 
    all_chansParis = []; 
    allPatsChina = []; 
    allPatsParis = []; 
    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji};
        if ~isempty(EEG)
            if subji < 31
                chansChina = [{EEG.chanlocs.fsLabel}]'; 
                all_chansChina = [all_chansChina ; chansChina(:)]
                allPatsChina = [allPatsChina subji]; 
            else
                chansParis = [{EEG.chanlocs.fsLabel}]'; 
                all_chansParis = [all_chansParis ; chansParis(:)]
                allPatsParis = [allPatsParis subji]; 
            end
        end
    end
    
    data(filei, :) = [length(allPatsChina) length(allPatsParis)]; 
end

%% barplot
% Plot with a pie chart
bar(1:5, data, 'stacked');
%xlabel('Entries');
ylabel('Patients');
set(gca, 'xticklabels', {'AMY', 'TMP', 'PFC', 'HPC', 'OFC'}); 
title('Patients with > 1 electrode in each ROI');
legend('Guangzhou', 'Paris')

exportgraphics(gca, 'myP.png', 'Resolution',300)

%%
% Example Data
observed_counts = (data(5,:)); % 30 subjects with electrodes, 70 without
expected_proportions = [0.5, 0.5]; % Null hypothesis: Equal chance

% Convert proportions to expected counts
expected_counts = sum(observed_counts) * expected_proportions;

% Perform chi-square goodness-of-fit test
[h, p, stats] = chi2gof(observed_counts, 'Expected', expected_counts);

% Display results
disp(['Chi-square statistic: ', num2str(stats.chi2stat)]);
disp(['p-value: ', num2str(p)]);
if h == 1
    disp('Reject the null hypothesis: The proportion is significantly different from chance.');
else
    disp('Fail to reject the null hypothesis: No significant difference.');
end


%% SELECT ALL ELECTRODES IN ALL SUBJECTS
clear, clc
paths = load_paths_EXT; 
allFiles = { ['allS_' 'AMY' '_C'],  ['allS_' 'TMP' '_C'],  ['allS_' 'PFC' '_C'],  ['allS_' 'HPC' '_C'],  ['allS_' 'OFC' '_C']}
all_chans = []; 
all_sub = []; 
for filei = 1:5
    file2load =allFiles{filei}; 
    
    load ([paths.results.power file2load]); 
    
    %%count number of electrodes in Guangzhou and Paris
    clearvars -except ALLEEG paths file2load filei data allFiles all_chans all_sub
    
    for subji = 1:length(ALLEEG)
        EEG = ALLEEG{subji};
        if ~isempty(EEG)
            chans= EEG.chanlocs; 
            all_chans = [all_chans; chans(:)]
            all_sub = [all_sub ; repelem(subji, length(chans))']; 
        end
    end
   
end

allChans = struct2cell(all_chans)';
allChans = [allChans num2cell(all_sub)]; 

%% sort by subject

ids2sort = [allChans{:, 7}]; 
[id1 id2 ] = sort(ids2sort); 
allChansSorted = allChans(id2, :); 
%%
allCHS = allChansSorted(:, [1 2 4 7]); 
allCHST = cell2table(allCHS, 'VariableNames', {'Channel', 'Label', 'MNI Coordinates', 'Subjects'});
allCHST = allCHST(:, {'Subject', 'Channel', 'MNI Coordinates','Label'} ); 

writetable(allCHST, 'electrode_info.csv')

%%


eventChannel = 'POL DC12';TTL = strmatch(eventChannel, {EEG.chanlocs.labels}, 'exact');chanids = [TTL];
eegplot(EEG.data(chanids,:), 'srate', EEG.srate, 'eloc_file',EEG.chanlocs(chanids), ...
    'winlength', 50, 'spacing', 10000000, 'events', EEG.event);


%%

responseEvent = '_V';
myEvents = string({EEG.event.type})'; 
respEV = contains(myEvents, responseEvent);
trialEvent = '_C';
trialEV = contains(myEvents, trialEvent);

tr_events = {EEG.event(trialEV | respEV).latency}'
tr_eventsLabels = {EEG.event(trialEV | respEV).type}'

%tr_events(1, :) = []; 
%tr_eventsLabels(1, :) = []; 


%%
numLate = double(string(tr_events)); 
oneDiff = numLate(5) - numLate(4)
diffE = diff(numLate);

%%
clear, clc
% Define the seizure onset zones for each patient
seizureOnsetZones = { 
            'Right temporal lobe, bilateral medial temporal lobes, insula, and right extra-orbital cortex'
            'Right medial temporal lobe'
            'Left temporo-occipital lobe'
            'Left frontal lobe'
            'Left medial temporal lobe'
            'Left temporal lobe'
            'Right temporal insular cortex'
            'Bilateral medial temporo-occipital junction regions'
            'Right temporo-occipital junction region'
            'Right limbic system'
            'Right temporal lobe'
            'Bilateral limbic system'
            'Right temporo-insular region'
            'Bilateral medial temporal lobes, right hippocampus, and right insula'
            'Left temporo-insular lobe, orbital frontal gyrus, precentral sulcus, and parietal lobe.'
            'Medial temporal lobe combined with insula, opercula, and orbitofrontal region'
            'Right temporo-parietal junction region'
            'Left posterior temporal lobe'
            'Right hippocampus, right insula, right central operculum'
            'Bilateral parieto-occipital regions'
            'Right posterior temporal region and temporo-occipital junction cortex'
            'Left temporo-insular lobe'
            'Right anterior temporal lobe'
            'Right limbic system'
            'Right temporal insular cortex'
            'Right insula, right posterior orbital region and right medial temporal lobe'
            'Left medial temporal lobe and insula'
            'Left prefrontal lobe'
            'Left temporo-insular lobe'
            'Left temporal lobe'
            'Left temporal lobe'
            'Bilateral temporo-mesial'
            'Left temporal lobe'
            'Right temporo-parietal'
            'Right temporal lobe'
            'Left temporal lobe'
            'undetermined focus'
            'Temporo-polar'
            'Right temporal'
            'Bilateral temporo-basal'
            'Bilateral temporo-mesial'
            'Bilateral temporo-mesial'
            'Bilateral temporo-mesial'
            'Left temporo-mesial'
            'Bilateral temporo-mesial'
            'Polar and baso-medial region of the right temporal lobe'
            'Left temporo-mesial'
            'Bilateral temporo-mesial'
            'Right temporo-mesial'

};

totalPatients = length(seizureOnsetZones);
% Count the number of patients with electrodes in the right hemisphere
myCount = 0;
for i = 1:length(seizureOnsetZones)
    if contains(seizureOnsetZones{i}, 'insul', 'IgnoreCase', true)
        myCount = myCount + 1;
    end
end

%myCount = 13; 
myPerc = (myCount / totalPatients) * 100;

% Display the result
fprintf('Number of patients with electrodes in the right hemisphere: %d\n', myCount);
fprintf('Percentage of patients with electrodes in the right hemisphere: %.2f%%\n', myPerc);




%% load electrode localizations Paris and compare with Guangzhou

clear, clc
allElec = readcell('C:\Users\Neuropsychology\Desktop\extinction_project\fromSara\c_sub06\c_sub06_anatomical_localizationsOK.csv','Delimiter', ';')


%%
C = allElec; 
C_filled = C;
for i = 2:size(C_filled, 1)
    if ismissing(C_filled{i,1})
        C_filled{i,1} = C_filled{i-1,1};
    end
end

%%

C_filled(:, 1) = strrep(C_filled(:, 1), '_RH', ''''); 
C_filled1 = C_filled(:, [1 2 36])
count =1; 
for eleci = 1:154

    if ~ismissing(C_filled1{eleci, 2})
        tmp = join([string(C_filled1{eleci, 1}) string(C_filled1{eleci, 2})]); 
        tmp = strrep(tmp, ' ', ''); 
        myElecParis{count, 1} = tmp; 
        myElecParis{count, 2} = C_filled1{eleci, 3}; 
        count = count+1; 
    end
end


%% find equivalent in Guangzhou
load ('D:\extinction\raw_data\c_sub06\elec\g_wenxiyin_190522_Electrode.mat')

myElecGuangzhou = table2cell(elec_info); 
myElecGuangzhou = string(myElecGuangzhou(:, 1:2)); 



%% take only the same numbers 
clc
A = [myElecParis{:,1}]'; 
B = myElecGuangzhou(:,1); 
[C ids1 ids2] = intersect(A, B)

myElecParis = myElecParis (ids1, :); 
myElecGuangzhou = myElecGuangzhou(ids2, :); 



%% same for column 2
clc
A = myElecParis(:,2); 
B = myElecGuangzhou(:,2); 
sum(A==B)


%% now only with electrodes included in the analyses
allROIs = {'Amygdala' 'Hippocampus' 'orbitofrontal' 'inferiortemporal' 'middletemporal' 'superiortemporal' ...
            'transversetemporal' 'fusiform' 'temporalpole' 'bankssts' 'parahippocampal' 'entorhinal' 'occipital' ...
            '-cuneus' 'lingual' 'pericalcarine' 'caudalmiddlefrontal' 'parsopercularis' 'parsorbitalis' ...
            'superiorfrontal' 'parstriangularis' 'rostralmiddlefrontal' 'frontalpole'}; 
                
selChans = contains(myElecParis(:, 2), allROIs);


myEP_ROI = myElecParis(selChans, :); 
myEG_ROI = myElecGuangzhou(selChans, :); 

%% take only the same numbers 
clc
A = [myEP_ROI{:,1}]'; 
B = myEG_ROI(:,1); 
[C ids1 ids2] = intersect(A, B)

myEP_ROI = myElecParis (ids1, :); 
myEG_ROI = myElecGuangzhou(ids2, :); 


clc
A = myEP_ROI(:,2); 
B = myEG_ROI(:,2); 
sum(A==B)











%%










%%