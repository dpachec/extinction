%% Extinction behavioral analysis
%% INDIVIDUAL SUBJECT PLOTS: separately for each subject and type
%%
% 1. what trial number (position in presentation)?
% 2. which Phase?
% 3. which context was used?
% 4. what was the role of the video (A,B,C1,C2)
% 5. which item was shown?
% 6. which type of item was shown? % cs+/cs+=1;cs+/cs-=2;cs-/cs-=3;
% 7. what response was given?
% 8. cs (0/1) current cs+/cs-
% 9. us 0/1 (y/n)
%%%% SR logfile 10000
% 10. sample point trialonset
% 11. sample point videoonset
% 12. sample point cueonset
% 13. sample point us onset
% 14. sample point response 
clear 

type2u = [3]; % (1, 2, 3) ; Cs+Cs+, Cs+Cs-, Cs-Cs-

types={'cs+cs+','cs+cs-','cs-cs-'};

paths = load_paths_EXT;


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21'}'; % no p_sub08



responses = [];
figure(1); set(gcf, 'Position', [10 20 1850 1350])
tiledlayout(7, 8,'TileSpacing','compact', 'Padding','none');

for subji=1:numel(allsubs)

    nexttile
    sub=allsubs{subji};
    info_file=strcat(paths.trlinfo,sub,'_trlinfo');
    load(info_file)  
    trialinfo = trlinfo; 
    
    label={'Dangerous','Safe'};
    
    numberOfResponses{subji, :}= length(find(~isnan(trialinfo(:,7))));
    type_item=mean(trialinfo(trialinfo(:,5)==type2u,6));
    tmp=trialinfo(trialinfo(:,5)==type2u,7);
    responses=[tmp',nan(1,length(responses)-numel(tmp))];
    tmp_us=trialinfo(trialinfo(:,5)==type2u,9);
    us_vec=[tmp_us',nan(1,length(tmp_us)-numel(tmp_us))];


    responses(responses>4)=4; % sub c7 c21 c30 and p7 have 1 rating (only %1 of value 5
    avg_response_type =mean(responses,2, 'omitnan') 
    
    
    hold on
    rectangle('Position',[0 0 24 4],'FaceColor',[1 .9 .9])
    rectangle('Position',[24 0 24 4],'FaceColor',[0.9 .9 .9])
    rectangle('Position',[48 0 16 4],'FaceColor',[1 .8 1])
    rectangle('Position',[64 0 16 4],'FaceColor',[1 1 .8])
    
    
    stem(responses)
    title(sub, 'Interpreter','none');
    h=gca;
    h.YTick=[0.1;4-(0.2*4)];
    h.YTickLabel=label;
    h.YTickLabelRotation=90;

  
end

  %filename = [paths.results.behavior types{type2u} '_behav.png']
  filename = ['myP.png']
  exportgraphics(gcf, filename, 'Resolution',300)
  

%% AVERAGE REPONSES across TRIALS IN EACH BLOCK 
clear, clc
paths = load_paths_EXT;
path_trlinfo= paths.trlinfo;


% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
%            'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
%            'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
%             'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
%             'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
%             'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
%             'p_sub19', 'p_sub20', 'p_sub21'}'; % no p_sub08

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11', 'c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27', 'c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21'}'; 


for subji=1:numel(allsubs)
    % plot responses for each type across experiment
    load(strcat(path_trlinfo,allsubs{subji},'_trlinfo.mat'))
    trialinfo=trlinfo;
    allTRINFO{subji, :} = trlinfo; 
    nanResponses(subji, :) = length(find(isnan(trialinfo(:, 7))));
    types={'cs+/cs+','cs+/cs-','cs-/cs-'};
    items=1:3;
    
    for i=1:3
        type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
        tmp=trialinfo(trialinfo(:,5)==i,7);
        %responses(i,:)=[tmp',nan(1,64-numel(tmp))]; % does not affect results
        responses(i,:)=tmp;
    end
    
    % average responses for each  n x type x response
    for ty=1:3
        %avg_response_type(subji,ty,:)=nanmean(responses(type_item==ty,:),1) 
        avg_response_type(subji,ty,:)=nanmean(responses(find(type_item==ty),:),1);
    end
end


%%
d2check = squeeze(avg_response_type(37,:,:)); 
length(find(isnan(d2check(1,:))))
length(find(isnan(d2check(2,:))))
length(find(isnan(d2check(3,:))))



%% average rating each block

block_def=[1 24;25 48;49 64];
fig_stuff=subplot(1,1,1);
cmap_default=fig_stuff.ColorOrder;
green= colormap(brewermap([1],'Greens'))
green = green*.9;
red = cmap_default(2,:);
yellow = cmap_default(3,:); 
lCs = {red; yellow; green};
hold  on
for b=1:3
    response_avgblock(:,b)=nanmean(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),1),3);
    response_stdblock(:,b)=(nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1))./sqrt(numel(allsubs));
    response_avgblocksub(:,:,b)=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);
end

for ty=1:3
    errorbar(1:3,response_avgblock(ty,:),response_stdblock(ty,:), 'Color', lCs{ty}, 'LineWidth',2);
end

legend('cs+/cs+','cs+/cs-','cs-/cs-')
h=gca;
h.XLim=[0.5,4.5];
h.XTick=[1:3];
h.XTickLabel={'ACQ','EXT','TEST'}
h.FontSize = 18; 


%filename = [paths.results.behavior 'average_per_type_in_block.png']
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

%% mean across all experiment for CS++ and CS--

threeConditions = nanmean(avg_response_type, 3); 

CSPP = threeConditions(:, 1); 
CSPM = threeConditions(:, 2); 
CSMM = threeConditions(:, 3); 

threeV = [CSPP;  CSPM ; CSMM]

threeConditionsN = normalize(threeV); 

threeV = reshape(threeConditionsN, [], 3)

mean(threeV)
mean(threeConditions)





%% Reshape data 4 Anova
%clc 
sub2exc = [27 37]; % This subject needs to be excluded because of a technical problem, no ratings for the test period
respAVG = response_avgblocksub; 

respAVG(sub2exc, :, :) = []; 
d4anova = respAVG(:);

%rm_anova2
nSubj = size(respAVG, 1); 
subID = repmat((1:nSubj)', 9, 1); 
trial_type = repmat([ones(1,nSubj) , ones(1,nSubj)*2 , ones(1,nSubj)*3]', 3, 1);
block_n = [ones(1,nSubj*3) , ones(1,nSubj*3)*2 , ones(1,nSubj*3)*3]';


allDtog = [subID, trial_type, block_n, d4anova];


anovaStats = rm_anova2(d4anova,subID,trial_type,block_n,{'trial_type', 'block_number'})


%% WITHIN SUBJECTS DESING IN MATLAB
clc
% organize the data in a table
respAN = reshape(respAVG, 48, []); 
T = array2table(respAN);
T.Properties.VariableNames = {'t1b1' 't1b2' 't1b3' 't2b1' 't2b2' 't2b3' 't3b1' 't3b2' 't3b3' };
% create the within-subjects design
withinDesign = table([1 1 1 2 2 2 3 3 3]', [1 2 3 1 2 3 1 2 3]', 'VariableNames',{'block_n', 'trial_type'});
withinDesign.trial_type = categorical(withinDesign.trial_type);
withinDesign.block_n = categorical(withinDesign.block_n);
% create the repeated measures model and do the anova
rm = fitrm(T,'t1b1-t3b3 ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','trial_type*block_n'); % remove comma to see ranova's table
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));


%%
writematrix(allDtog, 'myM.xlsx')


%% Sum of Squares values from ANOVA table
SS_A = 41.8509;
SS_error_A = 103.4369;

SS_B = 5.9429;
SS_error_B = 30.9197;

SS_AB = 6.5836;
SS_error_AB = 44.4520;

% Compute Partial Eta Squared
eta_p2_A = SS_A / (SS_A + SS_error_A);
eta_p2_B = SS_B / (SS_B + SS_error_B);
eta_p2_AB = SS_AB / (SS_AB + SS_error_AB);

% Display results
fprintf('Partial Eta Squared for trial type (Factor A): %.3f\n', eta_p2_A);
fprintf('Partial Eta Squared for block number (Factor B): %.3f\n', eta_p2_B);
fprintf('Partial Eta Squared for Interaction (A × B): %.3f\n', eta_p2_AB);


%% one way anova for each phase separately
clc 
sub2exc = [27 37]; % This subject needs to be excluded because of a technical problem, no ratings for the test period
respAVG = response_avgblocksub; 
respAVG(sub2exc, :, :) = []; 

b2use = 3; 
respAVG1 = squeeze(respAVG(:, :, b2use )); 
data = [respAVG1(:, 1), respAVG1(:, 2), respAVG1(:, 3)];
nSubj = size(respAVG1(:, 1), 1); 
d4anova = data(:);
d4anova = d4anova(:); 
d4anova(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
d4anova(:,3) = [1:nSubj 1:nSubj 1:nSubj];

[allPs allFs] = RMAOV1(d4anova);



%% Reshape data 4 Anova ONLY CS++ vs CS--

clc 
sub2exc = [27 37]; % This subject needs to be excluded because of a technical problem, no ratings for the test period
respAVG = response_avgblocksub(:, [1 3], :); 

respAVG(sub2exc, :, :) = []; 
d4anova = respAVG(:);

%rm_anova2
nSubj = size(respAVG, 1); 
subID = repmat((1:nSubj)', 6, 1); 
trial_type = repmat([ones(1,nSubj), ones(1,nSubj)*3]', 3, 1);
block_n = [ones(1,nSubj*2) , ones(1,nSubj*2)*2 , ones(1,nSubj*2)*3]';



anovaStats = rm_anova2(d4anova,subID,trial_type,block_n,{'trial_type', 'block_number'})


%% Friedman
clc
data = squeeze(response_avgblocksub(:,3,:)); 
data([27 37], :) = []; 
bar(mean(data))

% Run the Friedman test:
[p, tbl, stats] = friedman(data, 1, 'off');
% The second argument (1) indicates that there is one observation per block per treatment.
% The 'off' option suppresses the display of the ANOVA table.

% Display the p-value:
fprintf('The p-value from the Friedman test is: %g\n', p);




%% posthocs


% % % % BLOCK 1
[h1,p1, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,2,1))); %CS+/CS+ vs CS+/CS- B1
t1 = ts.tstat; 
[h2,p2, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,3,1))); %CS+/CS+ vs CS-/CS- B1
t2 = ts.tstat; 
[h3,p3, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,1)),squeeze(response_avgblocksub(:,3,1))); %CS+/CS- vs CS-/CS- B1
t3 = ts.tstat; 

% % % % BLOCK 2
[h4,p4, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,2)),squeeze(response_avgblocksub(:,2,2))); %CS+/CS+ vs CS+/CS- B2
t4 = ts.tstat; 
[h5,p5, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,2)),squeeze(response_avgblocksub(:,3,2))); %CS+/CS+ vs CS-/CS- B2
t5 = ts.tstat; 
[h6,p6, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,2)),squeeze(response_avgblocksub(:,3,2))); %CS+/CS- vs CS-/CS- B2
t6 = ts.tstat; 

% % % % BLOCK 3
[h7,p7, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,3)),squeeze(response_avgblocksub(:,2,3))); %CS+/CS+ vs CS+/CS- B3
t7 = ts.tstat; 
[h8,p8, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,3)),squeeze(response_avgblocksub(:,3,3))); %CS+/CS+ vs CS-/CS- B3
t8 = ts.tstat; 
[h9,p9, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,3)),squeeze(response_avgblocksub(:,3,3))); %CS+/CS- vs CS-/CS- B3
t9 = ts.tstat; 
%%
data = squeeze(response_avgblocksub(:,2,3)); 

%data = squeeze(response_avgblocksub(:,1,:)); 
%data = data(:); 

% Perform the Lilliefors test
h = lillietest(data);
%h = jbtest(data);
%h = adtest(data);

histogram(data, 20)

if h == 0
    disp('Do not reject the null hypothesis: data appears to be normally distributed.');
else
    disp('Reject the null hypothesis: data does not appear to be normally distributed.');
end
set(gca, fontsize=16)
exportgraphics (gca, 'myP.png', 'Resolution', 300)

%% shapiroWilk

%data = squeeze(response_avgblocksub(:,3,3)); 

data = squeeze(response_avgblocksub(:,:,:)); 
data = data(:); 
histogram(data, 20)
[H, pValue, W] = swtest(data);

if pValue > 0.05
    disp('Do not reject the null hypothesis: data appears to be normally distributed.');
else
    disp('Reject the null hypothesis: data does not appear to be normally distributed.');
end







%% 
% % % % CS+/CS+ 
[h1,p1, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,1,2))); %B1 vs B2 CS+/CS+ 
t1 = ts.tstat; 
[h2,p2, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,1,3))); %B1 vs B3 CS+/CS+ 
t2 = ts.tstat; 
[h3,p3, ci, ts]=ttest(squeeze(response_avgblocksub(:,1,2)),squeeze(response_avgblocksub(:,1,3))); %B2 vs B3 CS+/CS+ 
t3 = ts.tstat; 

% % % % CS+/CS-
[h4,p4, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,1)),squeeze(response_avgblocksub(:,2,2))); %B1 vs B2 CS+/CS-
t4 = ts.tstat; 
[h5,p5, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,1)),squeeze(response_avgblocksub(:,2,3))); %B1 vs B3 CS+/CS- 
t5 = ts.tstat; 
[h6,p6, ci, ts]=ttest(squeeze(response_avgblocksub(:,2,2)),squeeze(response_avgblocksub(:,2,3))); %B2 vs B3 CS+/CS- 
t6 = ts.tstat; 

% % % % CS-/CS-
[h7,p7, ci, ts]=ttest(squeeze(response_avgblocksub(:,3,1)),squeeze(response_avgblocksub(:,3,2))); %B1 vs B2 CS-/CS- 
t7 = ts.tstat; 
[h8,p8, ci, ts]=ttest(squeeze(response_avgblocksub(:,3,1)),squeeze(response_avgblocksub(:,3,3))); %B1 vs B3 CS-/CS- 
t8 = ts.tstat; 
[h9,p9, ci, ts]=ttest(squeeze(response_avgblocksub(:,3,2)),squeeze(response_avgblocksub(:,3,3))); %B2 vs B3 CS-/CS-
t9 = ts.tstat; 




%% correct Bonferroni 18

p1 = p1*18; 
p2 = p2*18; 
p3 = p3*18; 
p4 = p4*9; 
p5 = p5*18; 
p6 = p6*18; 
p7 = p7*18; 
p8 = p8*18; 
p9 = p9*18; 



disp('corrected 18')


%% EFFECT SIZES
response_avgblocksub_here = response_avgblocksub; 
response_avgblocksub_here (sub2exc, :, :) = []; 
% CS+- EXT vs CS+- ACQU
x1 = squeeze(response_avgblocksub_here(:,2,1)); 
x2 = squeeze(response_avgblocksub_here(:,2,2)); 

% CS++ TEST vs CS+- TEST
%x1 = squeeze(response_avgblocksub_here(:,2,3)); 
%x2 = squeeze(response_avgblocksub_here(:,3,3)); 

[h4,p4, ci, ts]=ttest(x1,x2); %B1 vs B2 CS+/CS-
t4 = ts.tstat; 

%%
myD = meanEffectSize(x2, x1, Effect="cohen", Paired=true)


%% PROGRESSIVE LEARNING: FIRST LOAD DATA
clear
paths = load_paths_EXT;
path_trlinfo= paths.trlinfo;


% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
%            'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
%            'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
%             'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
%             'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
%             'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
%             'p_sub19', 'p_sub20', 'p_sub21'}'; % no p_sub08

allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11', 'c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27', 'c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07', 'p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21'}'; 


for subji=1:numel(allsubs)
    % plot responses for each type across experiment
    load(strcat(path_trlinfo,allsubs{subji},'_trlinfo.mat'))
    trialinfo=trlinfo;
    types={'cs+/cs+','cs+/cs-','cs-/cs-'};
    items=1:3;
    
    for i=1:3
        type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
        tmp=trialinfo(trialinfo(:,5)==i,7);
        %responses(i,:)=[tmp',nan(1,64-numel(tmp))]; % does not affect results and don't understand, removed
        responses(i,:)=tmp;
    end
    
    % average responses for each  n x type x response
    for ty=1:3
        %avg_response_type(subji,ty,:)=nanmean(responses(type_item==ty,:),1) 
        avg_response_type(subji,ty,:)=nanmean(responses(find(type_item==ty),:),1) 
    end
end



%% plot without smoothing
figure
hold on
fig_stuff=subplot(1,1,1);
cmap_default=fig_stuff.ColorOrder;
green= colormap(brewermap([1],'Greens'))
green = green*.9;
red = cmap_default(2,:);
yellow = cmap_default(3,:); 
newM(1, :) = red; 
newM(2, :) = yellow; 
newM(3, :) = green; 


for ty=1:3
    x1=1:size(avg_response_type,3)
    y1=squeeze(nanmean(avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(avg_response_type(:,ty,:),1))./sqrt(numel(allsubs));

    boundedline(x1, y1, b1, 'LineWidth', 2, 'cmap',newM(ty,:),'transparency',0.2,'alpha');
    %shadedErrorBar(x1, y1, b1, cmap_default(ty,:), .1)

    source_data_mean(ty, :) = y1; 
    source_data_std(ty,:) = b1; 

    %plot(squeeze(nanmean(avg_response_type(:,ty,:))))
    h=gca;
    h.YLim=[1 4];
    h.FontSize = 20;
end


plot([24 24],[1 4],':k', 'Linewidth', 2)
plot([48 48],[1 4],':k', 'Linewidth', 2)
%set(gca, 'xtick', [24 48], 'xticklabels', {'24' '48'})
set(gca, 'xtick', [12 36 56], 'xticklabels', {'ACQ' 'EXT' 'TEST'})

%%
writematrix(source_data_std, 'source_data_std.xlsx')

%% 
d2c =  squeeze(avg_response_type(:,1,:));

clear isNorMT; 
for triali = 1:64

    popTrial = d2c(:, triali); 
    [H, pValue, W] = swtest(popTrial);
    if pValue > 0.05
        disp('Do not reject the null hypothesis: data appears to be normally distributed.');
        isNorMT(triali) = 1; 
    else
        disp('Reject the null hypothesis: data does not appear to be normally distributed.');
        isNorMT(triali) = 0; 
    end
    

end


%% Plot with smoothing

window=1;

% smooth trajectories
move_avg=ones(1,window)/window;
for subji=1:numel(allsubs)
    for ty=1:3
     tmp_avg_response_type(subji,ty,:)=avg_response_type(subji,ty,:);
           
     % interpolate NaNs
     nan_ind=  find(isnan(tmp_avg_response_type(subji,ty,:)));
     
     while ~isempty(nan_ind)  
         for ind=1:numel(nan_ind)
             j=nan_ind(ind);
             if j>2 & j<(size(tmp_avg_response_type,3)-1)
                tmp_avg_response_type(subji,ty,j)=squeeze(nanmean(tmp_avg_response_type(subji,ty,j-1:j+1)));
             elseif j<=2
                tmp_avg_response_type(subji,ty,j)=squeeze(nanmean(tmp_avg_response_type(subji,ty,j+1)));
             elseif j>=(size(avg_response_type,3)-1)
                tmp_avg_response_type(subji,ty,j)=squeeze(nanmean(tmp_avg_response_type(subji,ty,j-1)));   
             end
         end
        nan_ind=  find(isnan(tmp_avg_response_type(subji,ty,:)));
     end
     
    filt_avg_response_type(subji,ty,:) = filtfilt(move_avg,1,squeeze(tmp_avg_response_type(subji,ty,:)));
    % plot average trajectories
    end
end


figure
hold on
fig_stuff=subplot(1,1,1);
cmap_default=fig_stuff.ColorOrder;
%cmap_default=cmap_default([1 2 4],:)


for ty=1:3
    x1=1:size(avg_response_type,3);
    y1=squeeze(nanmean(filt_avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(filt_avg_response_type(:,ty,:),1))./sqrt(numel(allsubs));
   
boundedline(x1, y1, b1, 'LineWidth', 2, 'cmap',cmap_default(ty,:),'transparency',0.2,'alpha');

%plot(squeeze(nanmean(avg_response_type(:,ty,:))))
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
h.YTickLabelRotation=90;
h.FontSize = 24;
end

plot([24 24],[1 4],':k', 'Linewidth', 2)
plot([48 48],[1 4],':k', 'Linewidth', 2)



%% plot all subjects without smoothing IMAGESC PLOTS

avg_response_type_NAN0 = avg_response_type; 
avg_response_type_NAN0(isnan(avg_response_type_NAN0))=0; 

figure(1); set(gcf, 'Position', [10 20 1850 1350])
tiledlayout(7, 8,'TileSpacing','compact', 'Padding','none');

for subji = 1:49
    nexttile
    imagesc(squeeze(avg_response_type_NAN0(subji, :, :))); colorbar
    title(allsubs{subji}, 'Interpreter','none')
end
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

%% plot all subjects with smoothing


figure(1); set(gcf, 'Position', [10 20 1850 1350])
tiledlayout(7, 8,'TileSpacing','compact', 'Padding','none');

for subji = 1:49
    nexttile
    imagesc(squeeze(filt_avg_response_type(subji, :, :))); 
    set(gca, 'clim', [0 4])
    colorbar
    title(allsubs{subji}, 'Interpreter','none')
end
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

%% plot new colors ; perform ANOVA and t-tests at every trial

filt_avg_response_type = avg_response_type; 

figure; set(gcf, 'Position', [50 50 500 400])
hold on
fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;

green= colormap(brewermap([1],'Greens'))
green = green*.9;
red = cmap_default(2,:);
yellow = cmap_default(3,:); 
newM(1, :) = red; 
newM(2, :) = yellow; 
newM(3, :) = green; 

for ty=1:3
    x1=1:size(filt_avg_response_type,3)
    y1=squeeze(nanmean(filt_avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(filt_avg_response_type(:,ty,:),1))./sqrt(numel(allsubs));
   
    boundedline(x1, y1, b1, 'LineWidth', 2, 'cmap',newM(ty,:),'transparency',0.2,'alpha');
    
    
    %plot(squeeze(nanmean(avg_response_type(:,ty,:))))
    h=gca;
    %h.YLim=[1,4];
    %h.YTick=[1;4];
    %h.YTickLabel={'Dangerous','Safe'};
    %h.YTickLabelRotation=90;
    h.FontSize = 18;
end


plot([24 24],[1 4],':k', 'Linewidth', 2)
plot([48 48],[1 4],':k', 'Linewidth', 2)
%set(gca, 'xtick', [24 48], 'xticklabels', {'24' '48'})
set(gca, 'xtick', [12 36 56], 'xticklabels', {'ACQ' 'EXT' 'TEST'})


% anova %here
clear tbl allP allF
for timei = 1:64
    d4ANOVA =  squeeze(filt_avg_response_type(:,:,timei)); 
    nSubj = size(d4ANOVA, 1); 
    d4ANOVA = d4ANOVA(:); 
    d4ANOVA(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
    d4ANOVA(:,3) = [1:nSubj 1:nSubj 1:nSubj];
    [p f] = RMAOV1(d4ANOVA);
    allP(timei) = p; 
    allF(timei) = f; 
end

hb = allP<0.05; 
clustinfo = bwconncomp(hb);
hb = double(hb); hb(hb == 0) = nan; hb(hb==1) = 1; 
clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(allF(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_ANOVA = allSTs(id);
else
    max_clust_obs_ANOVA = 0; 
end

% ttest 1
for i=1:size(filt_avg_response_type,3)
    [htest1(i),p1(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,1,i)),squeeze(filt_avg_response_type(:,2,i)));
end
t1 = [ts.tstat]; 

clustinfo1 = bwconncomp(htest1);
clear allSTs
for pxi = 1:length(clustinfo1.PixelIdxList)
   allSTs(pxi,:) = sum(t1(clustinfo1.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT1 = allSTs(id); 
else
    max_clust_obs_TT1 = 0; 
end
%remove non significant 
htest1 = zeros(1, length(htest1)); 
htest1(clustinfo1.PixelIdxList{id}) = 1; 
hb1 = double(htest1); hb1(hb1 == 0) = nan; hb1(hb1==1) = 1.15; 

[sxs iDp1ACQ] = min(p1(1:24)); 
t1ACQ = t1(iDp1ACQ); 

% ttest 2
for i=1:size(filt_avg_response_type,3)
    [htest2(i),p2(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,1,i)),squeeze(filt_avg_response_type(:,3,i)));
end
t2 = [ts.tstat]; 
clustinfo2 = bwconncomp(htest2);
clear allSTs
for pxi = 1:length(clustinfo2.PixelIdxList)
   allSTs(pxi,:) = sum(t2(clustinfo2.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT2 = allSTs(id); 
else
    max_clust_obs_TT2 = 0; 
end
%remove non significant 
htest2 = zeros(1, length(htest2)); 
htest2(clustinfo2.PixelIdxList{id}) = 1; 
hb2 = double(htest2); hb2(hb2 == 0) = nan; hb2(hb2==1) = 1.3; 

% ttest 3
for i=1:size(filt_avg_response_type,3)
    [htest3(i),p3(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,2,i)),squeeze(filt_avg_response_type(:,3,i)));
end
t3 = [ts.tstat]; 
clustinfo3 = bwconncomp(htest3);
clear allSTs
for pxi = 1:length(clustinfo3.PixelIdxList)
   allSTs(pxi,:) = sum(t3(clustinfo3.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT3 = allSTs(id); 
else
    max_clust_obs_TT3 = 0; 
end
%remove non significant 
htest3 = zeros(1, length(htest3)); 
htest3(clustinfo3.PixelIdxList{id}) = 1; 
hb3 = double(htest3); hb3(hb3 == 0) = nan; hb3(hb3==1) = 1.45; 


plot(hb, 'k', 'LineWidth',5)
plot(hb1, 'Color', '#EC5300', 'LineWidth',5); hold on %ORANGE
plot(hb2, 'Color', '#964B00', 'LineWidth',5); hold on %BROWN
plot(hb3,  'Color', [.5 .5 .5], 'LineWidth',5) % GREY


%filename = [paths.results.behavior 'average_per_type_filtered.png']
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

%%

allAb = max_clust_perm_TT3(abs(max_clust_perm_TT3) > abs(-2.56149440000433));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

%% plot new colors ; perform ANOVA and WILCOXON at every trial
clc
clearvars -except avg_response_type allsubs

filt_avg_response_type = avg_response_type; 

figure; set(gcf, 'Position', [50 50 500 400])
hold on
fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;

green= colormap(brewermap([1],'Greens'))
green = green*.9;
red = cmap_default(2,:);
yellow = cmap_default(3,:); 
newM(1, :) = red; 
newM(2, :) = yellow; 
newM(3, :) = green; 

for ty=1:3
    x1=1:size(filt_avg_response_type,3)
    y1=squeeze(nanmean(filt_avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(filt_avg_response_type(:,ty,:),1))./sqrt(numel(allsubs));
   
    boundedline(x1, y1, b1, 'LineWidth', 2, 'cmap',newM(ty,:),'transparency',0.2,'alpha');
    
    %plot(squeeze(nanmean(avg_response_type(:,ty,:))))
    h=gca;
    %h.YLim=[1,4];
    %h.YTick=[1;4];
    %h.YTickLabel={'Dangerous','Safe'};
    %h.YTickLabelRotation=90;
    h.FontSize = 18;
end


plot([24 24],[1 4],':k', 'Linewidth', 2)
plot([48 48],[1 4],':k', 'Linewidth', 2)
%set(gca, 'xtick', [24 48], 'xticklabels', {'24' '48'})
set(gca, 'xtick', [12 36 56], 'xticklabels', {'ACQ' 'EXT' 'TEST'})


% anova %here
clear tbl allP allF
for timei = 1:64
    d4ANOVA =  squeeze(filt_avg_response_type(:,:,timei)); 
    nSubj = size(d4ANOVA, 1); 
    d4ANOVA = d4ANOVA(:); 
    d4ANOVA(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
    d4ANOVA(:,3) = [1:nSubj 1:nSubj 1:nSubj];
    [p f] = RMAOV1(d4ANOVA);
    allP(timei) = p; 
    allF(timei) = f; 
end

hb = allP<0.05; 
clustinfo = bwconncomp(hb);
hb = double(hb); hb(hb == 0) = nan; hb(hb==1) = 1; 
clear allSTs
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(allF(clustinfo.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_ANOVA = allSTs(id);
else
    max_clust_obs_ANOVA = 0; 
end

% Wilcoxon 1
for i=1:size(filt_avg_response_type,3)
    %[htest1(i),p1(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,1,i)),squeeze(filt_avg_response_type(:,2,i)));
    differences = squeeze(filt_avg_response_type(:,1,i)) - squeeze(filt_avg_response_type(:,2,i)); 
    differences = differences(~isnan(differences));
    [p1(i) htest1(i) stats(i)] = signrank(differences);
end
Z1 = [stats.zval];

clustinfo1 = bwconncomp(htest1);
clear allSTs
for pxi = 1:length(clustinfo1.PixelIdxList)
   %allSTs(pxi,:) = length(clustinfo1.PixelIdxList{pxi});% 
   allSTs(pxi,:) = sum(Z1(clustinfo1.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT1 = allSTs(id); 
else
    max_clust_obs_TT1 = 0; 
end
%remove non significant 
htest1 = zeros(1, length(htest1)); 
htest1(clustinfo1.PixelIdxList{id}) = 1; 
hb1 = double(htest1); hb1(hb1 == 0) = nan; hb1(hb1==1) = 1.15; 


% Wilcoxon 2
clear differences
for i=1:size(filt_avg_response_type,3)
    %[htest2(i),p2(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,1,i)),squeeze(filt_avg_response_type(:,3,i)));
    differences = squeeze(filt_avg_response_type(:,1,i)) - squeeze(filt_avg_response_type(:,3,i)); 
    differences = differences(~isnan(differences));
    [p2_prev htest2_prev stats2_prev] = signrank(differences); 
    p2(i) = p2_prev; 
    htest2(i) = htest2_prev;
    if isfield(stats2_prev, 'zval')
        stats2(i) = stats2_prev;
    else
        %stats2(i) = nan;
    end
end

Z2 = [stats2.zval];

clustinfo2 = bwconncomp(htest2);
clear allSTs
for pxi = 1:length(clustinfo2.PixelIdxList)
   %allSTs(pxi,:) = length(clustinfo2.PixelIdxList{pxi});% 
   allSTs(pxi,:) = sum(Z2(clustinfo2.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT2 = allSTs(id); 
else
    max_clust_obs_TT2 = 0; 
end
%remove non significant 
htest2 = zeros(1, length(htest2)); 
htest2(clustinfo2.PixelIdxList{id}) = 1; 
hb2 = double(htest2); hb2(hb2 == 0) = nan; hb2(hb2==1) = 1.3; 

% Wilcoxon 3
clear differences
for i=1:size(filt_avg_response_type,3)
    %[htest3(i),p3(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type(:,2,i)),squeeze(filt_avg_response_type(:,3,i)));
    differences = squeeze(filt_avg_response_type(:,2,i)) - squeeze(filt_avg_response_type(:,3,i)); 
    differences = differences(~isnan(differences));
    [p3_prev htest3_prev stats3_prev] = signrank(differences); 
    p3(i) = p3_prev; 
    htest3(i) = htest3_prev;
    if isfield(stats3_prev, 'zval')
        stats3(i) = stats3_prev;
    else
        %stats2(i) = nan;
    end
end
Z3 = [stats3.zval];

clustinfo3 = bwconncomp(htest3);
clear allSTs
for pxi = 1:length(clustinfo3.PixelIdxList)
   %allSTs(pxi,:) = length(clustinfo3.PixelIdxList{pxi});% 
   allSTs(pxi,:) = sum(Z3(clustinfo3.PixelIdxList{pxi}));% 
end
if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    max_clust_obs_TT3 = allSTs(id); 
else
    max_clust_obs_TT3 = 0; 
end
%remove non significant 
htest3 = zeros(1, length(htest3)); 
htest3(clustinfo3.PixelIdxList{id}) = 1; 
hb3 = double(htest3); hb3(hb3 == 0) = nan; hb3(hb3==1) = 1.45; 


plot(hb, 'k', 'LineWidth',5)
plot(hb1, 'Color', '#EC5300', 'LineWidth',5); hold on %ORANGE
plot(hb2, 'Color', '#964B00', 'LineWidth',5); hold on %BROWN
plot(hb3,  'Color', [.5 .5 .5], 'LineWidth',5) % GREY


%filename = [paths.results.behavior 'average_per_type_filtered.png']
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)



%% permutations ANOVA

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    for subji = 1:size(filt_avg_response_type)
        filt_avg_response_type_PERM(subji, :, :) = filt_avg_response_type(subji, randperm(3), :); 
    end
       
    
    % anova %here
    clear tbl, clc
    for timei = 1:64
        d4ANOVA =  squeeze(filt_avg_response_type_PERM(:,:,timei)); 
        nSubj = size(d4ANOVA, 1); 
        d4ANOVA = d4ANOVA(:); 
        d4ANOVA(:,2) = [ones(1,nSubj) ones(1,nSubj)*2 ones(1,nSubj)*3];
        d4ANOVA(:,3) = [1:nSubj 1:nSubj 1:nSubj];
        [p f] = RMAOV1(d4ANOVA);
        allP(timei) = p; 
        allF(timei) = f; 
    end
    
    hPerm = allP < 0.05; tPerm = allF; 

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = max(abs((allSTs)));
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

clear p mcsP
mcsP = max_clust_sum_perm;
allAb = mcsP(abs(mcsP) > abs(max_clust_obs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% permutations TTESTS

nPerm = 1000; 
clear max_clust_sum_perm max_clust_perm_TT1 max_clust_perm_TT2 max_clust_perm_TT3
for permi = 1:nPerm
    progress_in_console(permi)

    for subji = 1:size(filt_avg_response_type)
        filt_avg_response_type_PERM(subji, :, :) = filt_avg_response_type(subji, randperm(3), :); 
    end
       

    % ttest 1
    for i=1:size(filt_avg_response_type_PERM,3)
        [htest1(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,1,i)),squeeze(filt_avg_response_type_PERM(:,2,i)));
    end
    t1 = [ts.tstat]; 
    hb1 = double(htest1); hb1(hb1 == 0) = nan; hb1(hb1==1) = 1.2; 
    clustinfo = bwconncomp(htest1);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t1(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT1(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT1(permi, :) = 0; 
    end
    
    % ttest 2
    for i=1:size(filt_avg_response_type_PERM,3)
        [htest2(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,1,i)),squeeze(filt_avg_response_type_PERM(:,3,i)));
    end
    t2 = [ts.tstat]; 
    hb2 = double(htest2); hb2(hb2 == 0) = nan; hb2(hb2==1) = 1.3; 
    clustinfo = bwconncomp(htest2);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t2(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT2(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT2(permi, :) = 0; 
    end
    
    % ttest 3
    for i=1:size(filt_avg_response_type_PERM,3)
        [htest3(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,2,i)),squeeze(filt_avg_response_type_PERM(:,3,i)));
    end
    t3 = [ts.tstat]; 
    hb3 = double(htest3); hb3(hb3 == 0) = nan; hb3(hb3==1) = 1.4; 
    clustinfo = bwconncomp(htest3);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(t3(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT3(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT3(permi, :) = 0; 
    end


    

end

clear p mcsP
allAb1 = max_clust_perm_TT1(abs(max_clust_perm_TT1) > abs(max_clust_obs_TT1));
p1 = 1 - ((nPerm-1) - (length (allAb1)))  / nPerm; 

%allAb1ACQ = max_clust_perm_TT1(abs(max_clust_perm_TT1) > abs(t1ACQ));
%p1ACQ = 1 - ((nPerm-1) - (length (allAb1ACQ)))  / nPerm; 

allAb2 = max_clust_perm_TT2(abs(max_clust_perm_TT2) > abs(max_clust_obs_TT2));
p2 = 1 - ((nPerm-1) - (length (allAb2)))  / nPerm; 

allAb3 = max_clust_perm_TT3(abs(max_clust_perm_TT3) > abs(max_clust_obs_TT3));
p3 = 1 - ((nPerm-1) - (length (allAb3)))  / nPerm; 



%% permutations WILCOXON
clc
nPerm = 1000; 
clear max_clust_sum_perm max_clust_perm_TT1 max_clust_perm_TT2 max_clust_perm_TT3
for permi = 1:nPerm
    progress_in_console(permi)

    for subji = 1:size(filt_avg_response_type)
        filt_avg_response_type_PERM(subji, :, :) = filt_avg_response_type(subji, randperm(3), :); 
    end
       
    % ttest 1
    for i=1:size(filt_avg_response_type_PERM,3)
        %[htest1(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,1,i)),squeeze(filt_avg_response_type_PERM(:,2,i)));
        differences = filt_avg_response_type_PERM(:,1,i) - filt_avg_response_type_PERM(:,2,i); 
        [p1_prev htest1_prev stats1_prev] = signrank(differences); 
        p1(i) = p1_prev; 
        htest1(i) = htest1_prev;
        if isfield(stats1_prev, 'zval')
            stats1(i) = stats1_prev;
        else
            %stats2(i) = nan;
        end
    end
    
    Z1 = [stats1.zval];
    hb1 = double(htest1); hb1(hb1 == 0) = nan; hb1(hb1==1) = 1.2; 
    clustinfo = bwconncomp(htest1);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(Z1(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT1(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT1(permi, :) = 0; 
    end
    
    % ttest 2
    for i=1:size(filt_avg_response_type_PERM,3)
        %[htest2(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,1,i)),squeeze(filt_avg_response_type_PERM(:,3,i)));
        differences = filt_avg_response_type_PERM(:,1,i) - filt_avg_response_type_PERM(:,3,i); 
        [p2_prev htest2_prev stats2_prev] = signrank(differences); 
        p2(i) = p2_prev; 
        htest2(i) = htest2_prev;
        if isfield(stats2_prev, 'zval')
            stats2(i) = stats2_prev;
        else
            %stats2(i) = nan;
        end
    end
    Z2 = [stats2.zval];
    hb2 = double(htest2); hb2(hb2 == 0) = nan; hb2(hb2==1) = 1.3; 
    clustinfo = bwconncomp(htest2);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(Z2(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT2(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT2(permi, :) = 0; 
    end
    
    % ttest 3
    for i=1:size(filt_avg_response_type_PERM,3)
        %[htest3(i),p(i) ci ts(i)]=ttest(squeeze(filt_avg_response_type_PERM(:,2,i)),squeeze(filt_avg_response_type_PERM(:,3,i)));
        differences = filt_avg_response_type_PERM(:,2,i) - filt_avg_response_type_PERM(:,3,i); 
        [p3_prev htest3_prev stats3_prev] = signrank(differences); 
        p3(i) = p3_prev; 
        htest3(i) = htest3_prev;
        if isfield(stats3_prev, 'zval')
            stats3(i) = stats3_prev;
        else
            %stats2(i) = nan;
        end
    end
    Z3 = [stats3.zval];
    hb3 = double(htest3); hb3(hb3 == 0) = nan; hb3(hb3==1) = 1.4; 
    clustinfo = bwconncomp(htest3);
    clear allSTs
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi,:) = sum(Z3(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        max_clust_perm_TT3(permi, :) = allSTs(id); 
    else
        max_clust_perm_TT3(permi, :) = 0; 
    end


    

end

clear p mcsP
allAb1 = max_clust_perm_TT1(abs(max_clust_perm_TT1) > abs(max_clust_obs_TT1));
p1 = 1 - ((nPerm-1) - (length (allAb1)))  / nPerm; 

%allAb1ACQ = max_clust_perm_TT1(abs(max_clust_perm_TT1) > abs(t1ACQ));
%p1ACQ = 1 - ((nPerm-1) - (length (allAb1ACQ)))  / nPerm; 

allAb2 = max_clust_perm_TT2(abs(max_clust_perm_TT2) > abs(max_clust_obs_TT2));
p2 = 1 - ((nPerm-1) - (length (allAb2)))  / nPerm; 

allAb3 = max_clust_perm_TT3(abs(max_clust_perm_TT3) > abs(max_clust_obs_TT3));
p3 = 1 - ((nPerm-1) - (length (allAb3)))  / nPerm; 




%% LME
clc
clear d4LME
d4LME = response_avgblocksub(:);

subj4LME = [1:50]; subj4LME = [subj4LME subj4LME subj4LME subj4LME subj4LME subj4LME subj4LME subj4LME subj4LME ]'; 
type = [[ones(1, 50) ones(1, 50)*2 ones(1, 50)*3]']; type = [type; type; type]
block = [[ones(1, 150) ones(1, 150)*2 ones(1, 150)*3]']; 

d4LME = [d4LME subj4LME block type]
tbl = table(d4LME(:,1), d4LME(:,2), d4LME(:,3), d4LME(:,4), ...
    'VariableNames',{'performance','subID', 'block','type'});

lme = fitlme(tbl,'performance ~ block + type + block*type + (1|subID)'); % random intercept model

lme




%% WITHIN SUBJECTS DESING IN MATLAB
clear d2ADD
d2ADD(:, 1) = squeeze(response_avgblocksub(:, 1, 1));
d2ADD(:, 2) = squeeze(response_avgblocksub(:, 1, 2));
d2ADD(:, 3) = squeeze(response_avgblocksub(:, 1, 3));

d4ANOVA = [[1:50]' ]; 
% organize the data in a table
T = array2table(d4ANOVA(:,2:end));
T.Properties.VariableNames = {'nnHBLNET' 'nnHALEX' 'nnHCAT'};
% create the within-subjects design
withinDesign = table([1 2 3]','VariableNames',{'Model'});
withinDesign.Model = categorical(withinDesign.Model);
% create the repeated measures model and do the anova
rm = fitrm(T,'nnHBLNET-nnHCAT ~ 1','WithinDesign',withinDesign);
AT = ranova(rm,'WithinModel','Model'); % remove comma to see ranova's table
%tbl = multcompare(rm, 'Model', 'ComparisonType', 'tukey-kramer'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'
tbl = multcompare(rm, 'Model', 'ComparisonType', 'bonferroni'); % see: help RepeatedMeasuresModel/multcompare;  'tukey-kramer' (default), 'dunn-sidak', 'bonferroni','scheffe'


% output a conventional anova table
disp(anovaTable(AT, 'Measure (units)'));


%%
d4ANOVA = [response_avgblocksub]; 
d4ANOVA(:,2) = [ones(1,50) ones(1,50)*2 ones(1,50)*3];
d4ANOVA(:,3) = [1:15 1:15 1:15];
[p f] = RMAOV1(d4ANOVA);
allP(freqi, timei) = p; 
allF(freqi, timei) = f; 


%%
d4anova = squeeze(response_avgblock(1,:,:));
[p, tbl, stats] = anova1(d4anova)

multcompare(stats)

%% split performance block 3

block_def=[1 24;25 48;49 64];
 b=3
 response_avgsumblock3=nanmean(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1);
 response_stdsumblock3=nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1)./sqrt(numel(allsubs));

all_response_avgblock3=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);

figure
h= subplot(1,3,1)

ba = bar(response_avgsumblock3); hold on;
ba.FaceColor = 'flat';
ba.CData(1,:) = lCs{1};
ba.CData(2,:) = lCs{2};
ba.CData(3,:) = lCs{3};

errorbar(1:3,response_avgsumblock3,response_stdsumblock3, 'k', LineWidth=2)
scatter(reshape(repmat((1:3),numel(allsubs),1),numel(allsubs)*3,1), ...
        reshape(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),numel(allsubs)*3,1), 10,'ko')
%plot(nansum(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3)')


set(h,'xtick', 1:3, 'xticklabel', {'cs+/cs+','cs+/cs-','cs-/cs-'}, 'FontSize', 14)

avg_response_typeC=(avg_response_type(:,:,block_def(b,1):block_def(b,2)));


%filename = [paths.results.behavior 'perf_at_test_per_type.png']
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

