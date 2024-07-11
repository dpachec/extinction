%% Extinction behavioral analysis
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

%later added
%15. number of item rep in each block
%16 number of us in total
%path_data='D:\Extinction\iEEG\rawdata\extinction_ieeg\';
%path_out='D:\Extinction\iEEG\data\preproc\trialinfo\';


clear 

paths = load_paths_EXT;


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21'}'; % no p_sub08


responses = [];


for subji=1:numel(allsubs)

    sub=allsubs{subji};
    info_file=strcat(paths.trlinfo,sub,'_trlinfo');
    load(info_file)  
    trialinfo = trlinfo; 
    allTrialinfo{subji,:} = trlinfo;

    label={'Dangerous','Safe'};
    
    
    types={'cs+/cs+','cs+/cs-','cs-/cs-'};
    items=1:3;


    for i=1:3
        type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
        tmp=trialinfo(trialinfo(:,5)==i,7);
        responses(i,:)=[tmp',nan(1,length(responses)-numel(tmp))];
        tmp_us=trialinfo(trialinfo(:,5)==i,9);
        us_vec(i,:)=[tmp_us',nan(1,length(tmp_us)-numel(tmp_us))]
    end
    % 
    % figure
    % for i=1:3
    %     subplot(2,3,i)
    %     stem(responses(i,:))
    %     title(types{type_item(i)})
    % end
    
    % % % uncomnment to plot one fig per subject
% % % % %     figure
% % % % %     for ty=1:3
% % % % %         avg_response_type(ty,:)=nanmean(responses(type_item==ty,:),1) 
% % % % %         avg_us(ty,:)=us_vec(type_item==ty,:);
% % % % %         subplot(2,2,ty)
% % % % %         hold on
% % % % %         rectangle('Position',[0 0 24 4],'FaceColor',[1 .9 .9])
% % % % %         rectangle('Position',[24 0 24 4],'FaceColor',[0.9 .9 .9])
% % % % %         rectangle('Position',[48 0 16 4],'FaceColor',[1 .8 1])
% % % % %         rectangle('Position',[64 0 16 4],'FaceColor',[1 1 .8])
% % % % %         
% % % % %         
% % % % %         stem(avg_response_type(ty,:))
% % % % %         title(types{ty})
% % % % %         h=gca;
% % % % %         h.YTick=[0.1;4-(0.2*4)];
% % % % %         h.YTickLabel=label;
% % % % %         h.YTickLabelRotation=90;
% % % % %     
% % % % %     end
% % % % %     
% % % % %     mkdir(paths.results.behavior)
% % % % %     filename = [paths.results.behavior sub '_behav.png']
% % % % %     exportgraphics(gca, filename, 'Resolution',300)
    %close all

end
clear avg_response_type tmp tmp_us

%% plot separately for each type

clear 

type2u = 3; % (1, 2, 3) ; Cs+Cs+, Cs+Cs-, Cs-Cs-

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
  

%% analyse behavior: average responses
clear
paths = load_paths_EXT;
path_trlinfo= paths.trlinfo;


allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' 'p_sub01','p_sub02', ...
            'p_sub03','p_sub04','p_sub05','p_sub06','p_sub07','p_sub09', 'p_sub10', ...
            'p_sub11','p_sub12','p_sub13','p_sub14','p_sub15', 'p_sub16','p_sub17','p_sub18', ... 
            'p_sub19', 'p_sub20', 'p_sub21'}'; % no p_sub08


for subji=1:numel(allsubs)
% plot responses for each type across experiment
load(strcat(path_trlinfo,allsubs{subji},'_trlinfo.mat'))
trialinfo=trlinfo;
types={'cs+/cs+','cs+/cs-','cs-/cs-'};
items=1:3;




for i=1:3
type_item(i)=mean(trialinfo(trialinfo(:,5)==i,6));
tmp=trialinfo(trialinfo(:,5)==i,7);
responses(i,:)=[tmp',nan(1,64-numel(tmp))];
end

% average responses for each  n x type x response
for ty=1:3
  avg_response_type(subji,ty,:)=nanmean(responses(type_item==ty,:),1) 
    
end
end

figure
hold on
fig_stuff=subplot(1,1,1)
cmap_default=fig_stuff.ColorOrder;


for ty=1:3
    x1=1:size(avg_response_type,3)
    y1=squeeze(nanmean(avg_response_type(:,ty,:)));
    b1=squeeze(nanstd(avg_response_type(:,ty,:),1))./sqrt(numel(allsubs));
   
    boundedline(x1, y1, b1, 'LineWidth', 2, 'cmap',cmap_default(ty,:),'transparency',0.2,'alpha');
    %shadedErrorBar(x1, y1, b1, cmap_default(ty,:), .1)
    
    %plot(squeeze(nanmean(avg_response_type(:,ty,:))))
    h=gca;
    h.YLim=[1,4];
    h.YTick=[1;4];
    h.YTickLabel={'Dangerous','Safe'};
    h.YTickLabelRotation=90;
    h.FontSize = 18;
end


plot([24 24],[1 4],':k', 'Linewidth', 2)
plot([48 48],[1 4],':k', 'Linewidth', 2)

%filename = [paths.results.behavior 'average_per_type.png']

filename = ['myP.png']
exportgraphics(gcf, filename, 'Resolution',300)


%% Plot with smoothing



window=2;
% smooth trajectories
move_avg=ones(1,window)/window;
for vp=1:numel(allsubs)
    for ty=1:3
     tmp_avg_response_type(vp,ty,:)=avg_response_type(vp,ty,:);
           
     % interpolate NaNs
     nan_ind=  find(isnan(tmp_avg_response_type(vp,ty,:)));
     
     while ~isempty(nan_ind)  
         for ind=1:numel(nan_ind)
             j=nan_ind(ind);
             if j>2 & j<(size(tmp_avg_response_type,3)-1)
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j-1:j+1)));
             elseif j<=2
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j+1)));
             elseif j>=(size(avg_response_type,3)-1)
             tmp_avg_response_type(vp,ty,j)=squeeze(nanmean(tmp_avg_response_type(vp,ty,j-1)));   
             end
         end
        nan_ind=  find(isnan(tmp_avg_response_type(vp,ty,:)));
     end
     
filt_avg_response_type(vp,ty,:) = filtfilt(move_avg,1,squeeze(tmp_avg_response_type(vp,ty,:)));
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




%%
for t=1:size(filt_avg_response_type,3)
[htest(t),p(t)]=ttest(squeeze(filt_avg_response_type(:,1,t)),squeeze(filt_avg_response_type(:,2,t)))
end

figure
for ty=1:3
subplot(2,2,ty)
hold on
rectangle('Position',[0 0 24 4],'FaceColor',[1 .9 .9])
rectangle('Position',[24 0 24 4],'FaceColor',[0.9 .9 .9])
rectangle('Position',[48 0 16 4],'FaceColor',[1 .8 1])


stem(squeeze(nanmean(filt_avg_response_type(:,ty,:))))
title(types{ty})
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
h.YTickLabelRotation=90;


end

figure

contrasts=[1,2;1,3;2,3];
contrast_label={'cs+/cs+ vs cs+/cs-','cs+/cs+ vs cs-/cs-','cs+/cs- vs cs-/cs-'};
for con=1:3
    con1=contrasts(con,1);
    con2=contrasts(con,2);

subplot(2,3,con)
hold on
rectangle('Position',[0 -2 24 4],'FaceColor',[1 .9 .9])
rectangle('Position',[24 -2 24 4],'FaceColor',[0.9 .9 .9])
rectangle('Position',[48 -2 16 4],'FaceColor',[1 .8 1])
[htest,p,~,tstat]=ttest(squeeze(filt_avg_response_type(:,con1,:)),squeeze(filt_avg_response_type(:,con2,:)))


h=gca;
stem(squeeze(nanmean(filt_avg_response_type(:,con1,:)-filt_avg_response_type(:,con2,:))))
hold on
plot(htest)
title(types{ty})
h=gca;
h.YLim=[-2,2];
h.YTick=[-2;2];
h.YTickLabel={'less safe','safer'};
title(contrast_label{con})


end

%% plot only new colors

figure
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
    x1=1:size(avg_response_type,3)
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


% anova 
clear tbl, clc
for timei = 1:64
    d4anova = squeeze(filt_avg_response_type(:,:,timei));    
    [p(timei), tbl] = anova1(d4anova,[],'off');
    t(timei) = tbl{2,5};
end
hb = p<0.05; 
clustinfo = bwconncomp(hb);
hb = double(hb); hb(hb == 0) = nan; hb(hb==1) = 1.1; 
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi,:) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end
[max2u id] = max(abs(allSTs));
max_clust_obs = allSTs(id)

% ttests
for t=1:size(filt_avg_response_type,3)
[htest(t),p(t)]=ttest(squeeze(filt_avg_response_type(:,1,t)),squeeze(filt_avg_response_type(:,2,t)))
end
hb1 = htest; hb1(hb1 == 0) = nan; hb1(hb1==1) = 1.4; 
% ttests
for t=1:size(filt_avg_response_type,3)
[htest(t),p(t)]=ttest(squeeze(filt_avg_response_type(:,1,t)),squeeze(filt_avg_response_type(:,3,t)))
end
hb2 = htest; hb2(hb2 == 0) = nan; hb2(hb2==1) = 1.3; 
% ttests
for t=1:size(filt_avg_response_type,3)
[htest(t),p(t)]=ttest(squeeze(filt_avg_response_type(:,2,t)),squeeze(filt_avg_response_type(:,3,t)))
end
hb3 = htest; hb3(hb3 == 0) = nan; hb3(hb3==1) = 1.5; 

plot(hb, 'k', 'LineWidth',5)
plot(hb1, 'Color', '#FFA500', 'LineWidth',5)
plot(hb2, 'Color', '#964B00', 'LineWidth',5)
plot(hb3,  'Color', [.5 .5 .5], 'LineWidth',5)


%filename = [paths.results.behavior 'average_per_type_filtered.png']
filename = 'myP.png'; 
exportgraphics(gcf, filename, 'Resolution',300)

%% permutations

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    for subji = 1:50
        filt_avg_response_type_PERM(subji, :, :) = filt_avg_response_type(subji, randperm(3), :); 
    end
    % anova 
    clear tbl p
    for timei = 1:64
        d4anova = squeeze(filt_avg_response_type_PERM(:,:,timei));    
        [p(timei), tbl] = anova1(d4anova,[],'off');
        t(timei) = tbl{2,5};
    end
    hPerm = p < 0.05; tPerm = t; 

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



%% average rating each block

block_def=[1 24;25 48;49 64];

lCs = {red; yellow; green};
figure
hold  on
for b=1:3
    response_avgblock(:,b)=nanmean(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),1),3);
   response_stdblock(:,b)=(nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1))./sqrt(numel(allsubs));
    response_avgblocksub(:,:,b)=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);

end
for ty=1:3
  errorbar(1:3,response_avgblock(ty,:),response_stdblock(ty,:), 'Color', lCs{ty}, 'LineWidth',2)
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

%% statistics average in each block 
clear response_avgblock
for b=1:3
    response_avgblock(b, :, :)=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);
   %response_stdblock(:,b)=(nanstd(nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3),1))./sqrt(numel(allsubs));
   % response_avgblocksub(:,:,b)=nanmean(avg_response_type(:,:,block_def(b,1):block_def(b,2)),3);
end

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

%% split performance block 4 in A & B



path_trlinfo=paths.trlinfo;

%selvps={'sub01','sub02','sub03','sub04','sub05','sub06','sub07','sub08','sub09','sub10','sub11','sub12'};


for vp=1:numel(allsubs)
% plot responses for each type across experiment

load(strcat(path_trlinfo,allsubs{vp},'_trlinfo.mat'))
trialinfo=trlinfo;




trialinfoA=trialinfo(trialinfo(:,2)==4&trialinfo(:,4)==1,:);
trialinfoB=trialinfo(trialinfo(:,2)==4&trialinfo(:,4)==2,:);
 

for i=1:3
type_itemA(i)=nanmean(trialinfoA(trialinfoA(:,5)==i,6));
responsesA(i,:)=trialinfoA(trialinfoA(:,5)==i,7);

type_itemB(i)=nanmean(trialinfoB(trialinfoB(:,5)==i,6));
responsesB(i,:)=trialinfoB(trialinfoB(:,5)==i,7);
end

% average responses for each  n x type x response

for ty=1:3
  avg_response_typeB(vp,ty,:)=nanmean(responsesB(type_itemB==ty,:),1) 
  avg_response_typeA(vp,ty,:)=nanmean(responsesA(type_itemA==ty,:),1) 

end

end

all_response_avgblock4A=nanmean(avg_response_typeA,3);
all_response_avgblock4B=nanmean(avg_response_typeB,3);

 response_avgsumblock4A=nanmean(nanmean(avg_response_typeA(:,:,:),3),1);
 response_stdsumblock4A=nanstd(nanmean(avg_response_typeA(:,:,:),3),1)./sqrt(numel(allsubs));

 response_avgsumblock4B=nanmean(nanmean(avg_response_typeB(:,:,:),3),1);
 response_stdsumblock4B=nanstd(nanmean(avg_response_typeB(:,:,:),3),1)./sqrt(numel(allsubs));


subplot(1,3,2)
bar(response_avgsumblock4A)
hold on
errorbar(1:3,response_avgsumblock4A,response_stdsumblock4A)
% individual color and size of scatter dots
%sz = repmat(1:numel(selvps),1,4).*10;
%c = repmat(linspace(1,10,numel(selvps)),1,4);

scatter(reshape(repmat((1:3),numel(allsubs),1),numel(allsubs)*3,1),reshape(nanmean(avg_response_typeA(:,:,:),3),numel(allsubs)*3,1))
%plot(nansum(avg_response_typeA(:,:,:),3)')

h=gca;
h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
title('ABA')

subplot(1,3,3)
bar(response_avgsumblock4B)
hold on
scatter(reshape(repmat((1:3),numel(allsubs),1),numel(allsubs)*3,1),reshape(nanmean(avg_response_typeB(:,:,:),3),numel(allsubs)*3,1))
errorbar(1:3,response_avgsumblock4B,response_stdsumblock4B)
%plot(nansum(avg_response_typeB(:,:,:),3)')

h=gca;
h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
title('ABB')


figure
hold on
errorbar(1:3,response_avgsumblock4A,response_stdsumblock4B)
errorbar(1:3,response_avgsumblock4B,response_stdsumblock4B)
errorbar(1:3,response_avgsumblock3,response_stdsumblock4B)
h=gca;
h.XTick=[1:3];
h.XLim=[0.5,3.5];

h.XTickLabel={'cs+/cs+','cs+/cs-','cs-/cs-'}
legend('ABA','ABB','ABC')
ylabel({'average response';'safe                                                  dangerous'})

%% plot time course test
%addpath(genpath('E:\matlab_tools\boundedline\kakearney-boundedline-pkg-8179f9a\boundedline'))


figure
fig_stuff=subplot(1,3,1)
cmap_default=fig_stuff.ColorOrder;
for ty=1:3
    x1=1:size(avg_response_typeA,3)
    y1=squeeze(nanmean(avg_response_typeA(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeA(:,ty,:),1))./sqrt(numel(allsubs));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end
title('Test in A')

h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};

subplot(1,3,2)
for ty=1:3
    x1=1:size(avg_response_typeB,3)
    y1=squeeze(nanmean(avg_response_typeB(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeB(:,ty,:),1))./sqrt(numel(allsubs));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end


h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
title('Test in B')

h.YTickLabel={'Dangerous','Safe'};
subplot(1,3,3)

for ty=1:3
    x1=1:size(avg_response_typeC,3)
    y1=squeeze(nanmean(avg_response_typeC(:,ty,:),1));
    b1=squeeze(nanstd(avg_response_typeC(:,ty,:),1))./sqrt(numel(allsubs));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(ty,:),'transparency',0.1,'alpha');
end
h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
title('Test in C')
h.YTickLabel={'Dangerous','Safe'};
legend('cs+/cs+','','cs+/cs-','','cs-/cs-','');

%%

figure
fig_stuff=subplot(1,4,1)
cmap_default=fig_stuff.ColorOrder;

contexts={'A','B','C'};
conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
    figure
for ty=1:3
fig_stuff=subplot(1,3,ty)

cmap_default=fig_stuff.ColorOrder;
    
for con=1:3
    
    eval(strcat('tmp=avg_response_type',contexts{con}),';')
    x1=1:size(tmp,3)
    y1=squeeze(nanmean(tmp(:,ty,:),1));
    b1=squeeze(nanstd(tmp(:,ty,:),1))./sqrt(numel(selvps)); 
boundedline(x1, y1, b1, 'cmap',cmap_default(con,:),'transparency',0.1,'alpha');

end
title(conds{ty})

h=gca;
h.YLim=[1,4];
h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
end

legend('Test in A','','Test in B','','Test in C','');

%% plot distributions


% difference A vs b in each subject

    figure

diffab=all_response_avgblock4A-all_response_avgblock4B;

conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
    subplot(1,3,1)
    hist(diffab)
    legend (conds)
    title('difference A vs B')
    
    
diffac=all_response_avgblock4A-all_response_avgblock3;
    subplot(1,3,2)
    
hist(diffac)
    legend (conds)
    title('difference A vs C')
    
    diffbc=all_response_avgblock4B-all_response_avgblock3;
    subplot(1,3,3)
    
hist(diffbc)
    legend (conds)
    title('difference B vs C')

%%

    figure

diffab=all_response_avgblock4A-all_response_avgblock4B;
diffac=all_response_avgblock4A-all_response_avgblock3;
diffbc=all_response_avgblock4B-all_response_avgblock3;




conds={'cs+/cs+','cs+/cs-','cs-/cs-'};
diff_cond={'difference A vs B','difference A vs C','difference B vs C'};
bin_vec=-3:0.5:3;
    subplot(1,3,1)
    hist([diffab(:,1),diffac(:,1),diffbc(:,1)],bin_vec)


    legend (diff_cond)
    title(conds{1})
    h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    
    subplot(1,3,2)
    hist([diffab(:,2),diffac(:,2),diffbc(:,2)],bin_vec)

    legend (diff_cond)
    title(conds{2})
        h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    subplot(1,3,3)
    
    hist([diffab(:,3),diffac(:,3),diffbc(:,3)],bin_vec)

    legend (diff_cond)
    title(conds{3})
    h=gca;
    h.YLim=[0 20];
    h.XLim=[-3 3];
    
    
    %% bounded lines for learning -extinction - test
    mean_response_type=squeeze(nanmean(filt_avg_response_type(:,:,1:block_def(3,2))));
    std_response_type=squeeze(nanstd(filt_avg_response_type(:,:,1:block_def(3,2))))./sqrt(numel(selvps));
    
    figure
   hold on
   
rectangle('Position',[0 0 24 4])
rectangle('Position',[24 0 24 4])
rectangle('Position',[48 0 16 4])

for con=1:3
    
    x1=1:size(mean_response_type,2);
    y1=mean_response_type(con,:);
    b1=std_response_type(con,:);
boundedline(x1, y1, b1, 'cmap',cmap_default(con,:),'transparency',0.1,'alpha');

end
    h=gca;
    h.XLim=[1 block_def(3,2)];
    h.YLim=[1 4];
    h.YTick=[1;4];
h.YTickLabel={'Dangerous','Safe'};
xlabel('Trial No')
legend('cs+/cs+','','cs+/cs-','','cs-/cs-','');    

%%
% learning
conds={'plusplus','plusminus','minusminus'};
export_spss=[squeeze(response_avgblocksub(:,1,1)),squeeze(response_avgblocksub(:,2,1)),squeeze(response_avgblocksub(:,3,1)),...
    squeeze(response_avgblocksub(:,1,2)),squeeze(response_avgblocksub(:,2,2)),squeeze(response_avgblocksub(:,3,2)),...
     squeeze(response_avgblocksub(:,1,3)),squeeze(response_avgblocksub(:,2,3)),squeeze(response_avgblocksub(:,3,3))];
names={strcat('learn1_',conds{1}),strcat('learn1_',conds{2}),strcat('learn1_',conds{3}),...
    strcat('learn2_',conds{1}),strcat('learn2_',conds{2}),strcat('learn2_',conds{3}),...
    strcat('test_',conds{1}),strcat('test_',conds{2}),strcat('test_',conds{3})}';



% test