%% analyze behav data

clear 

load_paths
mkdir(path_out)
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
            'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' }';
       
for subji=1:length(allsubs)
    
    clearvars -except allsubs path_data path_info path_out subji x_all x_nous y_all y_nous
    sel_sub=allsubs{subji};
    info_file=strcat(path_info,sel_sub,'_trlinfo');
    load(info_file)  

    %   get response curves for each condition
    conditions={'cs+cs+','cs+cs-','cs-cs-'};
    figure
    for condi=1:3 
        cond{condi}=trlinfo(trlinfo(:,6)==condi,:);
        subplot(2,3,condi)
        hold on 
        y=cond{condi}(:,7);
        y=smoothdata(y,'movmean',3);
        plot(y,'--x')
        y_all{subji,condi}=cond{condi}(:,7);
        x_all{subji,condi}=1:numel(y_all{subji,condi});
        %plot dots for every us
        if condi<=2
            ind_us=find(cond{condi}(:,9)==1);
            scatter(ind_us,ones(size(ind_us)),'r','o')
        end 
        title([sel_sub, conditions{condi}])
   end
     
   for condi=1:3
       cond{condi}=trlinfo(trlinfo(:,6)==condi,:);
       subplot(2,3,4)
       hold on 
       x=find(cond{condi}(:,9)==0);
       x_nous{subji,condi}=x;
       y=cond{condi}(x,7);
       y_nous{subji,condi}=y;
       y=smoothdata(y,'movmean',3);
       plot(x,y,'--x')
       legend(conditions)
       title([sel_sub, conditions{condi}])
   end  
  
end
 




% %  % % % add nans for missing trials in last sub (this is for subj ***)
% %  x_all{end,1}=[x_all{end,1},NaN];
% %   y_all{end,1}=[y_all{end,1};NaN];
% %  x_all{end,2}=[x_all{end,2},NaN];
% %   y_all{end,2}=[y_all{end,2};NaN];
% %   
% % 
% %  x_nous{end,2}=[x_nous{end,2};NaN];
% %   y_nous{end,2}=[y_nous{end,2};NaN]; 
  
 % plot average response with and without
 
 %% 
 
 for condi = 1:3
    for subji=1:size(y_all,1)
         % for all use linspaced x
         y_all_mat(subji,:)=y_all{subji,condi};
         sm_y_all_mat(subji,:) = smoothdata(y_all{subji,condi},'movmean',3);
     end    
     y_all_mean(condi,:)=nanmean(y_all_mat);
     y_all_std(condi,:)=nanstd(y_all_mat);

     sm_y_all_mean(condi,:)=nanmean(sm_y_all_mat);
     sm_y_all_std(condi,:)=nanstd(sm_y_all_mat);

 end
 
 figure
 plot(y_all_mean')

 
 

 %%
 
 for c=1:size(y_all,2)
 for sub=1:size(y_all,1)
 % for all use linspaced x
 y_all_mat(sub,:)=y_all{sub,c};
 sm_y_all_mat(sub,:) = smoothdata(y_all{sub,c},'movmean',3);

 % for no use use average x
 y_nous_mat(sub,:)=y_nous{sub,c};
  sm_y_nous_mat(sub,:)=smoothdata(y_nous{sub,c},'movmean',3);

  x_nous_mat(sub,:)=x_nous{sub,c};

 end    
 y_all_mean(c,:)=nanmean(y_all_mat);
  y_all_std(c,:)=nanstd(y_all_mat);
x_nous_mean{c}=nanmean(x_nous_mat);
 y_nous_mean{c}=nanmean(y_nous_mat);
  y_nous_std{c}=nanstd(y_nous_mat);
  
   sm_y_all_mean(c,:)=nanmean(sm_y_all_mat);
  sm_y_all_std(c,:)=nanstd(sm_y_all_mat);
x_nous_mean{c}=nanmean(x_nous_mat);
 sm_y_nous_mean{c}=nanmean(sm_y_nous_mat);
  sm_y_nous_std{c}=nanstd(sm_y_nous_mat);
 clear y_all_mat y_nous_mat x_nous_mat  sm_y_all_mat sm_y_nous_mat
 end
 
 
 % plot boundedline plots
 figure
hold on
fig_stuff=subplot(2,2,1)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=1:size(y_all_mean,2);
    y1=y_all_mean(c,:);
    b1=y_all_std(c,:)./sqrt(size(y_all,1));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
 title('all responses')
 end
 
 fig_stuff=subplot(2,2,2)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=x_nous_mean{c};
    y1=y_nous_mean{c};
    b1=y_nous_std{c}./sqrt(size(y_all,1));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
  title('no post us responses')

 end
 
 fig_stuff=subplot(2,2,3)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=1:size(y_all_mean,2);
    y1=sm_y_all_mean(c,:);
    b1=sm_y_all_std(c,:)./sqrt(size(y_all,1));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
 title('all responses smoothed')
 end
 
 fig_stuff=subplot(2,2,4)
cmap_default=fig_stuff.ColorOrder;

 for c=1:3
     x1=x_nous_mean{c};
    y1=sm_y_nous_mean{c};
    b1=sm_y_nous_std{c}./sqrt(size(y_all,1));
   
boundedline(x1, y1, b1, 'cmap',cmap_default(c,:),'transparency',0.1,'alpha');
  title('no post us responses smoothed')

 end
 

 
 
 
 
 
 %% get some FIGURES for behavioral data per participant
path_info='D:\Extinction\iEEG\data\preproc\data_check\behav\'; % define where you want to save your file:
path_in='D:\Extinction\iEEG\data\preproc\trialinfo\';

mkdir(path_info)
% 
% allsubs = {'p_sub01','p_sub02','p_sub03','p_sub04','p_sub05','p_sub06','p_sub07',...
%           'c_sub01','c_sub02','c_sub03','c_sub05','c_sub06','c_sub07','c_sub08','c_sub09','c_sub10',...
%           'c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16','c_sub17','c_sub18'};
 %allsubs = {'c_sub04','c_sub20'};
 
 allsubs = {'p_sub10'}
 %allsubs = {'c_sub19','c_sub21','c_sub22','c_sub23','c_sub24','c_sub25'};

for sub = 1:numel(allsubs)
    sel_sub=allsubs{sub};
    load(strcat(path_in,sel_sub,'_trlinfo'))
    
    % organize trialinfo in matrix
    trialinfo=trlinfo;
    
    res = 4; % rating 1-4
    label = {'Dangerous','Safe'}; % 
    types = {'cs+/cs+','cs+/cs-','cs-/cs-'}; % only 3 types in EXT loss (cs+/cs+ = type1, cs+/cs- = type 2,cs-/cs- = type3) 
    items = 1:3; % only 3 items in EXT loss (2, 4, 6)
    block_def=[1 24;25 48;49 64]; % #of trials in blocks
    
    % which item is which type?
    type_item = zeros(1,length(items));
    responses_item = zeros(3, block_def(end,end)); 
    RT_item_alltrials = zeros(3, block_def(end,end));
    for i = 1:length(items)
        type_item(i) = mean(trialinfo(trialinfo(:,5) == items(i),6)); % mean gives exact value of itemtype
        responses_item(i,:) = trialinfo(trialinfo(:,5) == i,7); % give response from column 18 at (column 17 equals 1 or 2 or 3 (itemtype)
    
    end
    
    % compute average over types
    for ty=1:numel(types)
    RT_type_alltrials(ty,:)= trialinfo(trialinfo(:, 6) == ty, 14)-trialinfo(trialinfo(:, 6) == ty, 12)./10;
    responses(ty,:)=trialinfo(trialinfo(:,6) == ty,7);
    end
    
    
    % gives for every itemtype the responses throughout the experiment 
    
    % 1st line: type 1, 2nd line: type 2, 3rd line: type 3
    % gives reaction time for every itemtype (1st line: type 1, 2nd line: type 2, 3rd line: type 3)
    
    avg_response_type = nanmean(responses,2); % mean resonses per item type
    % 1st line: type 1, 2nd line: type 2, 3rd line: type 3
    % replace following line above if interested in median responses per item_type
    % avg_response_type = median(responses,2, 'omitnan');
    avg_RT_type_alltrials = nanmean(RT_type_alltrials,2); % mean RT per item type
    std_RT_type_alltrials = nanstd(RT_type_alltrials,0,2); % std RT per item type
    
    % median and mean responses and RT per BLOCK for every itemtype
    response_medianblock = zeros(length(items),size(block_def, 1));
    response_avgblock = zeros(length(items),size(block_def, 1));
    response_stdblock = zeros(length(items),size(block_def, 1));
    RT_avgblock = zeros(length(items),size(block_def, 1));
    RT_stdblock = zeros(length(items),size(block_def, 1));
    
    for b = 1:length(types)
        for i = 1:size(block_def, 1) % 4columns
            response_medianblock(b,i) = median(responses(b,(block_def(i,1):block_def(i,2))), 2, 'omitnan');
            response_avgblock(b,i) = nanmean(responses(b,(block_def(i,1):block_def(i,2))), 2);
            response_stdblock(b,i)= nanstd(responses(b,(block_def(i,1):block_def(i,2))));
            RT_avgblock(b,i) = nanmean(RT_type_alltrials(b,(block_def(i,1):block_def(i,2))), 2);
            RT_stdblock(b,i)= nanstd(RT_type_alltrials(b,(block_def(i,1):block_def(i,2))));
        end
    end
    
    clear tmp    
    figure % initializes figure
    
        for i = 1:length(types)
            subplot(3,4,i)
            hold on 
            rectangle('Position',[0 0 24 res],'FaceColor',[1 .9 .9])
            rectangle('Position',[24 0 24 res],'FaceColor',[0.9 .9 .9])
            rectangle('Position',[48 0 16 res],'FaceColor',[1 .8 1])
               
            stem(responses(i,:))
            title(types{i})
            h = gca;
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
        end % creates first 3 figures with colored stem plot
    
        % mean responses per itemtype
        subplot(3,3,4)
        plot(avg_response_type)
        title('mean responses per type')
        hold on
        h = gca;
        h.XTick = 1:3;
        h.XTickLabel = types;
        h.YLim = [0 4];
        h.YTick = 0:res;
        h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
        h.YTickLabelRotation = 90;
        
        % median responses per itemtype per block
        for i = 1:length(types)
            subplot(3,3,5)
            hold on 
            plot(response_medianblock(i,:))
            title('median rating/block')
            h = gca;
            h.XTick = 1:4;
            h.YLim = [0 4];
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
            %errorbar(response_avgblock(i,:),response_stdblock(i,:))
        end
        legend(types)
    
        % mean rating for each block
        for ty = 1:length(types)
            subplot(3,3,6)
            hold on
            plot(response_avgblock(ty,:));
            %errorbar(1:4,response_avgblock(ty,:),response_stdblock(ty,:))
            title('mean rating/block')
            h = gca;
            h.XTick = 1:4;
            h.YLim = [0 4];
            h.YTick = 0:res;
            h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
            h.YTickLabelRotation = 90;
        end
        legend(types)
        
        % mean RT per itemtype
        subplot(3,3,7)
        hold on
        bar(1:length(types),avg_RT_type_alltrials)
        errorbar(1:length(types),avg_RT_type_alltrials,std_RT_type_alltrials,'.')
        title('mean RT per type')
        h = gca;
        h.XTick = 1:3;
        h.XTickLabel = types;
        
        % mean RT per itemtype
        for ty = 1:length(types)
            subplot(3,3,8)
            hold on
            plot(RT_avgblock(ty,:));
            %errorbar(1:4,RT_avgblock(ty,:),RT_stdblock(ty,:))
            title('mean RT/block')
            h = gca;
            h.XTick = 1:4;
        end
        legend(types)
           
   savefig (strcat(path_info, sel_sub, '_rating_summary.fig'));   
    
%     figure
%     % plot rating for each item
%     for i=1:numel(items)
%     subplot(4,2,i)
%      hold on 
%             rectangle('Position',[0 0 16 res],'FaceColor',[1 .9 .9])
%             rectangle('Position',[16 0 16 res],'FaceColor',[0.9 .9 .9])
%             rectangle('Position',[32 0 16 res],'FaceColor',[1 .8 1])
%             rectangle('Position',[48 0 16 res],'FaceColor',[1 1 .8])
%               
%             stem(responses_item(i,:))
%             title(types{type_item(i)})
%             h = gca;
%             h.YTick = 0:res;
%             h.YTickLabel = {' ', label{1}, ' ', ' ',label{2}};
%             h.YTickLabelRotation = 90;
%     end
    
end
 
 
 
 
 
 
 
 
 
 
%% 