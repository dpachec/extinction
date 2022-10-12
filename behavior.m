%% analyze behav data

clear 

load_paths
mkdir(path_out)
allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
           'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
           'c_sub17','c_sub18','c_sub19','c_sub20', 'c_sub21','c_sub22'};
       
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
 

 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
%% 