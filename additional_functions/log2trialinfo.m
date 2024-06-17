%% This scripts reads the raw log files and creates trialinfo files 
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

clear
paths = load_paths_EXT; 

% allsubs = {'c_sub01','c_sub02','c_sub03','c_sub04','c_sub05','c_sub06','c_sub07','c_sub08', ...
%            'c_sub09','c_sub10','c_sub11','c_sub12','c_sub13','c_sub14','c_sub15','c_sub16', ...
%            'c_sub17','c_sub18', 'c_sub19','c_sub20','c_sub21', 'c_sub22', 'c_sub23','c_sub24', ...
%             'c_sub25','c_sub26','c_sub27','c_sub28','c_sub29','c_sub30' }';

allsubs = {'p_sub21'}


for sub=1:length(allsubs)
    
    clearvars -except allsubs sub paths
    sel_sub=allsubs{sub};
    sel_sub_str=str2double(sel_sub(6:7));
    path_in=strcat(paths.raw_data,sel_sub,'/log/');
    all_logs= dir(path_in);
    all_logs={all_logs.name}; 
    allfiles={'A&B','C'};

for f=1:numel(allfiles) 
sel_file=allfiles{f}
    switch sel_file
        case 'A&B'
            allruns={'A','B'};
            allblocks=[1,2];
        case 'C'
             allruns={'C'};
             allblocks=3;
    end

for r=1:numel(allruns)

    sel_phase=allruns{r};
    sel_block=allblocks(r);
    file_ind=strncmp(strcat(num2str(sel_sub_str),'-phase',sel_file,'_'),all_logs,9);
    file_name=strcat(path_in,all_logs{file_ind});
    
    log = log2matSHORT(file_name); % info from logfile_x is now in this variable;  
    clear file_ind 

        start_trial = strncmp(strcat(sel_phase,'trial'),log{1,3},6);
        id_trials = find(start_trial);
        log_time_trials=log{1,4}(start_trial);
         phase_ind = ((r-1)*numel(id_trials))+1:r*numel(id_trials); % 1-72 or 73-144

        item_trials = find(strncmp('Item',log{1,3},4));
        item_trials = item_trials(phase_ind);            
        trlinfo(:,1) = log{1,1}(id_trials); % logfile trl number in first colum of trialinfo
        trlinfo(:,2)=ones(size(trlinfo)).*sel_block; % % number of task block
        log_time_item=log{1,4}(item_trials);
        
        % indices contexts, converts context into double, writes in trlinfo
        z = find(strncmp('context',log{1,3}, 7));
        z = z(phase_ind);
        log_time_video=log{1,4}(z);

        vid1 = cell(1, length(z));
        vid2 = cell(1, length(z));

            for i=1:numel(z)
                % context first index (context1_4)
                vid1{i}=  log{1,3}{z(i)}(8);
                % context second index
                vid2{i}=  log{1,3}{z(i)}(10);
            end
        vid=cellfun(@str2num,strcat(vid1,vid2),'UniformOutput',false);
        trlinfo(:,3) = [vid{:}];
        clear z vid vid1 vid2 

        % indices items, converts items into double, writes in trlinfo
        % indices types, converts types into double, writes in trlinfo
        % z = find(strncmp('Item',log{1,3}, 4));
        % z = z(phase_ind);
        % z not needed, we've got item_trials
        item = zeros(length(item_trials), 1);
        type = zeros(length(item_trials), 1);

        for i=1:numel(item_trials)
            item(i)=str2double(log{1,3}{item_trials(i)}(5));
        %    type(i)=str2double(log{1,3}{item_trials(i)}(end));
         type(i)=0; % type coding is off, recode at the end
        end
        trlinfo(:,5) = item;
        trlinfo(:,6) = type;
            clear item type i 

    % --------------------
    % until here, same A-D
    % --------------------

    if f == 1 % only in Phase A and B and all response trials now, therefore...

        % responses (1,2,3 or 4 ); participant probably didn't respond! (or responded with SPACE or responded too early (before response screen) or responded no matter what....
        % if no answer to response trial: NaN
        % ... probably switched responses...
        % always take first response
        y = zeros(length(trlinfo),1);
        RT = zeros(length(id_trials),1);
        for i = 1:length(id_trials)
            if (id_trials(i)+7)<numel(log{1,3})

           tmp_tr= log{1,3}(id_trials(i)+3:id_trials(i)+7);
           tmp_rt= log{1,4}(id_trials(i)+3:id_trials(i)+7);
            else (id_trials(i)+7)>numel(log{1,3})
           tmp_tr= log{1,3}(id_trials(i)+3:end);
           tmp_rt= log{1,4}(id_trials(i)+3:end);
            end

           tmp_allres=str2double(tmp_tr);
           res_ind= find(~isnan(tmp_allres));
           if isempty(res_ind)
                y(i)= NaN; % didn't respond when response_trial -> missing
                RT(i)=NaN;
           else
                y(i)=tmp_allres(res_ind(1)); 
                RT(i)=(tmp_rt(res_ind(1)));
           end

        end
        trlinfo(:,7) = y;

        % RT (response - item trial)
        % NaN: missing

        % cs1 vs. cs0
        cs = strncmp('cs',log{1,3}, 2);
        cs = find(cs);
         cs = cs(phase_ind);
        u = zeros(length(id_trials),1);
        for i=1:numel(cs)
            u(i)=str2double(log{1,3}{cs(i)}(3));
        end
        trlinfo(:,8) = u;         
        % keep y!
        % us_event (1) or no_us(0)
        w = zeros(length(trlinfo),1);
        a1 = strncmp(log{1,3},'us_event',7);
        a2 = strncmp(log{1,3},'no_us',7);
        b = find(a1+a2);
        b = b(phase_ind);
        log_time_us=log{1,4}(b);
        for i = 1:length(b)
            if strcmp('no_us',log{1,3}(b(i)))
                w(i) = 0;
            elseif strcmp('us_event',log{1,3}(b(i)))
                w(i) = 1;
            end
        end
        trlinfo(:,9) = w;
       clear w b a1 a2 u cs

    elseif f == 2
        % responses and RTs
        response = strncmp(strcat(sel_phase,'response'),log{1,3}, 5); 
        id_response = find(response);
        resp = zeros(length(id_response), 1);
        %first_button_press = zeros(length(log{1,3}));
        RT = zeros(length(id_response), 1);
        y = zeros(length(trlinfo),1);
        for i = 1:length(id_response)
           tmp_tr= log{1,3}(id_response(i)-1:id_response(i)+1);
           tmp_rt= log{1,4}(id_response(i)-1:id_response(i)+1);
           tmp_allres=str2double(tmp_tr);
           res_ind= find(~isnan(tmp_allres));
           if isempty(res_ind)
                y(i)= NaN; % didn't respond when response_trial -> missing
                RT(i)=NaN;
           else
                y(i)=tmp_allres(res_ind(1)); 
                RT(i)=(tmp_rt(res_ind(1)));
           end
        end
        trlinfo(:,7) = y;

        % all no_us -> not in trialinfo
        % pattern ABC -> A=1; B =2; C=3;
        h = find(strncmp('A',log{1,3}, 1));
        test = cell(length(id_trials),1);
        type = cell(length(id_trials),1);
        for i=1:numel(h)
          test{i}=log{1,3}{h(i)}(1:3);
          type{i}=log{1,3}{h(i)}(end);
        end          
        test = strrep(test, 'A', '1');
        test = strrep(test, 'B', '2');
        test = strrep(test, 'C', '3');
        test = (cellfun(@str2num,test))-120;
        type = cellfun(@str2double, type);
        trlinfo(:,4) = test;
               
        clear h test %type
        a2 = strncmp(log{1,3},'no_us',7);     
        log_time_us=log{1,4}(a2);
    end
    trlinfo(:,10)=log_time_trials;
   trlinfo(:,11)= log_time_video;
    trlinfo(:,12)=log_time_item;
    trlinfo(:,13)=log_time_us;
    trlinfo(:,14) = RT;

    clear log_time_trials log_time_us log_time_video log_time_item RTclear;
           
   log_allblock{sel_block}=trlinfo;
   clear trlinfo  
        
end
end
% concatenate all logs to define different types:
% cs+/cs+=1;cs+/cs-=2;cs-/cs-=3;
trlinfo=vertcat(log_allblock{:});

type_sum=[24,12,0];

for i=1:3
sum_us(i)=sum(trlinfo(trlinfo(:,5)==i,9));
sel_type=find(sum_us(i)==type_sum);
trlinfo(trlinfo(:,5)==i,6)=sel_type;
end


%there was an error in instructing patient c_sub01, responses in block a&b
%need to be switch
if strcmp(sel_sub, 'c_sub01')
check=0;
blockAB=trlinfo(:,2)<=2;
trlinfo(blockAB,7)=5-trlinfo(blockAB,7);
elseif strcmp(sel_sub, 'c_sub01')
trlinfo(:,7)=5-trlinfo(:,7);
end

save(strcat(paths.trlinfo,sel_sub,'_trlinfo'),'trlinfo')

end










%%