function [out_contrasts] = create_contrasts_EXT (cfg)
 
contr2save              =       cfg.contr2sav;
oneListIds              =       cfg.oneListIds;
oneListPow              =       cfg.oneListPow;
batch_bin               = 500; 

disp ('>>>>> creating contrasts');
 
allContrasts = []; allContrastIds = [];

for ci = 1:length(contr2save)
    eval(['count' contr2save{ci} '= 1;'])
    eval([contr2save{ci} '= [];'])
end

cr2c = cellfun(@(x) strsplit(x,'-'), contr2save, 'un', 0);
c2c = unique(cellfun(@(x) x{1}(end), cr2c, 'un', 0));

for i = 1:length(oneListIds)
    
    evei = double(string(oneListIds(i,:)));
        
% % % % % % ACQUISITION
if ~isempty(intersect(c2c, 'A'))
    if evei(2) == 1
        for j = 1:length(oneListIds) % repetitions are needed so should not start at i
           evej =  double(string(oneListIds(j,:)));
           if ~(evei(1) == evej(1))
                    if evej(2) == 1
                
                       if ~(evei(5)== evej(5)) & (evei(8)== evej(8)) % DISV
                            %disp(join(['DISVA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countDISVA')
                                new_disva{countDISVA} = [i, j];
                                countDISVA = countDISVA+1;
                             end
                       end
                       if ~(evei(5)== evej(5)) & ~(evei(8)== evej(8)) % DIDV
                            %disp(['DISVA> ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countDIDVA')
                                new_didva{countDIDVA} = [i, j];
                                countDIDVA = countDIDVA+1;
                             end
                       end
                       if (evei(5)== evej(5)) & evei(8)== 1  & evej(8)== 1% SICSPA
                            %disp(join(['SICSP > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countSICSPA')
                                new_sicspa{countSICSPA} = [i, j];
                                countSICSPA = countSICSPA+1;
                             end
                       end
                      if (evei(5)== evej(5)) & evei(8)== 0 & evej(8)== 0 % SICSMA
                             %disp(join(['SICSM > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countSICSMA')
                                new_sicsma{countSICSMA} = [i, j];
                                countSICSMA = countSICSMA+1;
                             end
                       end


                end
           end
        end
    end
end


% % % % % % ACQUISITION
if ~isempty(intersect(c2c, 'E'))
    if evei(2) == 2
        for j = 1:length(oneListIds) % repetitions are needed so should not start at i
           evej =  double(string(oneListIds(j,:)));
           if ~(evei(1) == evej(1))
                    if evej(2) == 2
                
                       if ~(evei(5)== evej(5)) & (evei(8)== evej(8)) % DISV
                            %disp(join(['DISVA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countDISVE')
                                new_disve{countDISVE} = [i, j];
                                countDISVE = countDISVE+1;
                             end
                       end
                       if ~(evei(5)== evej(5)) & ~(evei(8)== evej(8)) % DIDV
                            %disp(['DISVA> ' oneListIds{i} '//' oneListIds{j}]);   
                             if exist('countDIDVE')
                                new_didve{countDIDVE} = [i, j];
                                countDIDVE = countDIDVE+1;
                             end
                       end
                       if (evei(5)== evej(5)) & evei(8)== 1  & evej(8)== 1% SICSPA
                            %disp(join(['SICSP > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countSICSPE')
                                new_sicspe{countSICSPE} = [i, j];
                                countSICSPE = countSICSPE+1;
                             end
                       end
                      if (evei(5)== evej(5)) & evei(8)== 0 & evej(8)== 0 % SICSMA
                             %disp(join(['SICSM > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                             if exist('countSICSME')
                                new_sicsme{countSICSME} = [i, j];
                                countSICSME = countSICSME+1;
                             end
                       end


                end
           end
        end
    end
    end



end
    
 

for ci = 1:length(contr2save)
       
    eval(['if exist(''new_' lower(contr2save{ci}) ''')    new_' lower(contr2save{ci}) ' =  new_' lower(contr2save{ci}) ''' ; '  ...
          'new_' lower(contr2save{ci}) ' = vertcat(new_' lower(contr2save{ci}) '{:}); end' ])   

end 


for ci = 1:length(contr2save)
eval(  ['if exist(''new_' lower(contr2save{ci}) ''') & any(strcmp(contr2save, ''' contr2save{ci} ''')) ...' newline ...
           'disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
           ' nTrials = length(new_' lower(contr2save{ci}) '); ' newline ...
           ' if nTrials > batch_bin ' newline ...
               ' nBatches = ceil(nTrials/batch_bin);' newline ... 
               contr2save{ci} '{nBatches} = [];' newline ...
               ' for batchi = 1:nBatches ' newline ... 
                   'if strcmp(''' contr2save{ci}(end-1:end) ''', ''EE'')' newline ... 
                       ' if batchi == 1 ' newline ... 
                            't = batchi:batchi*batch_bin;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,1), :, :, 1:25); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,2), :, :, 1:25); ...' newline ... 
                        ' elseif (batchi)*batch_bin < nTrials ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,1), :, :, 1:25); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,2), :, :, 1:25); ...' newline ... 
                            ' else ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 1), :, :, 1:25); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 2), :, :, 1:25); ...' newline ... 
                            'end' newline ... 
                 'else' newline ... 
                      ' if batchi == 1 ' newline ... 
                        't = batchi:batchi*batch_bin;' newline ... 
                        contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,1), :, :, :); ...' newline ... 
                        contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t,2), :, :, :); ...' newline ... 
                    ' elseif (batchi)*batch_bin < nTrials ' newline ... 
                            't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,1), :, :, :); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow(new_' lower(contr2save{ci}) '(t,2), :, :, :); ...' newline ... 
                        ' else ' newline ... 
                            't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 1), :, :, :); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(t, 2), :, :, :); ...' newline ... 
                        'end' newline ... 
               'end' newline...               
            'end' newline ... 
            'else' newline ...
                'if strcmp(''' contr2save{ci}(end-1:end) ''', ''EE'')' newline ... 
                    contr2save{ci} '{1}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,1), :, :, 1:25); ...' newline ... 
                    contr2save{ci} '{1}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,2), :, :, 1:25); ...' newline ... 
                'else' newline ... 
                    contr2save{ci} '{1}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,1), :, :, :); ...' newline ... 
                    contr2save{ci} '{1}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,2), :, :, :); ...' newline ... 
                'end' newline ...
            'end' newline ... 
            contr2save{ci} '_ID = [oneListIds(new_' lower(contr2save{ci}) '(:,1)) oneListIds(new_' lower(contr2save{ci}) '(:,2))];' newline ...
            'end'    ]      );
       
end


        


for ci = 1:length(contr2save)
    
    
allContrastIds = contr2save'; 
for i = 1:length(allContrastIds)
    if exist( allContrastIds{i} ) 
        allContrasts{i} = eval([allContrastIds{i}]);
    end
end
allContrasts = allContrasts';
 
for i = 1:length(allContrastIds)
    %[allContrastIds{i} '_ID']
    % exist( [allContrastIds{i} '_ID']) 
    if exist( [allContrastIds{i} '_ID']) 
        allIDs{i} = eval([allContrastIds{i} '_ID']);
    else
        allIDs{i} = [];
    end
    
end
allIDs = allIDs';
 
 
out_contrasts.allContrasts              = allContrasts;
out_contrasts.allContrastIds            = allContrastIds ;
out_contrasts.allIDs                    = allIDs; 
 
 
end
 
%disp ('>> conditions extracted');
 
 
%%
 
 
 

