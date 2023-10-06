function [out_contrasts] = create_contrasts_EXT (cfg)
 
contr2save              =       cfg.contr2sav;
oneListIds              =       cfg.oneListIds;
if isfield(cfg, 'oneListPow')
    oneListPow              =       cfg.oneListPow;
end
if isfield(cfg, 'oneListTraces')
    oneListTraces           =       cfg.oneListTraces;
end
batch_bin               =       500; 

disp ('>>>>> creating contrasts');
 
allContrasts = []; allContrastIds = [];

for ci = 1:length(contr2save)
    eval(['count' contr2save{ci} '= 1;'])
    eval([contr2save{ci} '= [];'])
end

cr2c = cellfun(@(x) strsplit(x,'-'), contr2save, 'un', 0);
c2c = unique(cellfun(@(x) x{1}(end), cr2c, 'un', 0));


oneListIdsB = double(string(oneListIds));

for i = 1:length(oneListIdsB)
    
    evei = oneListIdsB(i,:);
        
    
   for j = 1:length(oneListIdsB) % repetitions are needed so should not start at i
       evej =  oneListIdsB(j,:);
       if ~(evei(1) == evej(1))
       
           % % % % % % ACQUISITION
            if evei(2) == 1 & evej(2) == 1
            
                   if ~(evei(5)== evej(5)) & (evei(8)== evej(8)) % DISV
                        %disp(join(['DISVA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDISVA')
                            new_disva(countDISVA,:) = [i, j];
                            countDISVA = countDISVA+1;
                         end
                   end
                   if ~(evei(5)== evej(5)) & ~(evei(8)== evej(8)) % DIDV
                        %disp(['DISVA> ' oneListIds{i} '//' oneListIds{j}]);   
                         if exist('countDIDVA')
                            new_didva(countDIDVA,:)  = [i, j];
                            countDIDVA = countDIDVA+1;
                         end
                   end
                   if (evei(5)== evej(5)) & evei(8)== 1  & evej(8)== 1% SICSPA
                        %disp(join(['SICSP > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countSICSPA')
                            new_sicspa(countSICSPA,:)  = [i, j];
                            countSICSPA = countSICSPA+1;
                         end
                   end
                  if (evei(5)== evej(5)) & evei(8)== 0 & evej(8)== 0 % SICSMA
                         %disp(join(['SICSM > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countSICSMA')
                            new_sicsma(countSICSMA,:)  = [i, j];
                            countSICSMA = countSICSMA+1;
                         end
                  end
                  if ~(evei(5)== evej(5)) & (evei(3)== evej(3))  % DISC
                         %disp(join(['DISC > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDISCA')
                            new_disca(countDISCA,:)  = [i, j];
                            countDISCA = countDISCA+1;
                         end
                  end
                  if ~(evei(5)== evej(5)) & ~(evei(3)== evej(3))  % DISCA
                         %disp(join(['DIDCA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDIDCA')
                            new_didca(countDIDCA,:)  = [i, j];
                            countDIDCA = countDIDCA+1;
                         end
                  end
                    if (evei(3)== evej(3)) & evei(8)== 1  & evej(8)== 1 % SCCSPA
                         %disp(join(['SCCSPA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countSCCSPA')
                            new_sccspa(countSCCSPA,:)  = [i, j];
                            countSCCSPA = countSCCSPA+1;
                         end
                    end
                    if ~(evei(3)== evej(3)) & evei(8)== 1  & evej(8)== 1 % DCCSPA
                         %disp(join(['DCCSPA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDCCSPA')
                            new_dccspa(countDCCSPA,:)  = [i, j];
                            countDCCSPA = countDCCSPA+1;
                         end
                    end
                    if (evei(3)== evej(3)) & evei(8)== 0  & evej(8)== 0 % SCCSMA
                         %disp(join(['SCCSMA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countSCCSMA')
                            new_sccsma(countSCCSMA,:)  = [i, j];
                            countSCCSMA = countSCCSMA+1;
                         end
                    end
                    if ~(evei(3)== evej(3)) & evei(8)== 0  & evej(8)== 0 % DCCSMA
                         %disp(join(['DCCSMA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDCCSMA')
                            new_dccsma(countDCCSMA,:)  = [i, j];
                            countDCCSMA = countDCCSMA+1;
                         end
                    end
                    if (evei(3)== evej(3)) % SCA
                         %disp(join(['SCA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countSCA')
                            new_sca(countSCA,:)  = [i, j];
                            countSCA = countSCA+1;
                         end
                    end
                    if ~(evei(3)== evej(3)) % DCA
                         %disp(join(['DCA > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                         if exist('countDCA')
                            new_dca(countDCA,:)  = [i, j];
                            countDCA = countDCA+1;
                         end
                    end

                
           
    


    % % % % % % EXTINCTION
    elseif evei(2) == 2 & evej(2) == 2

                if ~(evei(5)== evej(5)) & (evei(8)== evej(8)) % DISVE
                    %disp(join(['DISVE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDISVE')
                        new_disve(countDISVE,:)  = [i, j];
                        countDISVE = countDISVE+1;
                     end
               end
               if ~(evei(5)== evej(5)) & ~(evei(8)== evej(8)) % DIDVE
                    %disp(['DISVE> ' oneListIds{i} '//' oneListIds{j}]);   
                     if exist('countDIDVE')
                        new_didve(countDIDVE,:)  = [i, j];
                        countDIDVE = countDIDVE+1;
                     end
               end
               if (evei(5)== evej(5)) & evei(8)== 1  & evej(8)== 1% SICSPE
                    %disp(join(['SICSP > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSPE')
                        new_sicspe(countSICSPE,:)  = [i, j];
                        countSICSPE = countSICSPE+1;
                     end
               end
              if (evei(5)== evej(5)) & evei(8)== 0 & evej(8)== 0 % SICSME
                     %disp(join(['SICSM > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSME')
                        new_sicsme(countSICSME,:)  = [i, j];
                        countSICSME = countSICSME+1;
                     end
              end
              if (evei(5)== evej(5)) & evei(6)== 2
                     %disp(join(['SICSMPPE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSMPPE')
                        new_sicsmppe(countSICSMPPE,:)  = [i, j];
                        countSICSMPPE = countSICSMPPE+1;
                     end
              end
              if (evei(5)== evej(5)) & evei(6)== 3
                     %disp(join(['SICSMPME > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSMPME')
                        new_sicsmpme(countSICSMPME,:)  = [i, j];
                        countSICSMPME = countSICSMPME+1;
                     end
              end
              if ~(evei(5)== evej(5)) & (evei(3)== evej(3))  % DISC
                     %disp(join(['DISCE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDISCE')
                        new_disce(countDISCE,:)  = [i, j];
                        countDISCE = countDISCE+1;
                     end
              end
              if ~(evei(5)== evej(5)) & ~(evei(3)== evej(3))  % DISC
                     %disp(join(['DIDCE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDIDCE')
                        new_didce(countDIDCE,:)  = [i, j];
                        countDIDCE = countDIDCE+1;
                     end
              end
              if (evei(3)== evej(3)) & evei(8)== 1  & evej(8)== 1 % SCCSPE
                     %disp(join(['SCCSPE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSCCSPE')
                        new_sccspe(countSCCSPE,:)  = [i, j];
                        countSCCSPE = countSCCSPE+1;
                     end
                end
                if ~(evei(3)== evej(3)) & evei(8)== 1  & evej(8)== 1 % DCCSPE
                     %disp(join(['DCCSPE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDCCSPE')
                        new_dccspe(countDCCSPE,:)  = [i, j];
                        countDCCSPE = countDCCSPE+1;
                     end
                end
                if (evei(3)== evej(3)) & evei(8)== 0  & evej(8)== 0 % SCCSME
                     %disp(join(['SCCSME > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSCCSME')
                        new_sccsme(countSCCSME,:)  = [i, j];
                        countSCCSME = countSCCSME+1;
                     end
                end
                if ~(evei(3)== evej(3)) & evei(8)== 0  & evej(8)== 0 % DCCSME
                     %disp(join(['DCCSME > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDCCSME')
                        new_dccsme(countDCCSME) = [i, j];
                        countDCCSME = countDCCSME+1;
                     end
                end
                if (evei(3)== evej(3)) % SCE
                     %disp(join(['SCE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSCE')
                        new_sce(countSCE,:)  = [i, j];
                        countSCE = countSCE+1;
                     end
                end
                if ~(evei(3)== evej(3)) % DCE
                     %disp(join(['DCE > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDCE')
                        new_dce(countDCE,:)  = [i, j];
                        countDCE = countDCE+1;
                     end
                end
           end
       
       



   % % % % % % RENEWAL
    elseif evei(2) == 3 & evej(2) == 3

                if ~(evei(5)== evej(5)) & ( (evei(6)== 1 & evej(6)== 1) | (evei(6)== 3 & evej(6)== 3) ) % DISVR
                    %disp(join(['DISVR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDISVR')
                        new_disvr(countDISVR,:)  = [i, j];
                        countDISVR = countDISVR+1;
                     end
               end
               if ~(evei(5)== evej(5)) & ( (evei(6)== 1 & evej(6)== 3) | (evei(6)== 3 & evej(6)== 1) ) % DIDVR
                    %disp(['DISVR> ' oneListIds{i} '//' oneListIds{j}]);   
                     if exist('countDIDVR')
                        new_didvr(countDIDVR,:)  = [i, j];
                        countDIDVR = countDIDVR+1;
                     end
               end
               if (evei(5)== evej(5)) & evei(6)== 1  & evej(6)== 1 %SICSPR
                    %disp(join(['SICSPR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSPR')
                        new_sicspr(countSICSPR,:)  = [i, j];
                        countSICSPR = countSICSPR+1;
                     end
               end
              if (evei(5)== evej(5)) & (evei(6)~=1 & evej(6)~= 1) 
                     %disp(join(['SICSMR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSMR')
                        new_sicsmr(countSICSMR,:)  = [i, j];
                        countSICSMR = countSICSMR+1;
                     end
              end
              if (evei(5)== evej(5)) & evei(6)== 2  & evej(6)== 2
                     %disp(join(['SICSMPPR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSMPPR')
                        new_sicsmppr(countSICSMPPR,:)  = [i, j];
                        countSICSMPPR = countSICSMPPR+1;
                     end
              end
              if (evei(5)== evej(5)) & evei(6)== 3  & evej(6)== 3
                     %disp(join(['SICSMPMR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSICSMPMR')
                        new_sicsmpmr(countSICSMPMR,:)  = [i, j];
                        countSICSMPMR = countSICSMPMR+1;
                     end
              end
              if ~(evei(5)== evej(5)) & (evei(3)== evej(3))  % DISCR
                     %disp(join(['DISCR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDISCR')
                        new_discr(countDISCR,:)  = [i, j];
                        countDISCR = countDISCR+1;
                     end
              end
              if ~(evei(5)== evej(5)) & ~(evei(3)== evej(3))  % DIDCR
                     %disp(join(['DIDCR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDIDCR')
                        new_didcr(countDIDCR,:)  = [i, j];
                        countDIDCR = countDIDCR+1;
                     end
              end
                if (evei(3)== evej(3)) % SCR
                     %disp(join(['SCR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countSCR')
                        new_scr(countSCR,:)  = [i, j];
                        countSCR = countSCR+1;
                     end
                end
                if ~(evei(3)== evej(3)) % DCR
                     %disp(join(['DCR > ' oneListIds(i,:) '//' oneListIds(j, :)],'_'));   
                     if exist('countDCR')
                        new_dcr(countDCR,:)  = [i, j];
                        countDCR = countDCR+1;
                     end
                end
       end



       
       
       end
   end
end


if strcmp(cfg.tyRSA, 'pRSA')
    for ci = 1:length(contr2save)
        eval(  ['if exist(''new_' lower(contr2save{ci}) ''') & any(strcmp(contr2save, ''' contr2save{ci} ''')) ...' newline ...
               'disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
               ' nTrials = length(new_' lower(contr2save{ci}) '); ' newline ...
               ' if nTrials > batch_bin ' newline ...
                   ' nBatches = ceil(nTrials/batch_bin);' newline ... 
                   contr2save{ci} '{nBatches} = [];' newline ...
                   ' for batchi = 1:nBatches ' newline ... 
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
                'end' newline ... 
                'else' newline ...
                        contr2save{ci} '{1}(:,1,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,1), :, :, :); ...' newline ... 
                        contr2save{ci} '{1}(:,2,:,:,:) = oneListPow (new_' lower(contr2save{ci}) '(:,2), :, :, :); ...' newline ... 
                'end' newline ... 
                contr2save{ci} '_ID = [oneListIds(new_' lower(contr2save{ci}) '(:,1),:) oneListIds(new_' lower(contr2save{ci}) '(:,2),:)];' newline ...
                'end'    ]      );
           
    end
elseif strcmp(cfg.tyRSA, 'tRSA')
    for ci = 1:length(contr2save)
        eval(  ['if exist(''new_' lower(contr2save{ci}) ''') & any(strcmp(contr2save, ''' contr2save{ci} ''')) ...' newline ...
               'disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
               ' nTrials = length(new_' lower(contr2save{ci}) '); ' newline ...
               ' if nTrials > batch_bin ' newline ...
                   ' nBatches = ceil(nTrials/batch_bin);' newline ... 
                   contr2save{ci} '{nBatches} = [];' newline ...
                   ' for batchi = 1:nBatches ' newline ... 
                          ' if batchi == 1 ' newline ... 
                            't = batchi:batchi*batch_bin;' newline ... 
                            contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t,1), :, :); ...' newline ... 
                            contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t,2), :, :); ...' newline ... 
                        ' elseif (batchi)*batch_bin < nTrials ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(t,1), :, :); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(t,2), :, :); ...' newline ... 
                            ' else ' newline ... 
                                't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
                                contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t, 1), :, :); ...' newline ... 
                                contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t, 2), :, :); ...' newline ... 
                            'end' newline ... 
                'end' newline ... 
                'else' newline ...
                        contr2save{ci} '{1}(:,1,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(:,1), :, :); ...' newline ... 
                        contr2save{ci} '{1}(:,2,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(:,2), :, :); ...' newline ... 
                'end' newline ... 
                contr2save{ci} '_ID = [oneListIds(new_' lower(contr2save{ci}) '(:,1),:) oneListIds(new_' lower(contr2save{ci}) '(:,2),:)];' newline ...
                'end'    ]      );
    end



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
 
 
 

