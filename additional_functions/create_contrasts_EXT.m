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

%disp ('>>>>> creating contrasts');
 
allContrasts = []; allContrastIds = [];

for ci = 1:length(contr2save)
    eval([contr2save{ci} '= [];'])
end

oneListIdsA = double(string(oneListIds));
idsIt= [1:size(oneListIdsA,1)]';
idsIt(:,2:10) = oneListIdsA(:, 1:9); 
ids = permn(1:size(oneListIdsA, 1),2); % both directions 
allComb = [idsIt(ids(:, 1),:) idsIt(ids(:,2),:)];

DT  = (allComb(:,2) ~= allComb(:,12)); 
ACQ = (allComb(:,3) == 1 & allComb(:,13) == 1); 
EXT = (allComb(:,3) == 2 & allComb(:,13) == 2); 
REN = (allComb(:,3) == 3 & allComb(:,13) == 3); 
DI2  = (allComb(:,6) ~= allComb(:,16)); 
SI2  = (allComb(:,6) == allComb(:,16));
SV  = (allComb(:,9) == allComb(:,19)); 
DV  = (allComb(:,9) ~= allComb(:,19));
CSP = (allComb(:,9) == 1 & allComb(:,19) == 1);
CSM = (allComb(:,9) == 0 & allComb(:,19) == 0);
CSMPP = (allComb(:,7) == 2 & allComb(:,17) == 2);
CSMPM = (allComb(:,7) == 3 & allComb(:,17) == 3);
SC2 = (allComb(:,4) == allComb(:,14));
DC2 = (allComb(:,4) ~= allComb(:,14));
LT = (allComb(:,2) > 39 & allComb(:,2) < 79) | (allComb(:,2) > 108 & allComb(:,2) < 154); 

EPP = (allComb(:,7) == 1 & allComb(:,17) == 3) | (allComb(:,7) == 3 & allComb(:,17) == 1);
EPM = (allComb(:,7) == 1 & allComb(:,17) == 2) | (allComb(:,7) == 2 & allComb(:,17) == 1);




if exist('ALLA')
   idF = ACQ; 
   new_alla = [oneListIdsA(idF, 1)];
   new_alla(:, 2) = new_alla(:, 1) + 1; 
end
if exist('ALLE')
   idF = DT & EXT;
   new_alle = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SI')
   idF = DT & SI2; 
   new_si = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DI')
   idF = DT & DI2; 
   new_di = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SIA')
   idF = DT & SI2 & ACQ; 
   new_sia = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIA')
   idF = DT & DI2 & ACQ; 
   new_dia = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SIE')
   idF = DT & SI2 & EXT;
   new_sie = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIE')
   idF = DT & DI2 & EXT; 
   new_die = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SISV')
   idF = DT & SI2 & SV; 
   new_sisv = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SISVAE')
   idF = DT & SI2 & SV & (ACQ | EXT); 
   new_sisvae = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SISVA')
   idF = DT & SI2 & SV & ACQ; 
   new_sisva = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SISVE')
   idF = DT & SI2 & SV & EXT;
   new_sisve = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DISV')
   idF = DT & DI2 & SV; 
   new_disv = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DISVA')
   idF = DT & DI2 & SV & ACQ; 
   new_disva = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DISVAE')
   idF = DT & DI2 & SV & (ACQ | EXT); 
   new_disvae = [allComb(idF, 1) allComb(idF, 11)];
end

if exist('DIDVA')
   idF = DT & DI2 & DV & ACQ; 
   new_didva = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DISVE')
  idF = DT & DI2 & SV & EXT; 
   new_disve = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIDVE')
    idF = DT & DI2 & DV & EXT; 
   new_didve = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIDVEPP') % % TAG
   idF = DT & DI2 & EPP & EXT; 
   new_didvepp = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIDVEPM') % % TAG
   idF = DT & DI2 & EPM & EXT; 
   new_didvepm = [allComb(idF, 1) allComb(idF, 11)];
end

if exist('SICSPAE')
   idF = DT & SI2 & CSP & (ACQ | EXT); 
   new_sicspae = [allComb(idF, 1) allComb(idF, 11)];
end

if exist('SICSMAE')
   idF = DT & SI2 & CSM & (ACQ | EXT); 
   new_sicsmae = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSP')
   idF = DT & SC2 & CSP; 
   new_sccsp = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSM')
   idF = DT & SC2 & CSM; 
   new_sccsm = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSP')
   idF = DT & DC2 & CSP; 
   new_dccsp = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSM')
   idF = DT & DC2 & CSM; 
   new_dccsm = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPA')
   idF = DT & DC2 & CSP & ACQ; 
   new_dccspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSMA')
   idF = DT & DC2 & CSM & ACQ; 
   new_dccsma = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPE')
   idF = DT & DC2 & CSP & EXT; 
   new_dccspe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSME')
   idF = DT & DC2 & CSM & EXT; 
   new_dccsme = [allComb(idF, 1) allComb(idF, 11)];
end

if exist('DCCSPA')
   idF = DT & DC2 & CSP & ACQ; 
   new_dccspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPE')
   idF = DT & DC2 & CSP & EXT; 
   new_dccspe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSPAE')
   idF = DT & SC2 & CSP & (ACQ | EXT); 
   new_sccspae = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSMAE')
   idF = DT & SC2 & CSM  & (ACQ | EXT); 
   new_sccsmae = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPA')
   idF = DT & SI2 & CSP & ACQ;
   new_sicspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPALT')
   idF = DT & SI2 & CSP & ACQ & LT;
   new_sicspalt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMA')
   idF = DT & SI2 & CSM & ACQ;
   new_sicsma = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMALT')
   idF = DT & SI2 & CSM & ACQ & LT;
   new_sicsmalt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSPA')
   idF = DT & SC2 & CSP & ACQ;
   new_sccspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSMA')
   idF = DT & SC2 & CSM & ACQ;
   new_sccsma = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSPE')
   idF = DT & SC2 & CSP & EXT;
   new_sccspe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSME')
   idF = DT & SC2 & CSM & EXT;
   new_sccsme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DICSPA')
   idF = DT & DI2 & CSP & ACQ;
   new_dicspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPE')
   idF = DT & SI2 & CSP & EXT;
   new_sicspe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPELT')
   idF = DT & SI2 & CSP & EXT & LT;
   new_sicspelt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSME')
   idF = DT & SI2 & CSM & EXT;
   new_sicsme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMELT')
   idF = DT & SI2 & CSM & EXT & LT;
   new_sicsmelt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DICSME')
   idF = DT & DI2 & CSM & EXT;
   new_dicsme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMPP')
   idF = DT & SI2 & CSMPP & EXT;
   new_sicsmpp = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMPM')
   idF = DT & SI2 & CSMPM & EXT;
   new_sicsmpm = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SC')
   idF = DT & SC2; 
   new_sc = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DC')
   idF = DT & DC2; 
   new_dc = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCA')
   idF = DT & SC2 & ACQ;
   new_sca = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCA')
   idF = DT & DC2 & ACQ;
   new_dca = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCE')
   idF = DT & SC2 & EXT;
   new_sce = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCE')
   idF = DT & DC2 & EXT;
   new_dce = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCR')
   idF = DT & SC2 & REN;
   new_scr = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCR')
   idF = DT & DC2 & REN;
   new_dcr = [allComb(idF, 1) allComb(idF, 11)];
end




if strcmp(cfg.tyRSA, 'POW')
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
elseif strcmp(cfg.tyRSA, 'TR') | strcmp(cfg.tyRSA, 'PHA')  | strcmp(cfg.tyRSA, 'PLV') 
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

ids = cellfun(@(x) double(string(x)),allIDs, 'un', 0);
out_contrasts.allIDs                    = ids; 
 
 
end
 
%disp ('>> conditions extracted');
 
 
%%
 
 
 

