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
TEST = (allComb(:,3) == 3 & allComb(:,13) == 3); 
DI2  = (allComb(:,6) ~= allComb(:,16)); 
SI2  = (allComb(:,6) == allComb(:,16));
SV  = (allComb(:,9) == allComb(:,19)); 
DV  = (allComb(:,9) ~= allComb(:,19));
CSP = (allComb(:,9) == 1 & allComb(:,19) == 1);
CSM = (allComb(:,9) == 0 & allComb(:,19) == 0);
CSPP = (allComb(:,7) == 1 & allComb(:,17) == 1);
CSPM = (allComb(:,7) == 2 & allComb(:,17) == 2);
CSMM = (allComb(:,7) == 3 & allComb(:,17) == 3);
SC2 = (allComb(:,4) == allComb(:,14));
DC2 = (allComb(:,4) ~= allComb(:,14));
LTA = (allComb(:,2) > 36); 
ETE = (allComb(:,2) < 111); 
SPHA = allComb(:,3) == allComb(:,13);
EPP = (allComb(:,7) == 1 & allComb(:,17) == 3) | (allComb(:,7) == 3 & allComb(:,17) == 1);
EPM = (allComb(:,7) == 1 & allComb(:,17) == 2) | (allComb(:,7) == 2 & allComb(:,17) == 1);


if exist('SICSPT')
   idF = DT & SI2 & CSP & TEST;
   new_sicspt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMT')
   idF = DT & SI2 & CSM & TEST;
   new_sicspt = [allComb(idF, 1) allComb(idF, 11)];
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
if exist('DISVAESPHA')
   idF = DT & DC2 & SV & SPHA & (ACQ | EXT); 
   new_disvaespha = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DISVAE')
   idF = DT & DI2 & SV & (ACQ | EXT); 
   new_disvae = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DIDV')
   idF = DT & DI2 & DV; 
   new_didv = [allComb(idF, 1) allComb(idF, 11)];
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
if exist('DCCSPPE')
   idF = DT & DC2 & CSPP & EXT;
   new_dccsppe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPME')
   idF = DT & DC2 & CSPM & EXT;
   new_dccspme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSMME')
   idF = DT & DC2 & CSMM & EXT;
   new_dccsmme = [allComb(idF, 1) allComb(idF, 11)];
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
   idF = DT & SI2 & CSP & ACQ & LTA;
   new_sicspalt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMA')
   idF = DT & SI2 & CSM & ACQ;
   new_sicsma = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMALT')
   idF = DT & SI2 & CSM & ACQ & LTA;
   new_sicsmalt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPEET')
   idF = DT & SI2 & CSP & EXT & ETE;
   new_sicspeet = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMEET')
   idF = DT & SI2 & CSM & EXT & ETE;
   new_sicsmeet = [allComb(idF, 1) allComb(idF, 11)];
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
if exist('SCCSPPE')
   idF = DT & SC2 & CSPP & EXT;
   new_sccsppe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPPE')
   idF = DT & DC2 & CSPP & EXT;
   new_dccsppe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSPME')
   idF = DT & SC2 & CSPM & EXT;
   new_sccspme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSPME')
   idF = DT & DC2 & CSPM & EXT;
   new_dccspme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SCCSMME')
   idF = DT & SC2 & CSMM & EXT;
   new_sccsmme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCCSMME')
   idF = DT & DC2 & CSMM & EXT;
   new_dccsmme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DICSPA')
   idF = DT & DI2 & CSP & ACQ;
   new_dicspa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPE')
   idF = DT & SI2 & CSP & EXT;
   new_sicspe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSME')
   idF = DT & SI2 & CSM & EXT;
   new_sicsme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DICSME')
   idF = DT & DI2 & CSM & EXT;
   new_dicsme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMPP')
   idF = DT & SI2 & CSPM & EXT;
   new_sicsmpp = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMPM')
   idF = DT & SI2 & CSMM & EXT;
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
if exist('DCSPHA')
   idF = DT & DC2 & SPHA; 
   new_dcspha = [allComb(idF, 1) allComb(idF, 11)];
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
if exist('SCT')
   idF = DT & SC2 & TEST;
   new_sct = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('DCT')
   idF = DT & DC2 & TEST;
   new_dct = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPPT')
   idF = DT & SI2 & CSPP & TEST;
   new_sicsppt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPMT')
   idF = DT & SI2 & CSPM & TEST;
   new_sicspmt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMMT')
   idF = DT & SI2 & CSMM & TEST;
   new_sicsmmt = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPPE')
   idF = DT & SI2 & CSPP & EXT;
   new_sicsppe = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPME')
   idF = DT & SI2 & CSPM & EXT;
   new_sicspme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMME')
   idF = DT & SI2 & CSMM & EXT;
   new_sicsmme = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SISCOE')
   idF = DT & SI2 & EXT  & (CSMM | CSPP); 
   new_siscoe =  [allComb(idF, 1) allComb(idF, 11)];
   
end
if exist('SIDCOE')
   idF = DT & SI2 & EXT & CSPM; 
   new_sidcoe =  [allComb(idF, 1) allComb(idF, 11)];   
end
if exist('SICSPPA')
   idF = DT & SI2 & CSPP & ACQ;
   new_sicsppa = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSPMA')
   idF = DT & SI2 & CSPM & ACQ;
   new_sicspma = [allComb(idF, 1) allComb(idF, 11)];
end
if exist('SICSMMA')
   idF = DT & SI2 & CSMM & ACQ;
   new_sicsmma = [allComb(idF, 1) allComb(idF, 11)];
end



%if strcmp(cfg.tyRSA, 'POW')
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
% elseif strcmp(cfg.tyRSA, 'TR') | strcmp(cfg.tyRSA, 'PHA')  | strcmp(cfg.tyRSA, 'PLV') 
%     for ci = 1:length(contr2save)
%         eval(  ['if exist(''new_' lower(contr2save{ci}) ''') & any(strcmp(contr2save, ''' contr2save{ci} ''')) ...' newline ...
%                'disp ([''new_' lower(contr2save{ci})  ' ' ' '' num2str(length(new_' lower(contr2save{ci}) '))]);' newline...
%                ' nTrials = length(new_' lower(contr2save{ci}) '); ' newline ...
%                ' if nTrials > batch_bin ' newline ...
%                    ' nBatches = ceil(nTrials/batch_bin);' newline ... 
%                    contr2save{ci} '{nBatches} = [];' newline ...
%                    ' for batchi = 1:nBatches ' newline ... 
%                           ' if batchi == 1 ' newline ... 
%                             't = batchi:batchi*batch_bin;' newline ... 
%                             contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t,1), :, :); ...' newline ... 
%                             contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t,2), :, :); ...' newline ... 
%                         ' elseif (batchi)*batch_bin < nTrials ' newline ... 
%                                 't = ((batchi-1)*batch_bin)+1:batchi*batch_bin;' newline ... 
%                                 contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(t,1), :, :); ...' newline ... 
%                                 contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(t,2), :, :); ...' newline ... 
%                             ' else ' newline ... 
%                                 't = ((batchi-1)*batch_bin)+1:nTrials ;' newline ... 
%                                 contr2save{ci} '{batchi}(:,1,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t, 1), :, :); ...' newline ... 
%                                 contr2save{ci} '{batchi}(:,2,:,:) = oneListTraces (new_' lower(contr2save{ci}) '(t, 2), :, :); ...' newline ... 
%                             'end' newline ... 
%                 'end' newline ... 
%                 'else' newline ...
%                         contr2save{ci} '{1}(:,1,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(:,1), :, :); ...' newline ... 
%                         contr2save{ci} '{1}(:,2,:,:) = oneListTraces(new_' lower(contr2save{ci}) '(:,2), :, :); ...' newline ... 
%                 'end' newline ... 
%                 contr2save{ci} '_ID = [oneListIds(new_' lower(contr2save{ci}) '(:,1),:) oneListIds(new_' lower(contr2save{ci}) '(:,2),:)];' newline ...
%                 'end'    ]      );
%     end
% 
% 
% end


        


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
 
 
 

