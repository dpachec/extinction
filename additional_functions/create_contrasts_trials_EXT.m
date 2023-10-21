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

oneListIdsA = [[1:size(oneListIds,1)]' double(string(oneListIds))];


CSP = oneListIdsA(:,9) == 1; 
CSM = oneListIdsA(:,9) == 0;
ACQ = oneListIdsA(:,3) == 1; 
EXT = oneListIdsA(:,3) == 2; 
REN = oneListIdsA(:,3) == 3;
CSMPP2 = oneListIdsA(:,7) == 2;
CSMPM2 = oneListIdsA(:,7) == 3;


IT1 = oneListIdsA(:,7) == 1;
IT2 = oneListIdsA(:,7) == 2;
IT3 = oneListIdsA(:,7) == 3;



if exist('ALLE')
   idF = EXT; 
   new_alle = [oneListIdsA(idF, 1)];
   new_alle(:, 2) = new_alle(:, 1) + 1; 
end
if exist('CSPA')
   idF = CSP & ACQ; 
   new_cspa = [oneListIdsA(idF, 1)];
   new_cspa(:, 2) = new_cspa(:, 1) + 1; 
end
if exist('CSMA')
   idF = CSM & ACQ; 
   new_csma = [oneListIdsA(idF, 1)];
   new_csma(:, 2) = new_csma(:, 1) + 1; 
end
if exist('CSPE')
   idF = CSP & EXT; 
   new_cspe = [oneListIdsA(idF, 1)];
   new_cspe(:, 2) = new_cspe(:, 1) + 1; 
end

if exist('CSME')
   idF = CSM & EXT; 
   new_csme = [oneListIdsA(idF, 1)];
   new_csme(:, 2) = new_csme(:, 1) + 1; 
end

if exist('CSMPP')
   idF = CSMPP2 & EXT; 
   new_csmpp = [oneListIdsA(idF, 1)];
   new_csmpp(:, 2) = new_csmpp(:, 1) + 1; 
end

if exist('CSMPM')
   idF = CSMPM2 & EXT; 
   new_csmpm = [oneListIdsA(idF, 1)];
   new_csmpm(:, 2) = new_csmpm(:, 1) + 1; 
end
if exist('IT1CSPA')
   idF = IT1 & ACQ; 
   new_it1cspa = [oneListIdsA(idF, 1)];
   new_it1cspa(:, 2) = new_it1cspa(:, 1) + 1; 
end
if exist('IT2CSPA')
   idF = IT2 & ACQ; 
   new_it2cspa = [oneListIdsA(idF, 1)];
   new_it2cspa(:, 2) = new_it2cspa(:, 1) + 1; 
end
if exist('IT2CSME')
   idF = IT2 & ACQ; 
   new_it2csme = [oneListIdsA(idF, 1)];
   new_it2csme(:, 2) = new_it1csme(:, 1) + 1; 
end
if exist('IT3CSME')
   idF = IT3 & ACQ; 
   new_it3csme = [oneListIdsA(idF, 1)];
   new_it3csme(:, 2) = new_it3csme(:, 1) + 1; 
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
 
 
 

