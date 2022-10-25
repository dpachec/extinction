function [oneListTraces oneListIds] = epoch_WM (EEG, eLim)

disp ('>>>>> epoching...');

EEG = add_EEGLAB_fields(EEG);
countE =1; 


for i = 1:length(EEG.event)
    eve = strsplit(EEG.event(i).type);
    EEG_b = pop_epoch( EEG, {EEG.event(i).type}, eLim, 'newname', 'verbose', 'no', 'Continuous EEG Data epochs', 'epochinfo', 'yes');
    traces{countE} = EEG_b.data;
    oneListIds{countE,:} = EEG.event(i).type;
    countE = countE+1;
    
end

%find empty epochs at the beginning and end of the exp
ids = find(cellfun(@isempty,traces));
oneListIds(ids) = [];
traces(ids) = []; 

oneListTraces = cat(3,traces{:});  


disp ('>>>>> epoching done');


