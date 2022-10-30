function [EEG] = normalize_WM(EEG)

    mT = mean(EEG.power,1, 'omitnan');
    stdT = std(EEG.power,[], 1, 'omitnan');
    EEG.power = bsxfun(@rdivide, EEG.power- mT, stdT);  
    

% %     for chani = 1:size(EEG.data, 1)
% %         dataChan = squeeze(EEG.data(:, chani, :, :));
% %         mT = mean(dataChan,2, 'omitnan');
% %         stdT = std(dataChan,[], 2, 'omitnan');
% %         EEG.data(chani, :, :) = bsxfun(@rdivide, dataChan - mT, stdT);  
% %     end

 
end
























                