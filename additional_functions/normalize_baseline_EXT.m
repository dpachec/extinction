function [EEG] = normalize_baseline_EXT(EEG, bline)

    if ndims(EEG.data) == 3
        for chani = 1:size(EEG.data,1)
            for triali = 1:size(EEG.data, 3)
                 data = EEG.data(chani, :, triali); 
                 mT = mean(data(bline),'omitnan');
                 stdT = std(data(bline),[], 'omitnan');
                 dataNorm = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);  
    
    
                EEG.data(chani, :, triali) = dataNorm; 
    
            end
        end
    elseif ndims(EEG.data) == 4
        for chani = 1:size(EEG.data,2)
            for triali = 1:size(EEG.data, 1)
                 data = EEG.data(triali, chani, :, :); 
                 mT = mean(data(:,:,:,bline),4,'omitnan');
                 stdT = std(data(:,:,:,bline),[], 4, 'omitnan');
                 dataNorm = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);  
    
    
                EEG.data(triali,chani, :, :) = dataNorm; 
    
            end
        end
    end


 
end
























                