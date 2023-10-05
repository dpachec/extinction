function [rsaZ] = normalize_rsaZ_EXT(rsaZ, bline)


    for chani = 1:size(rsaZ,1)
        for triali = 1:size(rsaZ, 2)
             data = rsaZ(chani, triali, :, :); 
             mT = mean(data(:,:,:,bline),4,'omitnan');
             stdT = std(data(:,:,:,bline),[],4, 'omitnan');
             dataNorm = bsxfun(@rdivide, bsxfun(@minus, data, mT), stdT);  


            rsaZ(chani, triali, :, :) = dataNorm; 

        end
    end


 
end
























                