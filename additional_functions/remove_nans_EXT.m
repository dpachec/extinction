function[data] = remove_nans_EXT(data)

count = 1; 
for triali = 1:size(data, 1)
    dtr = squeeze(data(triali, :, :,:)); 
    if find(isnan(dtr))
        id2rem(count,:) = triali; 
        count = count+1;
    end
end

if exist('id2rem')
    data(id2rem, :, :, :) = []; 
end