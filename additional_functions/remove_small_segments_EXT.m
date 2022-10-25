function [data] = remove_small_segments_EXT(data, minSegLength)

    mks = ~(data); 
    clustinfo = bwconncomp(mks, 8); %see help bwconncomp for 8
    clust_info1 = cellfun(@numel,clustinfo.PixelIdxList);
    id2exc = find(clust_info1< minSegLength);
    
    for i=1:length(id2exc)
        data(clustinfo.PixelIdxList{id2exc(i)})=1;
    end

end



