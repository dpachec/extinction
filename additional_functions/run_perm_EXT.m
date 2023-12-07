%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, cond1(:, 3:22, 3:22), cond2(:, 3:22, 3:22));
%junts = cat(1, cond1(:, 4:18, 4:18), cond2(:, 4:18, 4:18));
%junts = cat(1, cond1(:, 4:13, 4:13), cond2(:, 4:13, 4:13));
%junts = cat(1, cond1(:, 51:end, 51:end), cond2(:, 51:end, 51:end));
%junts = cat(1, cond1(:, 51:151, 51:151), cond2(:, 51:151, 51:151));
%junts = cat(1, cond1(:, 26:225, 26:225), cond2(:, 26:225, 26:225));
%junts = cat(1, cond1(:, 26:125, 26:125), cond2(:, 26:125, 26:125));

clear max_clust_sum_perm
for permi = 1:nPerm
    
    [M,N] = size(realCondMapping);
    rowIndex = repmat((1:M)',[1 N]);
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);

    cond1P = junts(fakeCondMapping == 0, :,:);
    cond2P = junts(fakeCondMapping == 1, :,:);

    diffC = cond1P - cond2P; 
    [hP pP ciP tsP] = ttest(diffC); 
    hP = squeeze(hP); hP(isnan(hP)) = 0; tP = squeeze(tsP.tstat); 
    clear allSTs  
    clustinfoP = bwconncomp(hP);
    for pxi = 1:length(clustinfoP.PixelIdxList)
       allSTs(pxi) = sum(tP(clustinfoP.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2uP idP] = max(abs(allSTs));
        max_clust_sum_perm(permi,:) = allSTs(idP); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
if exist('tObs')
    allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
    p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm
else 
    disp('no clusters')
    p = 1; 
end
