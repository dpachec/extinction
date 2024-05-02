%%
%% Temporal RSA IN LOOP 
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast
clear , clc

listF2sav = {

'POW_PFC_C_3-54_1_0_50-10_1_SCA-DCA';

};   


t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        EEG = ALLEEG{subji};
        
        
        if ~isempty(EEG)

            
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); 

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = Ev2; 

            if strcmp(cfg.tyRSA, 'TR')
                EEG = normalize_baseline_EXT(EEG, [2501:3000]); 
                %EEG = normalize_EXT(EEG);  %across trials
                EEG = downsample_EEG_EXT(EEG); 
                cfg.oneListTraces = permute(EEG.data(:, 251:550,:), [3 1 2]); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'POW')
                EEG = extract_power_EXT(EEG, 0.01); 
                %EEG = normalize_baseline_EXT(EEG, [251:300]); 
                EEG = normalize_EXT(EEG);  %across trials
                cfg.oneListPow = EEG.power(:, :, : ,251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'PHA')
                EEG = normalize_EXT(EEG);  %across trials
                phaTS = extract_pha_EXT(EEG, cfg);
                cfg.oneListTraces = phaTS(:, :, 251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT3(out_contrasts, cfg);
                toc
            elseif strcmp(cfg.tyRSA, 'PLV')
                EEG = normalize_EXT(EEG);  %across trials
                phaTS = extract_pha_EXT(EEG, cfg);
                cfg.oneListTraces = phaTS(:, :, 251:550); 
                out_contrasts = create_contrasts_EXT(cfg);
                tic
                out_rsa(subji, :, :, :) = rsa_EXT5(out_contrasts, cfg);
                toc
            end
        
            ids{subji,:} = out_contrasts.allIDs; nnans{subji,:} = EEG.nan; 
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids', 'nnans');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end





%% plot TG
clear, clc
paths = load_paths_EXT; 
 
f2sav =  'POW_OCC_C_3-54_1_0_50-10_1_SCA-DCA';



sub2exc = [];


load ([ paths.results.rsa f2sav '.mat']);

out_rsa(sub2exc, :, :,:) = []; 
ids = rem_nan_subj_EXT(out_rsa); 

cond1 = squeeze(out_rsa(:, 1, 1:23, 1:23)); 
cond2 = squeeze(out_rsa(:, 2, 1:23, 1:23)); 
 %cond1 = squeeze(out_rsa(:, 1, 1:20, 1:20)); 
 %cond2 = squeeze(out_rsa(:, 2, 1:20, 1:20)); 


cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 
diff = cond1-cond2; 

[cond1 cond2] = rem_half_matrix(cond1, cond2);

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
%h(1:251,1:26) = 0;  % % % no clusters before baseline

clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    %[max2u id] = max((allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;


clim = [-.03 .03];
plot_TG_map(m1, m2, h, t, f2sav, clim)
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)




%% PERMUTATIONS
nPerm = 1000;

nSubj =  size(cond1, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, cond1(:, 3:20, 3:20), cond2(:, 3:20, 3:20));
%junts = cat(1, cond1(:, 3:22, 3:22), cond2(:, 3:22, 3:22));

[M,N] = size(realCondMapping);
rowIndex = repmat((1:M)',[1 N]);
    
clear max_clust_sum_perm
for permi = 1:nPerm
    
    [~,randomizedColIndex] = sort(rand(M,N),2);
    newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
    fakeCondMapping = realCondMapping(newLinearIndex);
    fakeCondMapping = fakeCondMapping(:);

    cond1P = junts(fakeCondMapping == 0, :,:);
    cond2P = junts(fakeCondMapping == 1, :,:);

    diffC = cond1P - cond2P; 
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm

%% PERMUTATIONS IN LOOP 
clear , clc


nPerm = 10000;



listF2sav = {

'POW_PFC_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_PFC_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_PFC_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_PFC_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_PFC_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_PFC_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_PFC_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_PFC_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_PFC_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_PFC_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_PFC_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_PFC_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_PFCO_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_PFCO_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_PFCO_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_PFCO_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_PFCO_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_PFCO_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_PFCO_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_PFCO_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_PFCO_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_PFCO_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_PFCO_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_PFCO_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_HPC_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_HPC_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_HPC_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_HPC_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_HPC_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_HPC_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_HPC_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_HPC_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_HPC_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_HPC_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_HPC_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_HPC_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_OFC_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_OFC_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_OFC_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_OFC_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_OFC_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_OFC_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_OFC_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_OFC_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_OFC_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_OFC_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_OFC_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_OFC_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_AMY_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_AMY_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_AMY_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_AMY_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_AMY_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_AMY_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_AMY_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_AMY_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_AMY_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_AMY_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_AMY_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_AMY_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_OCC_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_OCC_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_OCC_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_OCC_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_OCC_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_OCC_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_OCC_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_OCC_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_OCC_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_OCC_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_OCC_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_OCC_C_3-54_1_0_50-10_1_SCT-DCT';

'POW_TMP_C_3-54_1_0_50-10_1_SISVA-DISVA';
'POW_TMP_C_3-54_1_0_50-10_1_SISVE-DISVE';
'POW_TMP_C_3-54_1_0_50-10_1_DISVA-DIDVA';
'POW_TMP_C_3-54_1_0_50-10_1_DISVE-DIDVE';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPA-SICSMA';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPE-SICSME';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPE-SICSMPP';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPE-SICSMPM';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPPT-SICSPMT';
'POW_TMP_C_3-54_1_0_50-10_1_SICSPPT-SICSMMT';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSPA-DCCSPA';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSPE-DCCSPE';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSMA-DCCSMA';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSME-DCCSME';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSPA-SCCSMA';
'POW_TMP_C_3-54_1_0_50-10_1_SCCSPE-SCCSME';
'POW_TMP_V_3-54_1_0_50-10_1_SCA-DCA';
'POW_TMP_V_3-54_1_0_50-10_1_SCE-DCE';
'POW_TMP_V_3-54_1_0_50-10_1_SCT-DCT';
'POW_TMP_C_3-54_1_0_50-10_1_SCA-DCA';
'POW_TMP_C_3-54_1_0_50-10_1_SCE-DCE';
'POW_TMP_C_3-54_1_0_50-10_1_SCT-DCT';

};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths nPerm allFilesP
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);

    load ([ paths.results.rsa f2sav '.mat']);

    % % % % % % compute real t
    ids = rem_nan_subj_EXT(out_rsa); 
    f2s = strsplit(f2sav, '_'); 
    
    if strcmp(f2s{3}, 'C')
        cond1 = squeeze(out_rsa(:, 1, 1:19, 1:19)); 
        cond2 = squeeze(out_rsa(:, 2, 1:19, 1:19)); 
    else
        cond1 = squeeze(out_rsa(:, 1, 1:22, 1:22)); 
        cond2 = squeeze(out_rsa(:, 2, 1:22, 1:22)); 
    end
    cond1(ids, :, :) = []; 
    cond2(ids, :, :) = []; 
    diff = cond1-cond2; 
    
    [cond1 cond2] = rem_half_matrix(cond1, cond2);
    
    m1 = squeeze(mean(cond1, 'omitnan')); 
    m2 = squeeze(mean(cond2, 'omitnan')); 
    
    [h p ci ts] = ttest(cond1, cond2); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    %h(1:26,1:2) = 0;  % % % no clusters before baseline
        
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max((allSTs));
        tObs = allSTs(id); 
    end
    
    
    nSubj =  size(cond1, 1);
    realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';
    
    if strcmp(f2s{3}, 'C')
        junts = cat(1, cond1(:, 3:19, 3:19), cond2(:, 3:19, 3:19));
    else
        junts = cat(1, cond1(:, 3:19, 3:19), cond2(:, 3:19, 3:19));
    end
    [M,N] = size(realCondMapping);
    rowIndex = repmat((1:M)',[1 N]);

    clear max_clust_sum_perm
    for permi = 1:nPerm
        
       
        [~,randomizedColIndex] = sort(rand(M,N),2);
        newLinearIndex = sub2ind([M,N],rowIndex,randomizedColIndex);
        fakeCondMapping = realCondMapping(newLinearIndex);
        fakeCondMapping = fakeCondMapping(:);
    
        cond1P = junts(fakeCondMapping == 0, :,:);
        cond2P = junts(fakeCondMapping == 1, :,:);
    
        diffC = cond1P - cond2P; 
        [h p ci ts] = ttest(diffC); 
        h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
        clear allSTs  
        clustinfo = bwconncomp(h);
        for pxi = 1:length(clustinfo.PixelIdxList)
           allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
        end
        
        if exist('allSTs')
            [max2u id] = max(abs(allSTs));
            %[max2u id] = max(allSTs);
            max_clust_sum_perm(permi,:) = allSTs(id); 
        else
            max_clust_sum_perm(permi,:) = 0; 
        end
    
    end
    
    
    disp('done')
    
    if exist('tObs') 
        allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
        allFilesP{listi, 1} = f2sav;
        allFilesP{listi, 2} = 1 - ((nPerm-1) - (length (allAb)))  / nPerm;
    else
        allFilesP{listi, 1} = f2sav;
        allFilesP{listi, 2} = 1;
    end


end

save ('allFilesP', 'allFilesP')


%% plot all subjects


for subji = 1:33

figure(); 
tiledlayout(1,2)
nexttile
imagesc(squeeze(cond1(subji, :, :))); axis square; 
set(gca, 'clim', [-.1 .1])
nexttile
imagesc(squeeze(cond2(subji, :, :))); axis square; 
set(gca, 'clim', [-.1 .1])




end


%% plot TG in LOOP
clear , clc


listF2sav = {


'POW_PFC_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_PFC_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

'POW_HPC_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_HPC_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

'POW_OFC_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_OFC_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

'POW_AMY_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_AMY_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

'POW_OCC_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_OCC_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

'POW_TMP_C_3-54_1_0_50-10_1_DCCSPA-DCCSMA';
'POW_TMP_C_3-54_1_0_50-10_1_DCCSPE-DCCSME';

};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);
                    

    load ([ paths.results.rsa f2sav '.mat']);
    
    ids = rem_nan_subj_EXT(out_rsa); 
    cond1 = squeeze(out_rsa(:, 1, 1:26, 1:26)); 
    cond2 = squeeze(out_rsa(:, 2, 1:26, 1:26)); 
    %cond1 = squeeze(out_rsa(:, 1, 1:250, 1:250)); 
    %cond2 = squeeze(out_rsa(:, 2, 1:250, 1:250)); 
    %cond1 = squeeze(out_rsa(:, 1, 1:125, 1:125)); 
    %cond2 = squeeze(out_rsa(:, 2, 1:125, 1:125)); 
    cond1(ids, :, :) = []; 
    cond2(ids, :, :) = []; 
    diff = cond1-cond2; 
    
    [cond1 cond2] = rem_half_matrix(cond1, cond2);
    
    m1 = squeeze(mean(cond1, 'omitnan')); 
    m2 = squeeze(mean(cond2, 'omitnan')); 
    
    [h p ci ts] = ttest(cond1, cond2); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    %h(1:26,1:2) = 0;  % % % no clusters before baseline
    %h(1:251,1:26) = 0;  % % % no clusters before baseline
    
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        %[max2u id] = max((allSTs));
        tObs = allSTs(id); 
    end
    
    
    %h = zeros(size(cond1, 2),size(cond1, 2)); 
    %h(clustinfo.PixelIdxList{id}) = 1;
    
    
    
    clim = [-.03 .03];
    plot_TG_map(m1, m2, h, t, f2sav, clim)
    exportgraphics(gcf, [paths.results.rsaPlots f2sav '.png'], 'Resolution',150)
    close all 

end


%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)


%% take mean in cluster 

for subji = 1:32

    c1 = squeeze(cond1(subji, :,:));
    c2 = squeeze(cond2(subji, :,:));

    mc1(subji,:) = mean(c1(clustinfo.PixelIdxList{4}));
    mc2(subji,:) = mean(c2(clustinfo.PixelIdxList{4}));

end



%% plot one bar
clc
ylim = [-0.1 0.1];
xlim = [0 3];
 
data.data = [mc1 mc2]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1 2], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 2],'XTickLabel',{'   '},     'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

%% plot 3 bars
clear, clc
load OCC_CSPE-CSMPP-CSMPM
ylim = [-0.1 0.1];
xlim = [0 4];
 
data.data = [mcCSPE mcCSMPP mcCSMPM]; 
figure(2); set(gcf,'Position', [0 0 500 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1:3], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 2);
set(hb,'linestyle','none', 'lineWidth', 2);
set(gca,'XTick',[1 3],'XTickLabel',{'1' '2' '3'},'FontSize', 15, 'linew',2, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 2.5);

%%
clear d4ANOVA
d4ANOVA = [mcCSPE ; mcCSMPP ; mcCSMPM];
d4ANOVA(:,2) = [ones(1,32) ones(1,32)*2 ones(1,32)*3];
d4ANOVA(any(isnan(d4ANOVA), 2), :) = [];
d4ANOVA(:,3) = [1:32 1:32 1:32];

x = RMAOV1(d4ANOVA);
[h p ci ts] = ttest(mcCSPE, mcCSMPP)
%[h p ci ts] = ttest(mcCSPE, mcCSMPM)



%% plot 2 lines from TG 
clear
paths = load_paths_EXT; 
f2sav =   'POW_PFC_C_3-54_0_0_50-1_1_SICSPA-SICSMA';
load ([ paths.results.rsa f2sav '.mat']);

ids = rem_nan_subj_EXT(out_rsa); 
cond1 = squeeze(out_rsa(:, 1, :, :)); 
cond2 = squeeze(out_rsa(:, 2, :, :)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

diff = cond1-cond2; 

for subji = 1:size(cond1, 1)
   cond1B(subji, :) = diag(squeeze(cond1(subji, :, :)));
   cond2B(subji, :) = diag(squeeze(cond2(subji, :, :)));
end       
cond1 = cond1B; cond2 = cond2B; 


d2pm1	= squeeze(mean(cond1,'omitnan'));
d2pm2	= squeeze(mean(cond2,'omitnan'));
d2pstd1	= std(cond1, 'omitnan');
d2pstd2	= std(cond2, 'omitnan');
se1 = d2pstd1/sqrt(size(cond1, 1))
se2 = d2pstd2/sqrt(size(cond1, 1))

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(1, size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;
hb = h; hb(h==0) = nan; hb(hb==1) = -.03; 

times = (-.5:.01:2) + .25;
%times = 1:21
figure(); 
colors2use = brewermap([6],'*Set1')*0.75;
shadedErrorBar(times,  d2pm1, se1, {'Color',colors2use(1,:)}, 1); hold on; 
shadedErrorBar(times, d2pm2, se2,  {'Color',colors2use(2,:)}, 1); hold on; 
plot(times, hb, LineWidth=6)
plot(get(gca,'xlim'), [0 0],'k:', 'linewidth', 1);
plot([0 0],get(gca,'ylim'),'k:', 'linewidth', 1);
%set(gca, 'xlim', [-.25 1.4],'Fontsize', 18);%'ylim', [-.032 .035], 
title(f2sav, 'Interpreter','none')
exportgraphics(gcf, [paths.results.rsa  'myP.png'], 'Resolution',150)

%% permutations 2D (line plot)

nPerm = 1000; 
clear max_clust_sum_perm
for permi = 1:nPerm
    %c1B = squeeze(cond1(:, 51:220)); 
    %c2B = squeeze(cond2(:, 51:220));
    c1B = squeeze(cond1(:, 51:151)); 
    c2B = squeeze(cond2(:, 51:151));
    c1B(c1B == 0) = nan; 
    c2B(c2B == 0) = nan; 
    for subji = 1:size(c1B, 1)
        if rand>.5
           tmp = c1B(subji, :);
           c1B(subji, :) = c2B(subji, :);
           c2B(subji, :) = tmp; 
        end
    end
    
    [hPerm p ci tsPerm] = ttest(c1B, c2B); 
    hPerm = squeeze(hPerm); tPerm = squeeze(tsPerm.tstat);

    clear allSTs  
    clustinfo = bwconncomp(hPerm);
    for pxi = 1:length(clustinfo.PixelIdxList)
        allSTs(pxi,:) = sum(tPerm(clustinfo.PixelIdxList{pxi}));% 
    end
    if exist('allSTs') & ~isempty(clustinfo.PixelIdxList)
        [max2u id] = min(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end
    

end

clear p ratings2u mcsP
ratings2u = tObs; 
mcsP = max_clust_sum_perm;

%allAb = mcsP(mcsP < ratings2u);
allAb = mcsP(abs(mcsP) > abs(ratings2u));
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% PLV 200
clear
paths = load_paths_EXT; 
f2sav =   'PLV_OCC_V_3-54_0_0_200-10_1_SCA-DCA';

load ([ paths.results.rsa f2sav '.mat']);

ids = rem_nan_subj_EXT(out_rsa); 

cond1 = squeeze(out_rsa(:, 1, 6, 6)); 
cond2 = squeeze(out_rsa(:, 2, 6, 6)); 

cond1(ids, :, :) = []; 
cond2(ids, :, :) = []; 

diff = cond1-cond2; 

m1 = squeeze(mean(cond1, 'omitnan')); 
m2 = squeeze(mean(cond2, 'omitnan')); 

[h p ci ts] = ttest(cond1, cond2); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);


allC = [cond1 cond2];
boxplot(allC)

disp (['t: ' num2str(t) ' //  p = ' num2str(p)])

%% CHECK CONTEXT DURING ACQ AND EXT > GENERALIZATION
clear, clc
paths = load_paths_EXT; 

myR = 'PFC';

f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SCA-DCA'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVA-DIDVA'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 
f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SCE-DCE'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVE-DIDVE'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 



% % % % %  remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)
    %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
    %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end
end

% cond1A = squeeze(out_rsa_ACQ(:, 1, 1:13, 1:13)); cond1A(ids, :, :) = []; 
% cond2A = squeeze(out_rsa_ACQ(:, 2, 1:13, 1:13)); cond2A(ids, :,:) = []; 
% cond1E = squeeze(out_rsa_EXT(:, 1, 1:13, 1:13)); cond1E(ids, :, :) = []; 
% cond2E = squeeze(out_rsa_EXT(:, 2, 1:13, 1:13)); cond2E(ids, :, :) = []; 
cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); cond1A(ids, :, :) = []; 
cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); cond2A(ids, :,:) = []; 
cond1E = squeeze(out_rsa_EXT(:, 1, :,:)); cond1E(ids, :, :) = []; 
cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); cond2E(ids, :, :) = []; 

diffA = cond1A-cond2A; 
diffE = cond1E-cond2E; 

[diffA diffE] = rem_half_matrix(diffA, diffE);
m1 = squeeze(mean(diffA, 'omitnan')); 
m2 = squeeze(mean(diffE, 'omitnan')); 

[h p ci ts] = ttest(diffA, diffE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
%h(1:26,1:2) = 0;  % % % no clusters before baseline
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end


%h = zeros(size(cond1, 2),size(cond1, 2)); 
%h(clustinfo.PixelIdxList{id}) = 1;

plot_TG_map(m1, m2, h, t, f2sav, [-.02 .02]); 
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)




%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(diffA, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, diffA(:, 3:20, 3:20), diffE(:, 3:20, 3:20));
%junts = cat(1, diffA(:, 4:13, 4:13), diffE(:, 4:13, 4:13));
%junts = cat(1, diffA(:, 26:125, 26:125), diffE(:, 26:125, 26:125));

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
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));

p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% CHECK CONTEXT (OR ITEM) DURING ACQ AND EXT IN LOOP 
clear, clc
paths = load_paths_EXT; 

listF2sav = {

                 {['POW_TMP_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_TMP_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_OFC_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_OFC_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_HPC_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_HPC_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_AMY_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_AMY_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_PFC_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_PFC_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_PFCO_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_PFCO_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };
                 {['POW_OCC_C_3-54_1_0_50-10_1_DISVA-DIDVA']  ['POW_OCC_C_3-54_1_0_50-10_1_DISVE-DIDVE'] };


%                  {['POW_TMP_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_TMP_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_OFC_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_OFC_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_HPC_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_HPC_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_AMY_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_AMY_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_PFC_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_PFC_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_PFCO_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_PFCO_C_3-54_1_0_50-10_1_SISVE-DISVE'] };
%                  {['POW_OCC_C_3-54_1_0_50-10_1_SISVA-DISVA']  ['POW_OCC_C_3-54_1_0_50-10_1_SISVE-DISVE'] };

%                  {['POW_TMP_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_TMP_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_OFC_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_OFC_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_HPC_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_HPC_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_AMY_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_AMY_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_PFC_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_PFC_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_PFCO_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_PFCO_C_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_OCC_C_3-54_1_0_50-10_1_SCA-DCA']  ['POW_OCC_C_3-54_1_0_50-10_1_SCE-DCE'] };

%                  {['POW_TMP_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_TMP_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_OFC_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_OFC_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_HPC_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_HPC_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_AMY_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_AMY_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_PFC_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_PFC_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_PFCO_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_PFCO_V_3-54_1_0_50-10_1_SCE-DCE'] };
%                  {['POW_OCC_V_3-54_1_0_50-10_1_SCA-DCA']  ['POW_OCC_V_3-54_1_0_50-10_1_SCE-DCE'] };



};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths pM
        
    
    f2sav = listF2sav{listi}{1};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_ACQ = out_rsa; 
    f2sav = listF2sav{listi}{2};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_EXT = out_rsa; 

    
    
    
    % % % % %  remove hack 
    ids = []; 
    for subji = 1:size(out_rsa, 1)
        %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
        %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
        cond1 = squeeze(out_rsa(subji, 1, :, :)); 
        cond2 = squeeze(out_rsa(subji, 2, :, :)); 
        if cond1(1) == 0
            ids = [ids subji];
        end
    end
    
    % cond1A = squeeze(out_rsa_ACQ(:, 1, 1:13, 1:13)); cond1A(ids, :, :) = []; 
    % cond2A = squeeze(out_rsa_ACQ(:, 2, 1:13, 1:13)); cond2A(ids, :,:) = []; 
    % cond1E = squeeze(out_rsa_EXT(:, 1, 1:13, 1:13)); cond1E(ids, :, :) = []; 
    % cond2E = squeeze(out_rsa_EXT(:, 2, 1:13, 1:13)); cond2E(ids, :, :) = []; 
    cond1A = squeeze(out_rsa_ACQ(:, 1, :, :)); cond1A(ids, :, :) = []; 
    cond2A = squeeze(out_rsa_ACQ(:, 2, :, :)); cond2A(ids, :,:) = []; 
    cond1E = squeeze(out_rsa_EXT(:, 1, :,:)); cond1E(ids, :, :) = []; 
    cond2E = squeeze(out_rsa_EXT(:, 2, :, :)); cond2E(ids, :, :) = []; 
    
    diffA = cond1A-cond2A; 
    diffE = cond1E-cond2E; 
    
    [diffA diffE] = rem_half_matrix(diffA, diffE);
    m1 = squeeze(mean(diffA, 'omitnan')); 
    m2 = squeeze(mean(diffE, 'omitnan')); 
    
    [h p ci ts] = ttest(diffA, diffE); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        tObs = allSTs(id); 
    end
    
    
    %h = zeros(size(cond1, 2),size(cond1, 2)); 
    %h(clustinfo.PixelIdxList{id}) = 1;
    
    plot_TG_map(m1, m2, h, t, f2sav, [-.02 .02]); 
    exportgraphics(gcf, [paths.results.rsa  '_'  listF2sav{listi}{1}(1:9) '_.png'], 'Resolution',150)

    close all; 


end




%% CHECK VALENCE DURING ACQ AND EXT 
clear, clc
paths = load_paths_EXT; 

myR = 'PFCO';
pM = 2; 

f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SCA-DCA'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SICSPA-SICSMA'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVA-DIDVA'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_ACQ = out_rsa; 
f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SCE-DCE'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_SICSPE-SICSME'];
%f2sav =   ['POW_' myR '_C_3-54_1_0_50-10_1_DISVE-DIDVE'];
load ([ paths.results.rsa f2sav '.mat']);
out_rsa_EXT = out_rsa; 




% % % % %  remove hack 
ids = []; 
for subji = 1:size(out_rsa, 1)
    %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
    %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
    cond1 = squeeze(out_rsa(subji, 1, :, :)); 
    cond2 = squeeze(out_rsa(subji, 2, :, :)); 
    if cond1(1) == 0
        ids = [ids subji];
    end
end

condA = squeeze(out_rsa_ACQ(:, pM, 1:19, 1:19)); condA(ids, :, :) = []; 
condE = squeeze(out_rsa_EXT(:, pM, 1:19,1:19)); condE(ids, :, :) = []; 

[condA condE] = rem_half_matrix(condA, condE);
m1 = squeeze(mean(condA, 'omitnan')); 
m2 = squeeze(mean(condE, 'omitnan')); 

[h p ci ts] = ttest(condA, condE); 
h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
clustinfo = bwconncomp(h);
for pxi = 1:length(clustinfo.PixelIdxList)
   allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
end

if exist('allSTs')
    [max2u id] = max(abs(allSTs));
    tObs = allSTs(id); 
end

if exist('id')
    h = zeros(size(condA, 2),size(condA, 2)); 
    h(clustinfo.PixelIdxList{id}) = 1;
end

plot_TG_map(m1, m2, h, t, f2sav, [-.02 .02]); 
exportgraphics(gcf, [paths.results.rsa  '_myP.png'], 'Resolution',150)



%% PERMUTATIONS
nPerm = 1000; 

nSubj =  size(condA, 1);
realCondMapping = [zeros(1,nSubj); ones(1, nSubj)]';

junts = cat(1, condA(:, 3:20, 3:20), condE(:, 3:20, 3:20));

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
    [h p ci ts] = ttest(diffC); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat); 
    clear allSTs  
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(allSTs);
        max_clust_sum_perm(permi,:) = allSTs(id); 
    else
        max_clust_sum_perm(permi,:) = 0; 
    end

end


disp('done')

 
%allAb = max_clust_sum_perm(max_clust_sum_perm < tObs);
allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > abs(tObs));
%allAb = max_clust_sum_perm(abs(max_clust_sum_perm) > 36.9304012965767);
p = 1 - ((nPerm-1) - (length (allAb)))  / nPerm


%% CHECK VALENCE DURING ACQ AND EXT IN LOOP 
clear, clc
paths = load_paths_EXT; 

pM = 1; 


listF2sav = {
                 {['POW_TMP_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_TMP_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_OFC_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_OFC_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_HPC_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_HPC_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_AMY_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_AMY_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_PFC_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_PFC_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_PFCO_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_PFCO_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };
                 {['POW_OCC_C_3-54_1_0_50-10_1_SICSPA-SICSMA']  ['POW_OCC_C_3-54_1_0_50-10_1_SICSPE-SICSME'] };

};   

paths = load_paths_EXT; 
t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1 paths pM
        
    
    f2sav = listF2sav{listi}{1};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_ACQ = out_rsa; 
    f2sav = listF2sav{listi}{2};
    load ([ paths.results.rsa f2sav '.mat']);
    out_rsa_EXT = out_rsa; 
    
    
    % % % % %  remove hack 
    ids = []; 
    for subji = 1:size(out_rsa, 1)
        %cond1 = squeeze(out_rsa(subji, 1, 1:13, 1:13)); 
        %cond2 = squeeze(out_rsa(subji, 2, 1:13, 1:13)); 
        cond1 = squeeze(out_rsa(subji, 1, :, :)); 
        cond2 = squeeze(out_rsa(subji, 2, :, :)); 
        if cond1(1) == 0
            ids = [ids subji];
        end
    end
    
    condA = squeeze(out_rsa_ACQ(:, pM, :, :)); condA(ids, :, :) = []; 
    condE = squeeze(out_rsa_EXT(:, pM, :,:)); condE(ids, :, :) = []; 
    
    [condA condE] = rem_half_matrix(condA, condE);
    m1 = squeeze(mean(condA, 'omitnan')); 
    m2 = squeeze(mean(condE, 'omitnan')); 
    
    [h p ci ts] = ttest(condA, condE); 
    h = squeeze(h); h(isnan(h)) = 0; t = squeeze(ts.tstat);
    clustinfo = bwconncomp(h);
    for pxi = 1:length(clustinfo.PixelIdxList)
       allSTs(pxi) = sum(t(clustinfo.PixelIdxList{pxi}));% 
    end
    
    if exist('allSTs')
        [max2u id] = max(abs(allSTs));
        tObs = allSTs(id); 
    end
    
    
    %h = zeros(size(cond1, 2),size(cond1, 2)); 
    %h(clustinfo.PixelIdxList{id}) = 1;
    
    plot_TG_map(m1, m2, h, t, f2sav, [-.02 .02]); 
    exportgraphics(gcf, [paths.results.rsa  '_'  listF2sav{listi}{1}(1:9) '_' num2str(pM) '_.png'], 'Resolution',150)
    
    close all 


end





%% plot histogram
figure
%tObs =  -30.4546%-86.4470;
histogram(max_clust_sum_perm, 20); hold on; 
scatter(tObs,0, 100, 'filled','r');
set(gca, 'FontSize', 16)


%% COMPUTE SUCCESIVE TRIALS SIMILARITY 
%rsaTYPE_freqs_avTimeFeatVect_freqResolv(0-1)_win-width_mf_TG_contrast
clear , clc

listF2sav = {
                'POW_HPC_C_3-54_1_0_50-1_1_SCA-STR';
                               
        };   

t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
        
        EEG = ALLEEG{subji};
        
        
        if ~isempty(EEG)

            
            EEG = add_EEGLAB_fields(EEG); 
            EEG = rem_nan_trials_EXT(EEG); %can't remove nan trials here or the order is lost

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = Ev2; 

            EEG = extract_power_EXT(EEG, 0.01); 
            %EEG = normalize_baseline_EXT(EEG, [251:300]); 
            EEG = normalize_EXT(EEG);  %across trials
            cfg.oneListPow = EEG.power(:, :, : ,251:550); 
            %out_contrasts = create_contrasts_trials_EXT(cfg);
            out_contrasts = create_contrasts_EXT(cfg);
            tic
            out_rsa{subji,:} = rsa_EXT6(out_contrasts, cfg);
            toc
        
            ids{subji,:} = out_contrasts.allIDs;
            nnans{subji} = EEG.nan; 
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'ids');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end


%% plot 

out_rsa1 = out_rsa(~cellfun('isempty', (out_rsa)));
id2sav1 = ids(~cellfun('isempty', (ids)));


for subji = 1:length(out_rsa1)

    rsah = out_rsa1{subji};
    idh = id2sav1{subji}{1}; 
    idC = mod(idh(:, 3), 10); 

    c1 = mean(rsah(idC == 1, toi), 2); 
    c2 = mean(rsah(idC == 2, toi), 2); 
    c3 = mean(rsah(idC == 3, toi), 2); 
    c4 = mean(rsah(idC == 4, toi), 2); 

    allC{subji,:} = padcat(c1, c2, c3, c4);



end

%% COMPUTE SUCCESIVE TRIALS SIMILARITY AND RATINGS CORRELATIONS
clear, clc
paths = load_paths_EXT; 

f2sav =  'POW_OCC_C_3-54_0_0_50-1_1_ALLAE-STR';
%f2sav =  'POW_OFC_C_3-54_0_0_50-1_1_CSME-STR';

load ([ paths.results.rsa f2sav '.mat']);

load clustinfo_OCC_SICSPE-SICSME_50-1
%load clustinfo_OFC_DISVA-DIDVA_50-1

out_rsa = out_rsa(~cellfun('isempty', out_rsa));
ids = ids(~cellfun('isempty', ids));

for subji = 1:length(out_rsa)


    idsH = ids{subji}{1};
    %ratings = idsH(:, 7); 
    %ratings = idsH(:, 17) ;
    ratings = abs ( idsH(:, 7) - idsH(:,17) );
    %ratings = mean([idsH(:, 7) idsH(:,17)], 2); 
    

    orS = out_rsa{subji}; 
    %orS2 = squeeze(mean(mean(orS(:,:,3:23, 3:23), 4), 3));
    orS2 = squeeze(mean(mean(orS(:,:,26:125, 26:125), 4), 3))';
    
%     clear orS2
%     for triali = 1:size(orS, 2)
%         orS2a = squeeze(orS(:,triali,:,:));
%         orS2(triali,:) = mean(orS2a(clustinfo.PixelIdxList{11}), 'omitnan');
%         
% 
%         %orS2a = squeeze(orS(:,triali,:,:));
%         %dv = diag(orS2a); dv = dv(26:100);
%         %orS2(triali,:) = mean(dv);
%     end

    idM = isnan(ratings);
    idN = isnan(orS2); 
    orS2(idN | idM) = []; 
    ratings(idN | idM) = []; 
    allRS(subji, :) = corr(orS2, ratings, 'type', 'k' ); 
        
% % %         figure()
% % %         %plot(ratings, orS2); hold on; 
% % %         scatter(orS2, ratings); hold on; 
% % %         h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% % %         C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% % %         allSlopes(subji, :) = C(2);
% % %         allIntercepts(subji, :) = C(1);
% % %         %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
% % %         

% %     figure()
% %     plot(orS2); hold on; 
% %     scatter(1:length(orS2), orS2); hold on; 
% %     h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
% %     C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
% %     allSlopes(subji, :) = C(2);
% %     allIntercepts(subji, :) = C(1);
% %     %set(gca, 'ylim', [1 4], 'xlim', [-2 2], 'Fontsize', 24)
% %     close all




end

%boxplot(allSlopes)
%[h p ci ts] = ttest(allSlopes); 

figure()
boxplot(allRS)
[h p ci ts] = ttest(allRS); 
disp(['T > ' num2str(ts.tstat) '  P > ' num2str(p)])



%% COMPUTE SIMILARITY OF EACH TRIAL TO ALL OTHER TRIALS
clear , clc

listF2sav = {

'POW_PFC_V_3-54_1_0_50-10_1_ALLE-TR';


};   

t1 = datetime; 
for listi = 1:length(listF2sav)
    clearvars -except listF2sav listi t1
        
    f2sav       = listF2sav{listi}; 
    cfg = getParams_EXT(f2sav);


    paths = load_paths_EXT; 
    
    ALLEEG = loadTracesEXT(cfg.roi, cfg.LT, paths); %LT = locked to
    
    
    for subji = 1:length(ALLEEG)
        
        EEG = ALLEEG{subji};

        if ~isempty(EEG)
            disp(['File > ' num2str(listi) '      ' listF2sav{listi} '    Subject > ' num2str(subji)]);
            
            EEG = add_EEGLAB_fields(EEG); 
            %EEG = rem_nan_trials_EXT(EEG); %can't remove nan trials here or the order is lost

            Ev = [{EEG.event.type}]';Ev1 = cellfun(@(x) strsplit(x, '_'), Ev, 'un', 0); 
            Ev2 = cat(1, Ev1{:});
            cfg.oneListIds = getIds_EXT(Ev2, cfg); 
            id2sav{subji, :} = double(string(Ev2(cfg.oneListIds, :)));

            EEG = extract_power_EXT(EEG, 0.01); 
            EEG = normalize_EXT(EEG);  %across trials
            cfg.oneListPow = EEG.power(cfg.oneListIds, :, : ,251:550); 

            neuralRDM = computeNeuralRDMs_EXT(cfg); 
            
            clear trialRDM
            for triali = 1:size(neuralRDM, 1)
                rowOfRDM = squeeze(neuralRDM(triali,:,:)); 
                rowOfRDM(triali,:) = []; 
                trialRDM(triali, :) = mean(rowOfRDM, 'omitnan');
            end
            out_rsa{subji,:} = trialRDM; 
    
        end
        
    end

    mkdir ([paths.results.rsa]);
    save([ paths.results.rsa f2sav '.mat'], 'out_rsa', 'id2sav');
    
    
    t2 = datetime; 
    etime(datevec(t2), datevec(t1))

end

%% plot 

clear , clc
paths = load_paths_EXT; 

load ([paths.results.rsa '/POW_OFC_C_3-54_1_0_50-10_1_ALLE-TR'])
toi = 3:20; 

out_rsa1 = out_rsa(~cellfun('isempty', (out_rsa)));
id2sav1 = id2sav(~cellfun('isempty', (id2sav)));


for subji = 1:length(out_rsa1)

    rsah = out_rsa1{subji};
    idh = id2sav1{subji}; 
    idC = mod(idh(:, 3), 10); 
    idCS = idh(:, 8); 

%     c1 = mean(rsah(idC == 1, toi), 2) - mean(rsah(idC ~= 1, toi), 2); 
%     c2 = mean(rsah(idC == 2, toi), 2)- mean(rsah(idC ~= 2, toi), 2);  
%     c3 = mean(rsah(idC == 3, toi), 2)- mean(rsah(idC ~= 3, toi), 2); 
%     c4 = mean(rsah(idC == 4, toi), 2)- mean(rsah(idC ~= 4, toi), 2); 

    for triali = 1:size(rsah, 1)
        idSC = idh(:,3) == idh(triali, 3); 
        idDC = idh(:,3) ~= idh(triali, 3); 
        rsat = rsah; 
        rsat(triali,:) = [];
        idSC(triali,:) = []; 
        idDC(triali,:) = []; 
        allCs(triali, :) = mean(rsat(idSC, toi), 'all', 'omitnan') - mean(rsat(idDC, toi), 'all', 'omitnan'); 

    end
    

    c1 = allCs(idC == 1,:); 
    c2 = allCs(idC == 2,:); 
    c3 = allCs(idC == 3,:); 
    c4 = allCs(idC == 4,:); 
    allC{subji,:} = [c1 c2 c3 c4];


end


M = max(cellfun(@length, allC));
allC3 = cellfun(@(x) [x; nan(M - numel(x), 1)], allC, 'un', 0);
allC4 = [allC3{:}]'; 


%figure
%plot(mean(allC4))
%plot(mean(allC4, 'omitnan'))

 

[aov, tbl, stats] = anova1(allC4);

multcompare(stats);


%% Correlate context specific activity and AMY theta power
clear , clc
paths = load_paths_EXT; 

toi  = 3:20;

load ([paths.results.rsa 'POW_AMY_C_3-54_1_0_50-10_1_ALLE-TR'])
load allPOWAMY


ids1 = ~cellfun('isempty', out_rsa); ids1(end+1:50) = 0; 
ids2 = ~cellfun('isempty', allPOWAMY);
idsBoth = ids1 & ids2; 

out_rsa1 = out_rsa(idsBoth);
id2sav1 = id2sav(idsBoth); 
allPOWAMY1 = allPOWAMY(idsBoth); 


for subji = 1:length(out_rsa1)
    
    amyPowV = allPOWAMY1{subji}; 
    ctxPFC = mean(out_rsa1{subji}(:, toi), 2);
    idh = id2sav1{subji}; 

    idNNaN = isnan(amyPowV) | isnan(ctxPFC); 

    amyPowV(idNNaN) = []; 
    ctxPFC(idNNaN) = []; 
    idh(idNNaN,:) = []; 
    
    %amyPowV = amyPowV(idh(:, 6) == 1); 
    %ctxPFC = ctxPFC(idh(:, 6) == 1); 

    % compute same minus different 
    for triali = 1:size(idh, 1)
        idtrSC = idh(:, 3) == idh(triali, 3); 
        idtrDC = idh(:, 3) ~= idh(triali, 3); 
        idtrSC(triali) = []; 
        idtrDC(triali) = []; 
        tr_AMY(triali, :) = amyPowV(triali); 
        %tr_AMY(triali, :) = mean(amyPowV(idtrSC)) - mean(amyPowV(idtrDC)) ; 
        tr_PFC(triali, :) = mean(ctxPFC(idtrSC)) - mean(ctxPFC(idtrDC)) ; 

    end

    %[B, tF1] = rmoutliers(tr_AMY, 'percentiles', [10 90]); 
    %[B, tF2] = rmoutliers(tr_PFC, 'percentiles', [10 90]); 
%     [B, tF1] = rmoutliers(tr_AMY); 
%     [B, tF2] = rmoutliers(tr_PFC); 
%     tF3 = tF1 | tF2; 
%     tr_AMY(tF3) = []; 
%     tr_PFC(tF3) = []; 



    allRHO(subji, :) = corr(tr_AMY, tr_PFC, 'type', 's'); 
    %allRHO(subji, :) = corr(amyPowV, ctxPFC, 'type', 's'); 
% % 
    figure()
    scatter(tr_AMY, tr_PFC, 150, 'filled');
    %scatter(amyPowV, ctxPFC, 150, 'filled');
    h2 = lsline;h2.LineWidth = 2;h2.Color = [.5 .5 .5 ];
    C = [ones(size(h2.XData(:))), h2.XData(:)]\h2.YData(:);
    allSlopes(subji, :) = C(2);
    allIntercepts(subji, :) = C(1);
    set(gca, 'Fontsize', 24)

end


%%plot one bar
clc
data.data = [allRHO]; 

ylim = [-0.85 0.85];
xlim = [0 2];
 

figure(2); set(gcf,'Position', [0 0 300 600]); 
mean_S = mean(data.data, 1);
std_S = std(data.data, [], 1);
h = bar (mean_S);hold on;
hb = plot ([1], data.data); hold on; % > lines
set(hb, 'lineWidth', 1, 'Marker', '.', 'MarkerSize',30);hold on;
set(h,'FaceColor', 'none', 'lineWidth', 1.5);
set(hb,'linestyle','none', 'lineWidth', 1.5);
set(gca,'XTick',[1],'XTickLabel',{'   '}, 'FontSize', 18, 'linew',1.5, 'ylim', ylim, 'xlim', xlim);
plot(get(gca,'xlim'), [0 0],'k','lineWidth', 1.5);

[h p ci ts] = ttest(allRHO); 
disp (['t(' num2str(ts.df) ') = ' num2str(ts.tstat, 3) ',' ' p = ' num2str(p, 3)]);


exportgraphics(gcf, ['_myP.png'], 'Resolution',150)










%% 
























%%
















