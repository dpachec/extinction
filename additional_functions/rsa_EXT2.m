function [allRSAZ] =  rsa_WM2 (out_contrasts, cfg)


win_width = cfg.win_width; 
mf = cfg.mf; 
f = cfg.freqs; 
TG = cfg.TG; 
aVTime = cfg.avTFV; 
    

currentContrast = out_contrasts.allContrasts;
currentIds = out_contrasts.allContrastIds;


for i = 1:length(currentContrast)
    if ~isempty(out_contrasts.allContrasts{i})
        if strcmp(cfg.tyRSA, 'pRSA')
            nTimepoints = size (out_contrasts.allContrasts{i}{end}, 5); %%all2all file is stored within this cell array
            aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
        elseif strcmp(cfg.tyRSA, 'tRSA')
            nTimepoints = size (out_contrasts.allContrasts{i}{end}, 4); %%all2all file is stored within this cell array
            aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
        end
        
    end
end

for coni = 1:length(currentContrast)
    bins = aBins(coni);
    allIDs = out_contrasts.allIDs{coni};
    id = currentIds{coni};
    clear allRSA  %critical to not merge conditions
    for batchi = 1:length(currentContrast{coni}) % for every batch
        clear rsaZ
        %batchi
        if ~isempty(currentContrast{coni}) 
            all2all = currentContrast{coni}{batchi};
        end
    
        trialN = size(all2all, 1);    
        chanN = size(all2all, 3);
    
    
   

    
        for triali = 1:trialN
            if isempty(find(isnan(all2all(triali, :,:,:,:))))
                for timei = 1:bins 
                    timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                    for timej = 1:bins 
                        timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                    
                        x = all2all(triali, 1,:,f,timeBinsi);
                        x = x(:);
        
                        y = all2all(triali, 2,:,f,timeBinsj);
                        y = y(:);
                        
                        r = circ_corrcc(x, y);
                        rsaZ(triali, timei, timej) = atanh(r);
                    
                    end
                end
            else
                rsaZ(triali, :, :) = nan(18);
            end
        end
    
        %rsaZ(rsaZ==0) = nan; 
        rsaZ(isinf(rsaZ)) = nan;
      
        allRSA{batchi} = rsaZ; 
     
        
    end
    
 
 


    if TG==1
        %filename = ['s' num2str(sessi, '%02.f') '_' id '_gOBO'   '_rsa.mat'];
        rsaZ = cat(1, allRSA{:});
        % count nan trials 
        count = 0; 
        for triali = 1:size(rsaZ, 1)
            if isnan(rsaZ(triali, 1, 1))
                count = count+1; 
            end
        end
        disp (['number of nan trials = ' num2str(count)])
        allRSAZ(coni, :, :) = squeeze(mean(rsaZ, 'omitnan')); 
        %allRSAZ{coni} = rsaZ; 
        
    else %only store the diagonal 
        rsaZ = cat(1, allRSA{:});
        parfor triali = 1:size(rsaZ, 1)
            rsaN(triali, :) = diag(squeeze(rsaZ(triali, :, :)));
        end
        rsaZ = rsaN;
        filename = ['s' num2str(sessi, '%02.f') '_' id '_dOBO'   '_rsa.mat'];
        save (filename, 'rsaZ', 'allIDs'); %, 'timeBins'
    end




    end 
end


 
 
 
 
