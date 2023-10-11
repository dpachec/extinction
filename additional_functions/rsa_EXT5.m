function [allRSAZ] =  rsa_EXT5 (out_contrasts, cfg)


win_width = cfg.win_width; 
mf = cfg.mf; 
f = cfg.freqs; 
TG = cfg.TG; 
aVTime = cfg.avTFV; 
    

currentContrast = out_contrasts.allContrasts;
currentIds = out_contrasts.allContrastIds;


for i = 1:length(currentContrast)
    if ~isempty(out_contrasts.allContrasts{i})
        if strcmp(cfg.tyRSA, 'POW')
            nTimepoints = size (out_contrasts.allContrasts{i}{end}, 5); %%all2all file is stored within this cell array
            aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
        elseif strcmp(cfg.tyRSA, 'TR')| strcmp(cfg.tyRSA, 'PHA')  | strcmp(cfg.tyRSA, 'PLV') 
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
        %batchi
        clear rsaZ
        if ~isempty(currentContrast{coni}) 
            all2all = currentContrast{coni}{batchi};
        end
    
        nTrials = size(all2all, 1);    
        nChans = size(all2all, 3);
        

        if cfg.TG 
            for chani = 1:nChans
                for triali = 1:nTrials
                    for timei = 1:bins 
                        timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                        x = all2all(triali, 1,chani,timeBinsi);
                        x = x(:); 
                        for timej = timei:bins
                            timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                            y = all2all(triali, 2,chani,timeBinsj);
                            y = y(:); 
                            
                            diffPha = angle(hilbert(x)) - angle(hilbert(y));
                            PLV2U = abs(mean(exp(1i*(diffPha))));
                            
                            rsaZ(chani, triali, timei, timej) = PLV2U;
                        end
                    end
                end
            end

        else
            for chani = 1:nChans
                for triali = 1:nTrials
                    for timei = 1:bins 
                        timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                    
                        x = all2all(triali, 1,chani,timeBinsi);
                        x = x(:); 
        
                        y = all2all(triali, 2,chani,timeBinsi);
                        y = y(:); 
                        
                        diffPha = angle(hilbert(x)) - angle(hilbert(y));
                        PLV2U = abs(mean(exp(1i*(diffPha))));
                        
                        rsaZ(chani, triali, timei) = PLV2U;
    
                    end
                end
            end
        end
            
        %rsaZ(rsaZ==1) = nan; 
        %rsaZ = squeeze(mean(rsaZ, 'omitnan'));
        if exist('rsaZ')
            rsaZ = mean(rsaZ, 'omitnan');
            allRSA{batchi} = rsaZ; 
        end
        
     
        
    end
    
 
 

 
        if exist('rsaZ')
            rsaZ = cat(2, allRSA{:});
            if cfg.TG
                allRSAZ(coni, :, :) = squeeze(mean(rsaZ, 2,'omitnan')); 
            end
        else
            if cfg.TG
                allRSAZ(coni, :, :) = nan(bins); 
            end
        end
       

    end 
end


 
 
 
 

