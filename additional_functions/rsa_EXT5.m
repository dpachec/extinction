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
        %batchi
        clear rsaZ
        if ~isempty(currentContrast{coni}) 
            all2all = currentContrast{coni}{batchi};
        end
    
        nTrials = size(all2all, 1);    
        nChans = size(all2all, 3);
        nFreq = size(all2all, 4);
    
   

        for chani = 1:nChans
            for triali = 1:nTrials
                if isempty(find(isnan(all2all(triali, :,chani,:,:))))
                    for freqi = 1:nFreq
                        for timei = 1:bins 
                            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                        
                            x = all2all(triali, 1,chani,freqi,timeBinsi);
                            x = x(:); 
            
                            y = all2all(triali, 2,chani,freqi,timeBinsi);
                            y = y(:); 
                            
                            diffPha = angle(hilbert(x)) - angle(hilbert(y));
                            PLV2U = abs(mean(exp(1i*(diffPha))));
                            
                            rsaZ(chani, triali, freqi, timei) = PLV2U;
                            
                        end
                    end
                else
                        rsaZ(chani, triali, :, :) = nan(nFreq,bins);
                end
            end
        end
            
            % % normalize
            rsaZ = normalize_rsaZ_EXT(rsaZ, 1:5);

            rsaZ(rsaZ==1) = nan; 
            %rsaZ = squeeze(mean(rsaZ, 'omitnan'));
            rsaZ = mean(rsaZ, 'omitnan');
            allRSA{batchi} = rsaZ; 
        
     
        
    end
    
 
 


   
    
        rsaZ = squeeze(cat(2, allRSA{:}));
        % count nan trials 
        count = 0; 
        for triali = 1:size(rsaZ, 1)
            if isnan(rsaZ(triali, 10, 10))
                count = count+1; 
            end
        end
        disp (['number of nan trials = ' num2str(count)])
        allRSAZ(coni, :, :) = squeeze(mean(rsaZ, 'omitnan')); 


       

    end 
end


 
 
 
 

