function [allRSAZ] =  rsa_WM4 (out_contrasts, cfg)


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
        elseif strcmp(cfg.tyRSA, 'TR')
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
    
        trialN = size(all2all, 1);    
        chanN = size(all2all, 3);
        nFreq = size(all2all, 4);
    
   

        if cfg.fR
            for triali = 1:trialN
                if isempty(find(isnan(all2all(triali, :,:,:,:))))
                    for freqi = 1:nFreq
                        for timei = 1:bins 
                            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                        
                            x = all2all(triali, 1,:,freqi,timeBinsi);
                            x = x(:); 
            
                            y = all2all(triali, 2,:,freqi,timeBinsi);
                            y = y(:); 
                            
                            %r = corr(x, y, 'type', 's');
        
                            a = [x y]';
                            n = length(a);
                            a1 = tiedrank(a')'; 
                            r = 1-6*sum((a1(1,:)-a1(2,:)).^2)/(n*(n^2-1));
                            
                            rsaZ(triali, freqi, timei) = atanh(r);
                            
                        end
                    end
                else
                     rsaZ = nan(trialN, nTimepoints, nTimepoints); 
                end
            end
        else
            if cfg.TG
                parfor triali = 1:trialN
                    if isempty(find(isnan(all2all(triali, :,:,:,:))))
                        for timei = 1:bins 
                            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                            x = all2all(triali, 1,:,f,timeBinsi);
                            x = x(:); 
                            for timej = 1:bins
                                timeBinsj = (timej*mf) - (mf-1):(timej*mf - (mf-1) )+win_width-1;
                                y = all2all(triali, 2,:,f,timeBinsj);
                                y = y(:); 
                                
                                %r = corr(x, y, 'type', 's');
            
                                a = [x y]';
                                n = length(a);
                                a1 = tiedrank(a')'; 
                                r = 1-6*sum((a1(1,:)-a1(2,:)).^2)/(n*(n^2-1));
                                
                                rsaZ(triali, timei, timej) = atanh(r);
                            end
                        end
                    end
                end

            else
                for triali = 1:trialN
                    if isempty(find(isnan(all2all(triali, :,:,:,:))))
                        for timei = 1:bins 
                            timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                            x = all2all(triali, 1,:,f,timeBinsi);
                            x = x(:); 
                            y = all2all(triali, 2,:,f,timeBinsi);
                            y = y(:); 
                            a = [x y]';
                            n = length(a);
                            a1 = tiedrank(a')'; % tiedrank accepts matrix input, but make sure the matrix is in the correct orientation!!
                            r = 1-6*sum((a1(1,:)-a1(2,:)).^2)/(n*(n^2-1));                           
                            rsaZ(triali, timei) = atanh(r);
                        end
                    end
                end
            end
        end
    
        if exist('rsaZ')
            rsaZ(rsaZ==0) = nan; 
            rsaZ(isinf(rsaZ)) = nan;
            allRSA{batchi} = rsaZ; 
        end
     
        
    end
    
 
 


   
        if exist('allRSA')
            rsaZ = cat(1, allRSA{:});
            % count nan trials 
            count = 0; 
            for triali = 1:size(rsaZ, 1)
                if isnan(rsaZ(triali, 1, 1))
                    count = count+1; 
                end
            end
            disp (['number of nan trials = ' num2str(count)])
            if ndims(rsaZ) == 2
                allRSAZ(coni, :) = squeeze(mean(rsaZ, 'omitnan')); 
            elseif ndims(rsaZ) == 3
                allRSAZ(coni, :, :) = squeeze(mean(rsaZ, 'omitnan')); 
            end
        else 
            if cfg.fR
                allRSAZ(coni, :, :) = nan (nFreq, bins);
            end
            if cfg.TG
                allRSAZ(coni, :, :) = nan (bins, bins);
            end
        end
        

       

    end 
end


 
 
 
 

