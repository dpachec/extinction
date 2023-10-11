function [allRSAZ] =  rsa_WM (out_contrasts, cfg)


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
        elseif strcmp(cfg.tyRSA, 'TR')
            nTimepoints = size (out_contrasts.allContrasts{i}{end}, 4); %%all2all file is stored within this cell array
            aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
        end
        
    end
end

idEmpty = cell2mat(cellfun(@(x) isempty(x), currentContrast, 'un', 0));
currentContrast(idEmpty) = []; 
for coni = 1:length(currentContrast)
    bins = aBins(coni);
    allIDs = out_contrasts.allIDs{coni};
    id = currentIds{coni};
    clear allRSA %critical to not merge conditions
    for batchi = 1:length(currentContrast{coni}) % for every batch
        %batchi
        if ~isempty(currentContrast{coni}) 
            all2all = currentContrast{coni}{batchi};
        end
    
        trialN = size(all2all, 1);    
        chanN = size(all2all, 3);
    
    
        %disp (['Cond ' id '   ' num2str(size(all2all, 1)) ' trials']);
    
        if strcmp(cfg.tyRSA, 'POW')
            if aVTime
                    xM = zeros (trialN, bins,  chanN * length(f));
                    yM = zeros (trialN, bins,  chanN * length(f));
            else 
                    xM = zeros (trialN, bins,  chanN * length(f) * win_width);
                    yM = zeros (trialN, bins,  chanN * length(f) * win_width);
            end
                    
            parfor timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                if aVTime
                    x = mean(all2all(:, 1,:,f, timeBins), 5, 'omitnan');
                    x = reshape (x, [trialN, chanN * length(f)]);
                    y = mean(all2all(:, 2,:,f,timeBins), 5, 'omitnan');
                    y = reshape (y, [trialN, chanN * length(f)]);
                
                else
                    x = all2all(:, 1,:,f,timeBins);
                    x = reshape (x, [trialN, chanN * length(f)* win_width]);
    
                    y = all2all(:, 2,:,f,timeBins);
                    y = reshape (y, [trialN, chanN * length(f)* win_width]);
                end
        
                xM(:, timei, :) =  x;
                yM(:, timei, :) =  y;
                %disp(['size xM >>    ' num2str(size(xM))])
                
            end

        elseif strcmp(cfg.tyRSA, 'TR') % temporal RSA

                xM = zeros (trialN, bins,  chanN *  win_width);
                yM = zeros (trialN, bins,  chanN *  win_width);
                
            parfor timei = 1:bins 
                %timeBins(timei,:) = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                timeBins = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1;
                x = all2all(:, 1,:,timeBins);
                x = reshape (x, [trialN, chanN * win_width]);

                y = all2all(:, 2,:,timeBins);
                y = reshape (y, [trialN, chanN * win_width]);
            
                xM(:, timei, :) =  x;
                yM(:, timei, :) =  y;
                %disp(['size xM >>    ' num2str(size(xM))])
                
            end

        end
    
        
    
        rsaZ = zeros (trialN, bins, bins);
        for triali = 1:trialN
            %if isempty(find(isnan(xM(triali, :,:)))) & isempty(find(isnan(yM(triali, :,:)))) & isempty(find(xM(triali, :,:)==0)) & isempty(find(yM(triali, :,:)==0))
                mX= squeeze(xM(triali,:,:));
                mY= squeeze(yM(triali,:,:));
                r = corr (mX', mY','Type', 's'); 
                
                rsaZ(triali, :, :) = atanh(r);
            %else
            %    rsaZ(triali, :, :) = nan(size(bins));
            %end
        end
        
        %rsaZ(rsaZ==0) = nan; 
        rsaZ(isinf(rsaZ)) = nan;
      
        allRSA{batchi} = rsaZ; 
     
        
    end
    
 
    if exist('allRSA')
            rsaZ = cat(1, allRSA{:});
    end


    if TG==1
        if ~isempty(rsaZ)
            allRSAZ(coni, :, :) = squeeze(mean(rsaZ, 'omitnan')); 
        else 
            if cfg.fR
                allRSAZ(coni, :, :) = nan (nFreq, bins);
            end
            if cfg.TG
                allRSAZ(coni, :, :) = nan (bins);
            end        
        end
        
        
    else %only store the diagonal 
        if ~isempty(rsaZ)
            rsaZ = cat(1, allRSA{:});
            parfor triali = 1:size(rsaZ, 1)
                rsaN(triali, :) = diag(squeeze(rsaZ(triali, :, :)));
            end
            allRSAZ(coni, :) = squeeze(mean(rsaN, 'omitnan')); 
        else 
            allRSAZ(coni, :) = nan (1, bins);
        end
    end




    end 
end


 
 
 
 

