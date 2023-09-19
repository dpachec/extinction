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
        nTimepoints = size (out_contrasts.allContrasts{i}{end}, 5); %%all2all file is stored within this cell array
        aBins(i,:)  =  floor ( (nTimepoints/mf)- win_width/mf+1 );
    end
end

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
    
    
        if aVTime
                xM = zeros (trialN, bins,  chanN * length(f));
                yM = zeros (trialN, bins,  chanN * length(f));
        else 
                xM = zeros (trialN, bins,  chanN * length(f) * win_width);
                yM = zeros (trialN, bins,  chanN * length(f) * win_width);
        end
                
        for timei = 1:bins 
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
    
        
    
        rsaZ = zeros (trialN, bins, bins);
        %fprintf('\n'); fprintf('trial correlation:          '); 
        for triali = 1:trialN
            mX= squeeze(xM(triali,:,:));
            mY= squeeze(yM(triali,:,:));
            r = corr (mX', mY','Type', 's'); 
            %id0 = find(r==0);
            %r = tril(squeeze(r)); %symmetric so only half is saved
            %r(r == 0) = nan;r(id0)=0;
            rsaZ(triali, :, :) = atanh(r);
        end
        
        rsaZ(isinf(rsaZ)) = nan;
      
        allRSA{batchi} = rsaZ; 
     
        
    end
    
 
 


    if TG==1
        %filename = ['s' num2str(sessi, '%02.f') '_' id '_gOBO'   '_rsa.mat'];
        rsaZ = cat(1, allRSA{:});
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


 
 
 
 

