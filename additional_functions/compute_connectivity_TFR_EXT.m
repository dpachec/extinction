function [CON] = compute_connectivity_EXT(data, toi, toiB, foi, ids2comp)
        clear PLV_M1 PLV_M2 PLI_M1 PLI_M2 WPLI COH PLV_M1Bas PLV_M2Bas PLI_M1Bas PLI_M2Bas WPLIBas COHBas
        nTrials = size(data, 3); 
            for combi = 1:length(ids2comp(:, 1))
                
                disp(['Comb: ' num2str(combi)  ])
                 for triali = 1:nTrials
                    
                    amyD = squeeze(data(ids2comp(combi, 1), :, triali));
                    hpcD = squeeze(data(ids2comp(combi, 2), :, triali));
                    
                    if isempty(find(isnan(amyD))) & isempty(find(isnan(hpcD))) 
                       amyDF = eegfilt (amyD,1000, foi(1), foi(end));
                       hpcDF = eegfilt (hpcD,1000, foi(1), foi(end));
                       
                       %hilbAMY = hilbert(amyDF);
                       %hilbHPC = hilbert(hpcDF);
                       %angleAMY = angle(hilbAMY(toi));
                       %angleHPC = angle(hilbHPC(toi));
                       %diffPha = angleAMY - angleHPC;
                       %PLV_M1(triali, combi, chanj) = abs(mean(exp(1i*(diffPha)))); % just to test, works the same 
                       %PLI_M1(triali, combi, chanj) = abs(mean(sign(imag(exp(1i*diffPha))))); % just to test, works the same 

                       % % % %  baseline
                       %angleAMYBas = angle(hilbAMY(toiB));
                       %angleHPCBas = angle(hilbHPC(toiB));
                       %diffPhaBas = angleAMYBas - angleHPCBas;

                       hilbAMY = hilbert(amyDF);
                       hilbHPC = hilbert(hpcDF);
                       cdd = hilbAMY.* conj(hilbHPC);% cross-spectral density
                       cdd = cdd(toi);
                       cdi = imag(cdd);
                       PLV_M2(triali, combi,  :) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
                       PLI_M2(triali, combi,  :) = abs(mean(sign(imag(cdd)),2));
                       WPLI(triali, combi, :) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)

                       % % % baseline
                       cddBas = hilbAMY.* conj(hilbHPC);% cross-spectral density
                       cddBas = cddBas(toiB);
                       cdiBas = imag(cddBas);
                       PLV_M2Bas(triali, combi, :) = abs(mean(exp(1i*angle(cddBas)),2)); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
                       PLI_M2Bas(triali, combi, :) = abs(mean(sign(imag(cddBas)),2));
                       WPLIBas(triali,  combi, :) = abs( mean( abs(cdiBas).*sign(cdiBas) ,2) )./mean(abs(cdiBas),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
                       

                       % % compute coherence
                       spec1 = mean(hilbAMY(toi).*conj(hilbAMY(toi)),2);
                       spec2 = mean(hilbHPC(toi).*conj(hilbHPC(toi)),2);
                       specX = abs(mean(hilbAMY(toi).*conj(hilbHPC(toi)),2)).^2;
                       COH(triali, combi, :) = specX./ (spec1.*spec2);

                       % % compute coherence baseline 
                       spec1 = mean(hilbAMY(toiB).*conj(hilbAMY(toiB)),2);
                       spec2 = mean(hilbHPC(toiB).*conj(hilbHPC(toiB)),2);
                       specX = abs(mean(hilbAMY(toiB).*conj(hilbHPC(toiB)),2)).^2;
                       COHBas(triali, combi, :) = specX./ (spec1.*spec2);
                       

                    else 
                       PLV_M2(triali, combi, :) = nan;
                       PLI_M2(triali, combi, :) = nan;
                       WPLI(triali, combi, :) = nan;
                       COH(triali, combi, :) = nan;

                       PLV_M2Bas(triali, combi, :) = nan;
                       PLI_M2Bas(triali, combi, :) = nan;
                       WPLIBas(triali, combi, :) = nan;
                       COHBas(triali, combi, :) = nan;
                   end
               end
               
           end
           CON(1,:,:) = PLV_M2; 
           CON(2,:,:) = PLI_M2; 
           CON(3,:,:) = WPLI; 
           CON(4,:,:) = COH; 
           CON(5,:,:) = PLV_M2Bas; 
           CON(6,:,:) = PLI_M2Bas; 
           CON(7,:,:) = WPLIBas; 
           CON(8,:,:) = COHBas; 
           

end
                