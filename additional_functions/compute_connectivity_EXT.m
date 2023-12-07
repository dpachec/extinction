function [CON] = compute_connectivity_EXT(data, toi, toiB, foi, ids2comp, mode2u)
        %clear PLV_M1 PLV_M2 PLI_M1 PLI_M2 WPLI COH PLV_M1Bas PLV_M2Bas PLI_M1Bas PLI_M2Bas WPLIBas COHBas

if strcmp(mode2u, 'T')

        nTrials = size(data, 3); 
            for combi = 1:length(ids2comp(:, 1))
                
                disp(['Comb: ' num2str(combi) '/' num2str(length(ids2comp(:, 1))) ])
                 parfor triali = 1:nTrials
                    
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
                       PLV_M2(triali, combi,  :) = abs(mean(exp(1i*angle(cdd)))); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
                       PLI_M2(triali, combi,  :) = abs(mean(sign(imag(cdd))));
                       WPLI(triali, combi, :) = abs( mean( abs(cdi).*sign(cdi)))./mean(abs(cdi));% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)

                       % % % baseline
                       cddBas = hilbAMY.* conj(hilbHPC);% cross-spectral density
                       cddBas = cddBas(toiB);
                       cdiBas = imag(cddBas);
                       PLV_M2Bas(triali, combi, :) = abs(mean(exp(1i*angle(cddBas)))); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
                       PLI_M2Bas(triali, combi, :) = abs(mean(sign(imag(cddBas))));
                       WPLIBas(triali,  combi, :) = abs( mean( abs(cdiBas).*sign(cdiBas)) )./mean(abs(cdiBas));% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
                       

                       % % compute coherence
                       spec1 = mean(hilbAMY(toi).*conj(hilbAMY(toi)));
                       spec2 = mean(hilbHPC(toi).*conj(hilbHPC(toi)));
                       specX = abs(mean(hilbAMY(toi).*conj(hilbHPC(toi)))).^2;
                       COH(triali, combi, :) = specX./ (spec1.*spec2);

                       % % % DIFFERENT TESTS
% %                         spec1 = sum(hilbAMY(toi).*conj(hilbAMY(toi)),2);
% %                         spec2 = sum(hilbHPC(toi).*conj(hilbHPC(toi)),2);
% %                         specX = sum(hilbAMY(toi).*conj(hilbHPC(toi)),2);
% %                         COH(triali, combi, :)= abs(specX./sqrt(spec1.*spec2)).^2;
% % % 

                       % % compute coherence baseline 
                       spec1 = mean(hilbAMY(toiB).*conj(hilbAMY(toiB)));
                       spec2 = mean(hilbHPC(toiB).*conj(hilbHPC(toiB)));
                       specX = abs(mean(hilbAMY(toiB).*conj(hilbHPC(toiB)))).^2;
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
           


elseif strcmp(mode2u, 'TFR')


nTimepoints = 2500; %%all2all file is stored within this cell array
win_width   = 500; 
mf          = 100; 
bins        =  floor ( (nTimepoints/mf)- win_width/mf+1 );

nFreqs = 40; 
nTrials = size(data, 3); 
PLV_M2 = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
PLI_M2 = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
wPLI = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
COH = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
parfor combi = 1:length(ids2comp(:, 1))
    disp(['Comb: ' num2str(combi) '/' num2str(length(ids2comp(:, 1))) ])
    
    PLV_M2_p3 = zeros(nTrials,nFreqs, bins); 
    PLI_M2_p3 = zeros(nTrials,nFreqs, bins); 
    WPLI_M2_p3 = zeros(nTrials, nFreqs,bins); 
    COH_M2_p3 = zeros(nTrials, nFreqs,bins); 
    for triali = 1:nTrials
        amyDP = squeeze(data(ids2comp(combi, 1), :, triali));
        hpcDP = squeeze(data(ids2comp(combi, 2), :, triali));
        
        if isempty(find(isnan(amyDP))) & isempty(find(isnan(hpcDP))) 
            PLV_M2_p2 = zeros(nFreqs, bins); 
            PLI_M2_p2 = zeros(nFreqs, bins); 
            WPLI_M2_p2 = zeros(nFreqs,bins); 
            COH_M2_p2 = zeros(nFreqs,bins); 
            for freqi = 3:nFreqs
                amyDF = eegfilt (amyDP,1000, freqi-1, freqi+1);
                hpcDF = eegfilt (hpcDP,1000, freqi-1, freqi+1);
                hilbAMY = hilbert(amyDF);
                hilbHPC = hilbert(hpcDF);
                hilbAMY = hilbAMY(2501:5000);
                hilbHPC = hilbHPC(2501:5000);
                cddP = hilbAMY.* conj(hilbHPC);% cross-spectral density
     

                PLV_M2_p1 = zeros( 1, bins); 
                PLI_M2_p1 = zeros( 1, bins); 
                WPLI_M2_p1 = zeros( 1, bins); 
                COH_M2_p1 = zeros( 1, bins); 
                 for timei = 1:bins 
                     timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1; 
                     cdd = cddP(timeBinsi);
                     cdi = imag(cdd);
                     PLV_M2_p1( timei) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
                     PLI_M2_p1( timei) = abs(mean(sign(imag(cdd)),2));
                     WPLI_M2_p1(  timei) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
                   
                    % % compute coherence
                    spec1 = mean(hilbAMY(timeBinsi).*conj(hilbAMY(timeBinsi)),2);
                    spec2 = mean(hilbHPC(timeBinsi).*conj(hilbHPC(timeBinsi)),2);
                    specX = abs(mean(hilbAMY(timeBinsi).*conj(hilbHPC(timeBinsi)),2)).^2;
                    COH_M2_p1(  timei) = specX./ (spec1.*spec2);
                 end
                PLV_M2_p2(freqi, :) = PLV_M2_p1; 
                PLI_M2_p2(freqi, :) = PLI_M2_p1; 
                WPLI_M2_p2(freqi, :) = WPLI_M2_p1; 
                COH_M2_p2(freqi, :) = COH_M2_p1; 
            end
        else 
            PLV_M2_p2 = nan(nFreqs, bins); 
            PLI_M2_p2 = nan(nFreqs, bins); 
            WPLI_M2_p2 = nan(nFreqs,bins); 
            COH_M2_p2 = nan(nFreqs,bins); 
        end

    PLV_M2_p3(triali,:,:,:) = PLV_M2_p2; 
    PLI_M2_p3(triali,:,:,:) = PLI_M2_p2; 
    WPLI_M2_p3(triali,:,:,:) = WPLI_M2_p2; 
    COH_M2_p3(triali,:,:,:) = COH_M2_p2; 
    end

    CON1(:, :, combi, :) = PLV_M2_p3;
    CON2(:, :, combi, :) = PLI_M2_p3;
    CON3(:, :, combi, :) = WPLI_M2_p3;
    CON4(:, :, combi, :) = COH_M2_p3;
   
end

CON(1,:,:,:,:) = CON1; 
CON(2,:,:,:,:) = CON2; 
CON(3,:,:,:,:) = CON3; 
CON(4,:,:,:,:) = CON4; 



% % % % % 
% % % % % nTimepoints = 2500; %%all2all file is stored within this cell array
% % % % % win_width   = 500; 
% % % % % mf          = 100; 
% % % % % bins        =  floor ( (nTimepoints/mf)- win_width/mf+1 );
% % % % % 
% % % % % nFreqs = 30; 
% % % % % nTrials = size(data, 3); 
% % % % % PLV_M2 = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
% % % % % PLI_M2 = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
% % % % % wPLI = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
% % % % % COH = zeros(nTrials, nFreqs, length(ids2comp(:, 1)), bins);
% % % % % for combi = 1:length(ids2comp(:, 1))
% % % % %     disp(['Comb: ' num2str(combi) '/' num2str(length(ids2comp(:, 1))) ])
% % % % %     
% % % % %     for triali = 1:nTrials
% % % % %         amyDP = squeeze(data(ids2comp(combi, 1), :, triali));
% % % % %         hpcDP = squeeze(data(ids2comp(combi, 2), :, triali));
% % % % %         
% % % % %         if isempty(find(isnan(amyDP))) & isempty(find(isnan(hpcDP))) 
% % % % %             for freqi = 3:nFreqs
% % % % %                 amyDF = eegfilt (amyDP,1000, freqi-1, freqi+1);
% % % % %                 hpcDF = eegfilt (hpcDP,1000, freqi-1, freqi+1);
% % % % %                 hilbAMY = hilbert(amyDF);
% % % % %                 hilbHPC = hilbert(hpcDF);
% % % % %                 hilbAMY = hilbAMY(2501:5000);
% % % % %                 hilbHPC = hilbHPC(2501:5000);
% % % % %                 cddP = hilbAMY.* conj(hilbHPC);% cross-spectral density
% % % % %      
% % % % %                  for timei = 1:bins 
% % % % %                      timeBinsi = (timei*mf) - (mf-1):(timei*mf - (mf-1) )+win_width-1; 
% % % % %                      cdd = cddP(timeBinsi);
% % % % %                      cdi = imag(cdd);
% % % % %                      PLV_M2(triali, freqi, combi,  timei) = abs(mean(exp(1i*angle(cdd)),2)); % note: equivalent to PLV_SI(combi, chanj, :) = abs(mean(exp(1i*(diffPha)), 2)); %PLV
% % % % %                      PLI_M2(triali, freqi,  combi,  timei) = abs(mean(sign(imag(cdd)),2));
% % % % %                      WPLI(triali,  freqi, combi,  timei) = abs( mean( abs(cdi).*sign(cdi) ,2) )./mean(abs(cdi),2);% weighted phase-lag index (eq. 8 in Vink et al. NeuroImage 2011)
% % % % %                    
% % % % %                     % % compute coherence
% % % % %                     spec1 = mean(hilbAMY(timeBinsi).*conj(hilbAMY(timeBinsi)),2);
% % % % %                     spec2 = mean(hilbHPC(timeBinsi).*conj(hilbHPC(timeBinsi)),2);
% % % % %                     specX = abs(mean(hilbAMY(timeBinsi).*conj(hilbHPC(timeBinsi)),2)).^2;
% % % % %                     COH(triali,  freqi, combi, timei) = specX./ (spec1.*spec2);
% % % % %                  end
% % % % %             end
% % % % %         else 
% % % % %             PLV_M2(triali,  :, combi, :) = nan(nFreqs, bins);
% % % % %             PLI_M2(triali,  :, combi, :) = nan(nFreqs, bins);
% % % % %             WPLI(triali,  :, combi, :) = nan(nFreqs, bins);
% % % % %             COH(triali,  :, combi, :) = nan(nFreqs, bins);
% % % % %         end
% % % % %    end
% % % % %    
% % % % % end
% % % % % 
% % % % % 
% % % % % CON(1,:,:,:,:) = PLV_M2; 
% % % % % CON(2,:,:,:,:) = PLI_M2; 
% % % % % CON(3,:,:,:,:) = WPLI; 
% % % % % CON(4,:,:,:,:) = COH;            

end


end
                