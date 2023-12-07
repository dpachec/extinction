%% Simulate data - Simple system

clear
cfg = [];
cfg.method = 'ar';
cfg.ntrials = 200;
cfg.triallength = 1;
cfg.fsample = 200;
cfg.nsignal = 2;
%AR model at time lag 1 (matrix of size Signals, Signals, Model order)
cfg.params(:,:,1) = [0.55 0.2; %Singla 1 influences itself at t-1 with 0.55
                     0  0.55]; %DIagonal > self term. Off diagonal > Signal 2 influences signal 1 at t-1 with a beta of 0.2
cfg.params(:,:,2) = [-0.8 -0.1;%Singla 1 influences itself at t-1 with -.8
                     0 -0.8]; %off-diagonal entries simulate 3->1 and 3->2 influence at the second time delay
cfg.noisecov = [1 0;
                0 1];
data = ft_connectivitysimulation(cfg);

%% Simulate data - Node 3 gives common input to the other nodes (nodes 1 and 2) at a time delay of 1 and 2 samples

clear
cfg = [];
cfg.method = 'ar';
cfg.ntrials = 200;
cfg.triallength = 1;
cfg.fsample = 200;
cfg.nsignal = 3;
cfg.params(:,:,1) = [0.55 0 0.25;
                     0 0.55 0.25;
                     0 0 0.55]; %off-diagonal entries simulate 3->1 and 3->2 influence at the first time delay
cfg.params(:,:,2) = [-0.8 0 -0.1;
                     0 -0.8 -0.1;
                     0 0 -0.8]; %off-diagonal entries simulate 3->1 and 3->2 influence at the second time delay
cfg.noisecov = [1 0 0;
                0 1 0;
                0 0 1];
data = ft_connectivitysimulation(cfg);

%%
for i=1:size(data.trial,2)
  EEG.data(:,:,i) = single(data.trial{i});
end

EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculate power, coherence, and Granger causality based on parametric and
%non-parametric estimates as in Figure 9 b and c
%calculate the fourier coefficients (non-parametric derivation of power)
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'dpss';
cfg.output = 'fourier';
cfg.tapsmofrq = 3;
cfg.foilim = [0 100];
freq = ft_freqanalysis(cfg, data);
%freqdescriptives calculates the power spectrum
cfg = [];
fd = ft_freqdescriptives(cfg, freq);

%% Parametric (auto-regressive model based) derivation of AR coefficients
%multivariate analysis will compute the auto-regressive coefficients and associated noise covariance matrix

cfg = [];
cfg.order = 2; %model order of 2, this is known a priori (we simulated the data using a model order of 2)
cfg.method = 'bsmart'
mdata = ft_mvaranalysis(cfg, data);



%% calculate cross-spectral density and transfer functions associated with the auto-regressive model
cfg = [];
cfg.method = 'mvar';
cfg.foi = [0:100];
mfreq = ft_freqanalysis(cfg, mdata);


%% Phase-slope index calculation
cfg = [];
cfg.method = 'psi';
cfg.bandwidth = 4;
psi1 = ft_connectivityanalysis(cfg, freq);
%% Coherence calculation
cfg = [];
cfg.method = 'coh';
cfg.complex = 'abs';
coh1 = ft_connectivityanalysis(cfg, freq);
coh2 = ft_connectivityanalysis(cfg, mfreq);

%% Imaginary part of coherency
cfg = [];
cfg.method = 'coh';
cfg.complex = 'imag';
icoh1 = ft_connectivityanalysis(cfg, freq);

%% Partial coherence calculation
cfg = [];
cfg.method = 'coh';
cfg.partchannel = 'signal003';
cfg.complex = 'abs';
pcoh1 = ft_connectivityanalysis(cfg, freq);

%% Granger causality calculation
cfg = [];
cfg.method = 'granger';
cfg.granger.sfmethod = 'multivariate';
g1 = ft_connectivityanalysis(cfg, freq);
g1 = ft_checkdata(g1, 'cmbrepresentation', 'full');
g2 = ft_connectivityanalysis(cfg, mfreq);


%% PLOT ONLY 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now plot the output for the various connectivity measures as in Figure 9 b and c
%output variables 1 - nonparam 2 - param
figure;
plot(fd.freq, fd.powspctrm(1,:)); hold on
plot(fd.freq, fd.powspctrm(2,:),'r');
legend('Power ch 1', 'Power ch 2');
title('Nonparametric Power');

figure;
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(1,1,:)))); hold on;
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(2,2,:))),'r');
legend('Power ch 1', 'Power ch 2');
title('Parametric Power');

figure;
plot(g1.freq,squeeze(coh1.cohspctrm(1,2,:))); hold on;
plot(g1.freq,squeeze(abs(icoh1.cohspctrm(1,2,:))),'g');
legend('1-2', '1-2 imaginary');
title('Nonparametric Coherence spectrum');

figure;
plot(g1.freq,squeeze(coh2.cohspctrm(1,2,:))); hold on
legend('1-2');
title('Parametric Coherence spectrum');

figure;
plot(g1.freq,squeeze(psi1.psispctrm(1,2,:))); hold on;
legend('1->2');title('PSI nonparametric');

figure;
plot(g1.freq,squeeze(g1.grangerspctrm(1,2,:)));hold on
plot(g1.freq,squeeze(g1.grangerspctrm(2,1,:)),'r');
title('Granger nonparametric estimates');
legend('1->2','2->1')

figure;
plot(g1.freq,squeeze(g2.grangerspctrm(1,2,:)));hold on
plot(g1.freq,squeeze(g2.grangerspctrm(2,1,:)),'r');
legend('1->2','2->1');
title('Granger parametric estimates ');




%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%now plot the output for the various connectivity measures as in Figure 9 b and c
%output variables 1 - nonparam 2 - param
figure;plot(fd.freq, fd.powspctrm(1,:)); hold on
plot(fd.freq, fd.powspctrm(2,:),'r');
plot(fd.freq, fd.powspctrm(3,:),'k');
legend('Power ch 1', 'Power ch 2','Power ch 3');
title('Nonparametric Power');
figure;plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(1,1,:)))); hold on;
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(2,2,:))),'r');
plot(mfreq.freq, squeeze(abs(mfreq.crsspctrm(3,3,:))),'k');
legend('Power ch 1', 'Power ch 2','Power ch 3');
title('Parametric Power');
figure;plot(g1.freq,squeeze(coh1.cohspctrm(1,2,:))); hold on;
plot(g1.freq,squeeze(coh1.cohspctrm(1,3,:)),'r');
plot(g1.freq,squeeze(coh1.cohspctrm(2,3,:)),'k');
plot(g1.freq,squeeze(abs(icoh1.cohspctrm(1,2,:))),'g');
plot(g1.freq,squeeze(pcoh1.cohspctrm(1,2,:)),'m');
legend('1-2','1-3','2-3', '1-2 imaginary', '1-2 | 3');
title('Nonparametric Coherence spectrum');
figure;plot(g1.freq,squeeze(coh2.cohspctrm(1,2,:))); hold on
plot(g1.freq,squeeze(coh2.cohspctrm(1,3,:)),'r');
plot(g1.freq,squeeze(coh2.cohspctrm(2,3,:)),'k');
legend('1-2','1-3','2-3');
title('Parametric Coherence spectrum');
figure;plot(g1.freq,squeeze(psi1.psispctrm(1,2,:))); hold on;
plot(g1.freq,squeeze(psi1.psispctrm(1,3,:)),'r');
plot(g1.freq,squeeze(psi1.psispctrm(2,3,:)),'k');
legend('1->2','1->3','3->2');title('PSI nonparametric');
figure;plot(g1.freq,squeeze(g1.grangerspctrm(1,2,:)));hold on
plot(g1.freq,squeeze(g1.grangerspctrm(2,1,:)),'r');
plot(g1.freq,squeeze(g1.grangerspctrm(3,1,:)),'k');







































%%
% Simulated data with directed connections
% We will first simulate some data with a known connectivity structure built in. This way we know what to expect in terms of connectivity. 
% To simulate data we use ft_connectivitysimulation. We will use an order 2 multivariate autoregressive model. The necessary ingredients are 
% a set of NxN coefficient matrices, one matrix for each time lag. These coefficients need to be stored in the cfg.param field. 
% Next to the coefficients we have to specify the NxN covariance matrix of the innovation noise. This matrix needs to be stored in the 
% cfg.noisecov field. The model we are going to use to simulate the data is as follows:
% 
% x(t) = 0.8*x(t-1) - 0.5*x(t-2)
% y(t) = 0.9*y(t-1) + 0.5*z(t-1) - 0.8*y(t-2)
% z(t) = 0.5*z(t-1) + 0.4*x(t-1) - 0.2*z(t-2)

clear 
cfg             = [];
cfg.ntrials     = 500;
cfg.triallength = 1;
cfg.fsample     = 200;
cfg.nsignal     = 3;
cfg.method      = 'ar';

cfg.params(:,:,1) = [ 0.8    0    0 ; 
                        0  0.9  0.5 ;
                      0.4    0  0.5];
                      
cfg.params(:,:,2) = [-0.5    0    0 ; 
                        0 -0.8    0 ; 
                        0    0 -0.2];
                        
cfg.noisecov      = [ 0.3    0    0 ;
                        0    1    0 ;
                        0    0  0.2];

data              = ft_connectivitysimulation(cfg);

%% The simulated data consists of 3 channels in 500 trials. You can easily visualize the data for example in the first trial using
% 
figure
plot(data.time{1}, data.trial{1}) 
legend(data.label)
xlabel('time (s)')

%%
% or browse through the complete data using
% 
cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly' 
ft_databrowser(cfg, data);

%% 
% 
% Computation of the multivariate autoregressive model
% To be able to compute spectrally resolved Granger causality, or other frequency-domain directional measures of connectivity, 
% we have to fit an autoregressive model to the data. This is done using the ft_mvaranalysis function.
% % For the actual computation of the autoregressive coefficients FieldTrip makes use of an implementation from third party toolboxes. 
% At present ft_mvaranalysis supports the biosig and bsmart toolboxes for these computations.
% % In this tutorial we will use the bsmart toolbox. The relevant functions have been included in the FieldTrip release in the
% fieldtrip/external/bsmart directory.
% 
cfg         = [];
cfg.order   = 5;
cfg.toolbox = 'bsmart';
mdata       = ft_mvaranalysis(cfg, data);

%%
% The resulting variable mdata contains a description of the data in terms of a multivariate autoregressive model. For each time-lag up
% to the model order (which is 5 in this case), a 3*3 matrix of coefficients is outputted. 
% The noisecov-field contains covariance matrix of the model's residuals.
% 
% Exercise 1
% Compare the parameters specified for the simulation with the estimated coefficients and discuss.
% 



%%  Computation of the spectral transfer function
% From the autoregressive coefficients it is now possible to compute the spectral transfer matrix, for which we use ft_freqanalysis.
% 
cfg        = [];
cfg.method = 'mvar';
mfreq      = ft_freqanalysis(cfg, mdata);
 
% The resulting mfreq data structure contains the pairwise transfer function between the 3 channels for 101 frequencies.


%%  It is also possible to compute the spectral transfer function using non-parametric spectral factorization of the cross-spectral density matrix. 
% For this, we need a Fourier decomposition of the data.  
% Non-parametric computation of the cross-spectral density matrix
% Some connectivity metrics can be computed from a non-parametric spectral estimate 
% (i.e. after the application of the FFT-algorithm and conjugate multiplication to get cross-spectral densities), 
% such as coherence, phase-locking value and phase slope index. The following part computes the fourier-representation of the data using 
% ft_freqanalysis.
% It is not necessary to compute the cross-spectral density at this stage, because the function used in the next step, 
% ft_connectivityanalysis, contains functionality to compute the cross-spectral density from the fourier coefficients.
% 
cfg           = [];
cfg.method    = 'mtmfft';
cfg.taper     = 'dpss';
cfg.output    = 'fourier';
cfg.tapsmofrq = 2;
freq          = ft_freqanalysis(cfg, data);

% The resulting freq structure contains the spectral estimate for 3 tapers in each of the 500 trials (hence 1500 estimates), for each of the 3 channels 
% and for 101 frequencies.
% 
%% Computation and inspection of the connectivity measures
% The actual computation of the connectivity metric is done by ft_connectivityanalysis. This function is transparent to the type of input data, 
% i.e. provided the input data allows the requested metric to be computed, the metric will be calculated. 
% Here, we provide an example for the computation and visualization of the coherence coefficient.
% 
cfg           = [];
cfg.method    = 'coh';
coh           = ft_connectivityanalysis(cfg, freq);
cohm          = ft_connectivityanalysis(cfg, mfreq);
% Subsequently, the data can be visualized using ft_connectivityplot.
% 
cfg           = [];
cfg.parameter = 'cohspctrm';
cfg.zlim      = [0 1];
figure()
ft_connectivityplot(cfg, coh, cohm);

% 
% The coherence measure is a symmetric measure, which means that it does not provide information regarding the direction of information flow 
% between any pair of signals. In order to analyze directionality in interactions, measures based on the concept of granger causality can be computed. 
% These measures are based on an estimate of the spectral transfer matrix, which can be computed in a straightforward way from the multivariate 
% autoregressive model fitted to the data.
% 
cfg           = [];
cfg.method    = 'granger';
granger       = ft_connectivityanalysis(cfg, freq);
grangerm       = ft_connectivityanalysis(cfg, mfreq);

cfg           = [];
cfg.parameter = 'grangerspctrm';
cfg.zlim      = [0 1];
figure()
ft_connectivityplot(cfg, granger, grangerm);


% Instead of plotting it with ft_connectivityplot, you can use the following low-level MATLAB plotting code which gives a better understanding 
% of the numerical representation of the results.
figure
for row=1:3
for col=1:3
  subplot(3,3,(row-1)*3+col);
  plot(granger.freq, squeeze(granger.grangerspctrm(row,col,:)))
  ylim([0 1])
end
end


cfgp = [];
cfgp.method = 'psi';
cfgp.bandwidth = 5;
psi = ft_connectivityanalysis(cfgp, freq);

cfg           = [];
cfg.parameter = 'psispctrm';
%cfg.zlim      = [0 1];
figure()
ft_connectivityplot(cfg, psi);


cfg           = [];
cfg.method    = 'plv';
plv           = ft_connectivityanalysis(cfg, freq);

cfg.parameter = 'plvspctrm';
cfg.zlim      = [0 1];
figure()
ft_connectivityplot(cfg, plv);

%% create EEGLAB struct with the data

for i=1:size(data.trial,2)
  EEG.data(:,:,i) = single(data.trial{i});
end

EEG.comments   = 'preprocessed with fieldtrip';
EEG.nbchan     = size(data.trial{1},1);
EEG.trials     = size(data.trial,2);
EEG.pnts       = size(data.trial{1},2);
EEG.srate      = data.fsample;
EEG.xmin       = data.time{1}(1);
EEG.xmax       = data.time{1}(end);
EEG.times      = data.time{1};


%% compute psi manually
clc, clear PSI
freqb = [1:1:100]'; freqb(:, 2) = 2:1:101'; freqb(:, 3) = 3:1:102'; freqb(:, 4) = 4:1:103'; freqb(:, 5) = 5:1:104'; 
%psiTmp = data2psiX(EEG.data, 200, freqb, 0);
data2 = EEG.data([1 3],:,:); data2 = reshape(data2, 2, [] )';
[psiTmp, stdpsi, psisum, stdpsisum]=data2psi(data2,251,300,freqb);
%and psi has then 3 indices with the last one refering to the row in freqs: 
psiTmp = psiTmp./(stdpsi+eps);

PSI = squeeze(psiTmp(2,1, :));


%%
figure
plot(PSI); hold on; 
plot(psi.freq(1:85), squeeze(psi.psispctrm(1,3,1:85)))
legend({'manual' 'fieldtrip'})

%% % trivial example for flow from channel 1 to channel 2. 
n=10000;x=randn(n+1,1);data=[x(2:n+1),x(1:n)]; 

% parameters for PSI-calculation 
segleng=100;epleng=200;

% calculation of PSI. The last argument is empty - meaning that 
% PSI is calculated over all frequencies
[psi, stdpsi, psisum, stdpsisum]=data2psi(data,segleng,epleng,[]);

% note, psi, as calculated by data2psi corresponds to \hat{\PSI} 
% in the paper, i.e., it is not normalized. The final is 
% the normalized version given by: 
psi./(stdpsi+eps)

%% To calculate psi in a band set, e.g., 
freqs=[5:10];
[psi, stdpsi, psisum, stdpsisum]=data2psi(data,segleng,epleng,freqs);
%with result:
psi./(stdpsi+eps)
% In this example, the flow is estimated to go from channel 
% 1 to channel 2 because the matrix element psi(1,2) is positive. 

%% You can also calculate many bands at once, e.g. 
freqs=[[5:10];[6:11];[7:12]];
[psi, stdpsi, psisum, stdpsisum]=data2psi(data,segleng,epleng,freqs);
%and psi has then 3 indices with the last one refering to the row in freqs: 
psi./(stdpsi+eps)

%%
figure
plot(PLV); hold on; 
plot(PLV2); hold on; 
plot(plv.freq(1:85), squeeze(plv.plvspctrm(1,3,1:85)))
legend({'manual' 'fieldtrip'})
%% COMPUTE COHERENCE MANUALLY 


%%  compute for all trials together
clear COH PLV
clc
for freqi = 6:80
    
    s1 = squeeze(EEG.data(1,:,:));
    s2 = squeeze(EEG.data(3,:,:));
    s1B = reshape (s1, 1, []);
    s2B = reshape (s2, 1, []);
    signal1 = eegfilt (s1B,200, freqi-5, freqi+5);
    signal2 = eegfilt (s2B,200, freqi-5, freqi+5);
    hilb1 = hilbert(signal1);
    hilb2 = hilbert(signal2);
    spec1 = mean(hilb1.*conj(hilb1));
    spec2 = mean(hilb2.*conj(hilb2));
    specX = abs(mean(hilb1.*conj(hilb2))).^2;
    COH(freqi)= specX./ (spec1.*spec2);
    
    cdd = hilb1.* conj(hilb2);
    PLV(freqi) = abs(mean(exp(1i*angle(cdd)))); 

    
    
end

disp('done')


%%  compute for all trials together 2
clear COH2 PLV2
clc

% wavelet and FFT parameters
time          = 0:1/EEG.srate:1; %time = EEG.times; 
half_wavelet  = (length(time)-1)/2;
n_wavelet     = length(time);
n_data        = EEG.pnts*EEG.trials;
n_convolution = n_wavelet+n_data-1;
%num_cycles    = logspace(log10(3),log10(10),85);
num_cycles    = linspace(3,10,85);

% data FFTs
data_fft1 = fft(reshape(EEG.data(1,:,:),1,n_data),n_convolution);
data_fft2 = fft(reshape(EEG.data(3,:,:),1,n_data),n_convolution);


% initialize
spectcoher = zeros(length(freqi),EEG.pnts);

for freqi = 1:85
    % create wavelet and take FFT
    s = num_cycles(freqi)/(2*pi*freqi);
    wavelet_fft = fft( exp(2*1i*pi*freqi.*time) .* exp(-time.^2./(2*(s^2))),n_convolution);
    % phase angles from channel 1 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft1,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    %sig1 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    sig1 = convolution_result_fft; 
    % phase angles from channel 2 via convolution
    convolution_result_fft = ifft(wavelet_fft.*data_fft2,n_convolution);
    convolution_result_fft = convolution_result_fft(half_wavelet+1:end-half_wavelet);
    %sig2 = reshape(convolution_result_fft,EEG.pnts,EEG.trials);
    sig2 = convolution_result_fft; 

    spec1 = mean(sig1.*conj(sig1));
    spec2 = mean(sig2.*conj(sig2));
    specX = abs(mean(sig1.*conj(sig2))).^2;

    COH2(freqi)= specX./ (spec1.*spec2);

    cdd = sig1.* conj(sig2);
    PLV2(freqi) = abs(mean(exp(1i*angle(cdd)))); 

    
end

disp('done')
%%
figure
plot(COH); hold on; 
plot(COH2); 
plot(coh.freq(1:85), squeeze(coh.cohspctrm(1,3,1:85)))
%plot(coh.freq(1:85), squeeze(coh.cohspctrm(1,3,1:85).^2))
legend({'manual' 'fieldtrip'})

%%
figure
plot(PLV); hold on; 
plot(PLV2)
plot(plv.freq(1:85), squeeze(plv.plvspctrm(1,3,1:85)))
legend({'manual' 'fieldtrip'})

%% Exercise 2
% Discuss the differences between the granger causality spectra, and the coherence spectra.


%% Exercise 3
% Compute the following connectivity measures from the mfreq data, and visualize and discuss the results: partial directed coherence (pdc), directed transfer function (dtf), phase slope index (psi)
% 
% Simulated data with common pick-up and different noise levels
%  this is under progress
% 
% When working with electrophysiological data (EEG/MEG/LFP) the signals that are picked up by the individual channels invariably consist of instantaneous mixtures of the underlying source signals. This mixing can severely affect the outcome of connectivity analysis, and thus affects the interpretation. We will demonstrate this by simulating data in 2 channels, where each of the channels consists of a weighted combination of temporally white noise unique to each of the channels, and a common input of a band-limited signal (filtered between 15 and 25 Hz). We will compute connectivity between these channels, and show that the common input can give rise to spurious estimates of connectivity.
% 
% % create some instantaneously mixed data
% 
% define some variables locally
nTrials  = 100;
nSamples = 1000;
fsample  = 1000;

% mixing matrix
mixing   = [0.8 0.2 0;
              0 0.2 0.8];

data       = [];
data.trial = cell(1,nTrials);
data.time  = cell(1,nTrials);
for k = 1:nTrials
  dat = randn(3, nSamples);
  dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [15 25]);
  dat = 0.2.*(dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
  data.trial{k} = mixing * dat;
  data.time{k}  = (0:nSamples-1)./fsample;
end
data.label = {'chan1' 'chan2'}';

figure;plot(dat'+repmat([0 1 2],[nSamples 1]));
title('original ''sources''');

figure;plot((mixing*dat)'+repmat([0 1],[nSamples 1])); 
axis([0 1000 -1 2]);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('mixed ''sources''');
 

% do spectral analysis
cfg = [];
cfg.method    = 'mtmfft';
cfg.output    = 'fourier';
cfg.foilim    = [0 200];
cfg.tapsmofrq = 5;
freq          = ft_freqanalysis(cfg, data);
fd            = ft_freqdescriptives(cfg, freq);

figure;plot(fd.freq, fd.powspctrm);
set(findobj(gcf,'color',[0 0.5 0]), 'color', [1 0 0]);
title('powerpectrum');

% 
% compute connectivity
cfg = [];
cfg.method = 'granger';
g = ft_connectivityanalysis(cfg, freq);
cfg.method = 'coh';
c = ft_connectivityanalysis(cfg, freq);
% visualize the results
cfg = [];
cfg.parameter = 'grangerspctrm';
ft_connectivityplot(cfg, g);
cfg.parameter = 'cohspctrm';
ft_connectivityplot(cfg, c);


% Exercise 4
% Simulate new data using the following mixing matrix:
% [0.9 0.1 0;0 0.2 0.8] 
% and recompute the connectivity measures. Discuss what you see.
% 
% Exercise 5
% Play a bit with the parameters in the mixing matrix and see what is the effect on the estimated connectivity.
% Exercise 6
% Simulate new data where the 2 mixed signals are created from 4 underlying sources, and where two of these sources are common input to both signals, and where these two sources are temporally shifted copies of one another.
% Hint: the mixing matrix could look like this:
% 
% [a b c 0; 0 d e f];
% and the trials could be created like this:
% 
% for k = 1:nTrials
%   dat = randn(4, nSamples+10);
%   dat(2,:) = ft_preproc_bandpassfilter(dat(2,:), 1000, [15 25]);
%   dat(3,1:(nSamples)) = dat(2,11:(nSamples+10)); 
%   dat = dat(:,1:1000);
%   dat = 0.2.*(dat-repmat(mean(dat,2),[1 nSamples]))./repmat(std(dat,[],2),[1 nSamples]);
%   data.trial{k} = mixing * dat;
%   data.time{k}  = (0:nSamples-1)./fsample;
% end
% Compute connectivity between the signals and discuss what you observe. In particular, also compute measures of directed interaction.
% 
% 
% Connectivity between MEG virtual channel and EMG
% The previous two examples were using simulated data, either with a clear directed connectivity structure, or with a trivial pick-up of a common source in two channels. We will now continue with connectivity analysis on real MEG data. The dataset is the same as the one used in the Analysis of corticomuscular coherence tutorial.
% 
% In short, the dataset consists of combined MEG and EMG recordings while the subject lifted his right hand. The coherence tutorial introduction contains a more elaborate description of the experiment and the dataset and a detailed analysis can be found in the corresponding paper 1). Due to the long distance between the EMG and the MEG, there is no volume conduction and hence no common pick-up. Hence this dataset lends itself well for connectivity analysis. But rather than doing an analysis between the EMG and one of the MEG channels (as in the original study), we will extract the cortical activity using a beamformer virtual channel.
% 
% 
% Compute the spatial filter for the region of interest
% We start with determining the motor cortex as region of interest. At the end of the coherence tutorial it is demonstrated how to make a 3-D reconstruction of the cortico-muscolar coherence (CMC) using the DICS algorithm. That source reconstruction serves as starting point for this analysis.
% 
% You can download the result from the DICS reconstruction from the FieldTrip ftp server (source.mat)
% 
% We will first determine the position on which the cortico-muscular coherence is the largest.
% 

cd(dccnpath('/home/common/matlab/fieldtrip/data/ftp/tutorial/connectivity'));
load source

[maxval, maxindx] = max(source.avg.coh);
maxpos = source.pos(maxindx,:);

% maxpos = 
%     4 -3 12
% The cortical position is expressed in individual subject head-coordinates and in centimeter. Relative to the center of the head (in between the ears) the position is 4 cm towards the nose, -3 towards the left side (i.e., 3 cm towards the right!) and 12 cm towards the vertex.
% 
% The ft_sourceanalysis methods are usually applied to the whole brain using a regular 3-D grid or using a triangulated cortical sheet. You can also just specify the location of a single or multiple points of interest with cfg.sourcemodel.pos and the LCMV beamformer will simply be performed at the location of interest.
% 
% The LCMV beamformer spatial filter for the location of interest will pass the activity at that location with unit-gain, while optimally suppressing all other noise and other source contributions to the MEG data. The LCMV implementation in FieldTrip requires the data covariance matrix to be computed with ft_timelockanalysis.
% 
% Rather than doing all the preprocessing again, you can download the preprocessed data from the FieldTrip ftp server (data.mat)
% 
load data
% 
%% compute the beamformer filter
cfg                   = [];
cfg.covariance        = 'yes';
cfg.channel           = 'MEG';
cfg.vartrllength      = 2;
cfg.covariancewindow  = 'all';
timelock              = ft_timelockanalysis(cfg, data);

cfg             = [];
cfg.method      = 'lcmv';
cfg.hdmfile     = 'SubjectCMC.hdm';
cfg.sourcemodel.pos = maxpos;
cfg.keepfilter  = 'yes';
source          = ft_sourceanalysis(cfg, timelock);
% The source reconstruction contains the estimated power and the source-level time-series of the averaged ERF, but here we are not interested in those. The cfg.keepfilter option results in the spatial filter being kept in the output source structure. That spatial can be used to reconstruct the single-trial time series as a virtual channel by multiplying it with the original MEG data.
% 
% 
% Extract the virtual channel time-series
%% construct the 3-D virtual channel at the location of interest
beamformer = source.avg.filter{1};

chansel = ft_channelselection('MEG', data.label); % find the names
chansel = match_str(data.label, chansel);         % find the indices

sourcedata = [];
sourcedata.label = {'x', 'y', 'z'};
sourcedata.time = data.time;
for i=1:length(data.trial)
  sourcedata.trial{i} = beamformer * data.trial{i}(chansel,:);
end
% The LCMV spatial filter is computed here without applying any time-domain filters. Consequently, it will have to suppress all noise in the data in all frequency bands. The spatial filter derived from the broadband data allows us to compute a broadband source level time-series.
% If you would know that the subsequent analysis would be limited to a specific frequency range in the data (e.g. everything above 30 Hz), you could first apply a filter using ft_preprocessing (e.g. cfg.hpfilter=yes and cfg.hpfreq=30) prior to computing the covariance and the spatial filter.
% 
% The sourcedata structure resembles the raw-data output of ft_preprocessing and consequently can be used in any follow-up function. You can for example visualize the single-trial virtual channel time-series using ft_databrowser:
% 
cfg = [];
cfg.viewmode = 'vertical';  % you can also specify 'butterfly'
ft_databrowser(cfg, sourcedata);

% 
% Notice that the reconstruction contains three channels, for the x-, the y- and the z-component of the equivalent current dipole source at the location of interest.
% 
% 
% Project along the strongest dipole direction
% The interpretation of connectivity is facilitated if we can compute it between two plain channels rather than between one channel versus a triplet of channels. Therefore we will project the time-series along the dipole direction that explains most variance. This projection is equivalent to determining the largest (temporal) eigenvector and can be computationally performed using the singular value decomposition (svd).
% 
%% construct a single virtual channel in the maximum power orientation
timeseries = cat(2, sourcedata.trial{:});

[u, s, v] = svd(timeseries, 'econ');

% whos u s v
%   Name           Size              Bytes  Class     Attributes
% 
%   s              3x3                  72  double              
%   u              3x3                  72  double              
%   v         196800x3             4723200  double            
% Matrix u contains the spatial decomposition, matrix v the temporal and on the diagonal of matrix s you can find the eigenvalues. See "help svd" for more details.
% 
% We now recompute the virtual channel time-series, but now only for the dipole direction that has the most power.
% 
% this is equal to the first column of matrix V, apart from the scaling with s(1,1)
timeseriesmaxproj = u(:,1)' * timeseries;
virtualchanneldata = [];
virtualchanneldata.label = {'cortex'};
virtualchanneldata.time = data.time;
for i=1:length(data.trial)
  virtualchanneldata.trial{i} = u(:,1)' * beamformer * data.trial{i}(chansel,:);
end
% Rather than using a sourcemodel in the beamformer that consists of all three (x, y, z) directions, you can also have the beamformer compute the filter for only the optimal source orientation. This is implemented using the cfg.lcmv.fixedori='yes' option.
% Recompute the spatial filter for the optimal source orientation and using that spatial filter (a 1151 vector) recompute the time-series.
% 
% Investigate and describe the difference between the two time-series. What is the difference between the two dipole orientations?
% 
% Note that one orientation is represented in the SVD matrix "u" and the other is in the source.avg.ori field.
% 
% 
% Combine the virtual channel with the EMG
% The raw data structure containing one (virtual) channel can be combined with the two EMG channels from the original preprocessed data.
% 
%% select the two EMG channels
cfg = [];
cfg.channel = 'EMG';
emgdata = ft_selectdata(cfg, data);

%% combine the virtual channel with the two EMG channels
cfg = [];
combineddata = ft_appenddata(cfg, virtualchanneldata, emgdata);

% save combineddata combineddata
% 
% Compute the connectivity
% The resulting combined data structure has three channels: the activity from the cortex, the left EMG and the right EMG. We can now continue with regular channel-level connectivity analysis.
% 
%% compute the spectral decomposition
cfg            = [];
cfg.output     = 'fourier';
cfg.method     = 'mtmfft';
cfg.foilim     = [5 100];
cfg.tapsmofrq  = 5;
cfg.keeptrials = 'yes';
cfg.channel    = {'cortex' 'EMGlft' 'EMGrgt'};
freq    = ft_freqanalysis(cfg, combineddata);

cfg = [];
cfg.method = 'coh';
coherence = ft_connectivityanalysis(cfg, freq);
% This computes the spectral decomposition and the coherence spectrum between all channel pairs, which can be plotted with
% 
cfg = [];
cfg.zlim = [0 0.2];
figure
ft_connectivityplot(cfg, coherence);
title('coherence')

% 
% To look in more detail into the numerical representation of the coherence results, you can use
% 
figure
plot(coherence.freq, squeeze(coherence.cohspctrm(1,2,:)))
title(sprintf('connectivity between %s and %s', coherence.label{1}, coherence.label{2}));
xlabel('freq (Hz)')
ylabel('coherence')

% 
% The spectrum reveals coherence peaks at 10 and 20 Hz (remember that the initial DICS localizer was done at beta). Furthermore, there is a broader plateau of coherence in the gamma range from 40-50 Hz.
% 
% The spectral decomposition was performed with mutitapering and 5 Hz spectral smoothing (i.e. 5Hz in both directions). Recompute the spectral decomposition and the coherence with a hanning taper. Recompute it with mutitapering and 10 Hz smoothing. Plot the three coherence spectra and look at the differences.
% Rather than looking at undirected coherence, the virtual channel level data can now also easily be submitted to directed connectivity measures. Compute the spectrally resolved granger connectivity and try to assess whether the directionality is from cortex to EMG or vice versa.
% 
% Summary and further reading
% This tutorial demonstrates how to compute connectivity measures between two time series. If you want to learn how to make a distributed representation of connectivity throughout the whole brain, you may want to continue with the corticomuscular coherence tutorial.
% 
% This tutorial was last tested by Robert with revision 6026 of FieldTrip (~20120611) on a 64-bit Mac OS X machine using MATLAB 2011b.






%% 
cfg = [];
cfg.method = 'mtmfft';
cfg.taper = 'hanning';
cfg.foi = 1:1:40; % 1-40 Hz
cfg.tapsmofrq = 1; % 1-Hz half-band
cfg.pad = 4; % pad to 4 s
cfg.output = 'fourier';
cfg.keeptrials = 'yes';
freq = ft_freqanalysis(cfg, data); % spectral decomposition

cfgp = [];
cfgp.method = 'psi';
cfgp.bandwidth = 5;
psi_enc = ft_connectivityanalysis(cfgp, freq);



cfg           = [];
cfg.parameter = 'psispctrm';
cfg.zlim      = [0 1];
figure()
ft_connectivityplot(cfg, psi_enc);





















%%
