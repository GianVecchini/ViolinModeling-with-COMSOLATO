%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                        %
%  VIOLIN PROJECT - EXPERIMENTAL MATLAB  %
%                                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions');

fs = 96000;
nfft = 8192;

%% Extract the file names and the number of them
dataPath='Data/violinFRF/hh';
[sigList,nSig]=extractfoldercontent(dataPath);

%% Extract the hammer and accelerometer signals
% We have 5 files, each file contains 2 signal, first the hammeR and second
% the accelerometer response
for i=1:nSig

    % Create the file path
    dataName = strcat(dataPath, '/', sigList(i)); 

    % Save the signals in data (thanks to the txt2mat function)
    data = txt2mat(dataName{1,1});
    
    % Preallocate the 3d matrix to save data:
    % - x: number of signal (1,...,5)
    % - y: save hammer and response ( 1 and 2)
    % - z: signal (length of the signal)
    if i == 1
        measurments = zeros(nSig,2,length(data(:,1)));
    end

    % Save in x = i the i-th signals
    measurments(i,1,:) = data(:,1); % Save in y = 1 the hammer
    measurments(i,2,:) = data(:,2); % Save in y = 2 the response
end

%% Compute the Average frf, coherence and freqArray
[frf, coherence, freqArray] = computefrfmobility(measurments,nSig,nfft,fs);

%% Plot results

% TRIM FRF MEASURES in the valid interval
fMin = 200;     % Hz
fMax = 6000;    % Hz

plotfrf(freqArray,frf,coherence,fMin,fMax,1,0)

% TRIM FRF MEASURES
idxTrim = (freqArray >= fMin) & (freqArray <= fMax);

freq_trim = freqArray(idxTrim);
frf_trim  = frf(idxTrim);
coh_trim  = coherence(idxTrim);

plotfrf(freq_trim,frf_trim,coh_trim,fMin,fMax,2,1)

% Interpolate the frf and the freq array
f_fine = freq_trim(1):0.1:freq_trim(end); 
frf_intr = interp1(freq_trim, frf_trim, f_fine);
intr_coh = interp1(freq_trim, coh_trim, f_fine);

%% FIND PEAKS
% Two methods:
% - the first just chooses the higher peaks
% - in the second method it chooses the higher peaks but also some peaks
% are added by hands
method = 2;

% Find peaks METHOD 1
if method == 1
    frf_linear = abs(frf_trim);

    % Find al peaks
    [cmsl_vals, cmsl_peaks, ~] = findpeaks(frf_linear', freq_trim', ...
        'Annotate', 'extents', ...
        'WidthReference', 'halfheight');

    % Find the 3db band width
    [~, ~, cmsl_bw] = findpeaks(abs(frf_trim).^2, freq_trim, 'WidthReference', 'halfheight');
    
    % Order based on the higher value
    [~, indexSorted] = sort(cmsl_vals, 'descend');
    
    % Select top N peaks
    npeaks = 10;
    target_indices = indexSorted(1:min(npeaks, length(indexSorted)));
    
    locs = cmsl_peaks(target_indices);  % Central Freqs [Hz]
    pks = cmsl_vals(target_indices);    % Real magnitude [m/s/N]
    bw = cmsl_bw(target_indices);       % band width [Hz]

    plotfrfPeaks(f_fine,frf_intr,intr_coh,fMin,fMax,3,1,locs,20*log10(pks))
end

% Find peaks METHOD 2
if method == 2
    frf_linear = abs(frf_trim); 
    
    [cmsl_vals, cmsl_peaks, ~] = findpeaks(frf_linear', freq_trim', ...
        'Annotate', 'extents', ...
        'WidthReference', 'halfheight');
    [~, ~, cmsl_bw] = findpeaks(abs(frf_trim).^2, freq_trim, 'WidthReference', 'halfheight');
    
    % Select top N peaks
    npeaks = 7;
    [~, indexSorted] = sort(cmsl_vals, 'descend');
    top_indices = indexSorted(1:min(npeaks, length(indexSorted)));
    top_indices = top_indices(:);
    
    % Select manually some peaks 
    manual_indices = [1, 2, 8, 9, 14, 18, 21, 24, 27, 41];
    % Filter them
    manual_indices = manual_indices(manual_indices <= length(cmsl_peaks));
    manual_indices = manual_indices(:);
    
    % Manually add the new peaks
    target_indices = unique([top_indices; manual_indices], 'stable');
    
    % Save data
    locs = cmsl_peaks(target_indices);
    pks = cmsl_vals(target_indices);   
    bw = cmsl_bw(target_indices);      
    
    % Plot
    plotfrfPeaks(f_fine, frf_intr, intr_coh, fMin, fMax, 3, 1, locs, 20*log10(pks))
end

%% CREATE FILTERBANK
nPointsOfMeasure = 72000;
impulse_t = 0:1/fs:(nPointsOfMeasure - 1)/fs;
impulse_s = double(impulse_t==0);

[simResponse, simulatedTime] = makefilterbank(locs, bw , pks, length(locs), nPointsOfMeasure, fs, impulse_s, impulse_t);

%% Plot
n_freqs = length(simResponse);
simfft = fft(simResponse, n_freqs);
simfft = simfft(1:length(simfft)/2);
n_freqs = n_freqs/2;
simFreq = linspace(0, fs/2, n_freqs);

fig = figure(4);
subplot(2,1,1)
semilogx(f_fine, 20*log10(abs(frf_intr)), 'LineWidth',2); hold on;
semilogx(simFreq, 20*log10(abs(simfft)), 'LineWidth',2); hold on;
plot(locs(:), 20*log10(pks(:)), 'o');
xline(locs, 'r');
xlim([fMin fMax])
grid on
hold off
title('measure VS filterbank');
xlabel('Frequency (Hz)');
ylabel('|FRF| [dB]');
legend('measured mobility', 'filterbank', 'Location','best');

subplot(2,1,2)
semilogx(f_fine,angle(frf_intr), 'LineWidth',2); hold on;
semilogx(simFreq, angle(simfft), 'LineWidth',2); hold on;
xline(locs, 'r');
xlim([fMin fMax]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});

grid on
hold off
ylabel('Phase [rad]');
xlabel('Frequency (Hz)');
exportgraphics(fig, 'Plots\LabMobilityRLC.png', 'Resolution', 300);


%% KARPLUS STRONG
[frequencies] = karplusStrong();
[simResponse_ks, simTime_ks] = makefilterbank(locs, bw , pks, length(locs), nPointsOfMeasure, fs, y_ks, t_ks);

%% play and plot result
sound(simResponse_ks, fs);
n_freqs_ks = length(simResponse_ks);
simfft_ks = fft(simResponse_ks, n_freqs_ks);
simfft_ks = simfft_ks(1:length(simfft_ks)/2);
n_freqs_ks = n_freqs_ks/2;
simFreq_ks = linspace(0, fs/2, n_freqs_ks);

fig = figure
semilogx(frequencies, 20*log10(abs(magnitude)), 'LineWidth',0.6,'Color', 'b'); hold on
semilogx(simFreq_ks, 20*log10(abs(simfft_ks./40)), 'LineWidth',0.8,'Color', 'g'); hold on;
%semilogx(f_fine, 20*log10(abs(frf_intr)), 'LineWidth',2); hold on;
semilogx(simFreq, 20*log10(abs(simfft)), 'LineWidth',2,'Color', 'k'); hold on;
xlim([100 5000]);
grid on
ylabel('Magnitude [dB)]');
xlabel('Frequencies [Hz]');
title('Lab Filterbank - Karplus Strong');
legend('Karplus-Strong input signal', 'Filterbank output', 'measured FRF', 'Filterbank response');
exportgraphics(fig, 'Plots\LabFilterbankOutput.png', 'Resolution', 300);

%% audiowrite
audiowrite("Sounds\LabFilterBank Output.wav", simResponse_ks(1:191995)*1.5, fs);
