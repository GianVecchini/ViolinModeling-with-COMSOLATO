%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                      %
%  VIOLIN PROJECT - SIMULATION MATLAB  %
%                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

addpath('Functions');

fs = 96000;
nfft = 8192;

%% Extract the file names and the number of them
filename = 'Data/ComsolMobility.csv';
data = readmatrix(filename); 

freq_comsol = data(:, 1); 
mobility_comsol = data(:, 2);

%% 2. TRIMMING (Definisci il range di interesse come nel Main precedente)

fig = figure(1);
% 1) FRF magnitude in dB
subplot(2,1,1)
semilogx(freq_comsol, 20*log10(mobility_comsol), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('|Mobility| [dB]')
title('Mobility Magnitude')
xlim([100 10000]);

% 2) FRF phase
subplot(2,1,2)
semilogx(freq_comsol, angle(mobility_comsol)*180/pi, 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
title('Magnitude Phase')
yticks([-180, -90, 0, 90, 180]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xlim([100 10000]);
%exportgraphics(fig, 'Plots\ComsolMobity.png', 'Resolution', 300);

%% TRIM FRF MEASURES
fMin = 200;     % Hz
fMax = 6000;    % Hz

idxTrim = (freq_comsol >= fMin) & (freq_comsol <= fMax);

freq_trim = freq_comsol(idxTrim);
frf_trim  = mobility_comsol(idxTrim);

% Interpolate
f_fine = freq_trim(1):0.1:freq_trim(end); 
frf_intr = interp1(freq_trim, frf_trim, f_fine);

%% FIND PEAKS
% Two methods:
% - the first just chooses the higher peaks
% - in the second method it chooses the higher peaks but also some peaks
% are added by hands
method = 2;

% Find peaks METHOD 1
if method == 1
    frf_linear = abs(frf_intr);

    % Find al peaks
    [cmsl_vals, cmsl_peaks, ~] = findpeaks(frf_linear', f_fine', ...
        'Annotate', 'extents', ...
        'WidthReference', 'halfheight');

    % Find the 3db band width
    [~, ~, cmsl_bw] = findpeaks(abs(frf_linear).^2, f_fine, 'WidthReference', 'halfheight');
    
    % Order based on the higher value
    [~, indexSorted] = sort(cmsl_vals, 'descend');
    
    % Select top N peaks
    npeaks = 30;
    target_indices = indexSorted(1:min(npeaks, length(indexSorted)));
    
    locs = cmsl_peaks(target_indices);  % Central Freqs [Hz]
    pks = cmsl_vals(target_indices);    % Real magnitude [m/s/N]
    bw = cmsl_bw(target_indices);       % band width [Hz]
end

% Find peaks METHOD 2
if method == 2
    frf_linear = abs(frf_trim); 
    
    [cmsl_vals, cmsl_peaks, ~] = findpeaks(frf_linear', freq_trim', ...
        'Annotate', 'extents', ...
        'WidthReference', 'halfheight');
    [~, ~, cmsl_bw] = findpeaks(abs(frf_trim).^2, freq_trim, 'WidthReference', 'halfheight');
    
    % Select top N peaks
    npeaks = 10;
    [~, indexSorted] = sort(cmsl_vals, 'descend');
    top_indices = indexSorted(1:min(npeaks, length(indexSorted)));
    top_indices = top_indices(:);
    
    % Select manually some peaks 
    manual_indices = [2, 4, 6, 7, 9, 10, 11, 17, 18, 24, 37,70];
    % Filter them
    manual_indices = manual_indices(manual_indices <= length(cmsl_peaks));
    manual_indices = manual_indices(:);
    
    % Manually add the new peaks
    target_indices = unique([top_indices; manual_indices], 'stable');
    
    % Save data
    locs = cmsl_peaks(target_indices);
    pks = cmsl_vals(target_indices);   
    bw = cmsl_bw(target_indices);      
end

% Plot
figure(2);
% 1) FRF magnitude in dB
subplot(2,1,1)
semilogx(f_fine, 20*log10(frf_intr), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('|FRF| [dB]')
title('FRF Magnitude (Mobility)')
xlim([fMin fMax]);
xline(locs, 'r'); hold on;
hold on;
for i = 1:length(locs)
    semilogx(locs(i), 20*log10(pks(i)), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
end

% 2) FRF phase
subplot(2,1,2)
semilogx(f_fine, angle(frf_intr)*180/pi, 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
title('FRF Phase')
xlim([fMin fMax]);
yticks([-180, -90, 0, 90, 180]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
xline(locs, 'r'); hold on;

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

fig=figure(4);
subplot(2,1,1)
semilogx(f_fine, 20*log10(abs(frf_intr)), 'LineWidth',2); hold on;
semilogx(simFreq, 20*log10(abs(simfft)), 'LineWidth',2); hold on;
plot(locs(:), 20*log10(pks(:)), 'o');
xline(locs, 'r');
xlim([fMin fMax])
grid on
hold off
title('measure VS filterbank');
ylabel('|FRF| [dB]');
xlabel('Frequency (Hz)');
legend('measured mobility', 'filterbank', 'Location', 'best');

subplot(2,1,2)
semilogx(f_fine,angle(frf_intr), 'LineWidth',2); hold on;
semilogx(simFreq, angle(simfft), 'LineWidth',2); hold on;
xline(locs, 'r');
xlim([fMin fMax]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
grid on
hold off
ylabel('Phase [rad]')
xlabel('Frequency (Hz)');
exportgraphics(fig, 'Plots\ComsolMobilityRLC.png', 'Resolution', 300);

%% KARPLUS STRONG
[frequencies, magnitude] = karplusStrong();
[simResponse_ks, simTime_ks] = makefilterbank(locs, bw , pks, length(locs), nPointsOfMeasure, fs, y_ks, t_ks);

%% play and plot results (measure)
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
title('Comsol Filterbank - Karplus Strong');
%legend('Karplus-Strong input signal', 'Filterbank output', 'measured FRF', 'Filterbank response');
legend('Karplus-Strong input signal', 'Filterbank output', 'Filterbank response', 'Location','Best');
exportgraphics(fig, 'Plots\ComsolFilterbankOutput.png', 'Resolution', 300);

%% audiowrite
audiowrite("Sounds\ComsolFilterBank Output.wav", simResponse_ks(1:191996)*1.5, fs);
