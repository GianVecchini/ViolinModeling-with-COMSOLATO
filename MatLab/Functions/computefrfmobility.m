function [frf, coherence, freqArray] = computefrfmobility(measures, nMeasures, nfft, fs)

% This function computes the FRF basing on Noise on the Output Estimator 
% and spectrum averaging procedure subsequent tests on one same 
% measurement set are averaged to obtain a more robust frf.

% INPUT PARAMETERS
%   - measures = the 3 dimensional array containing all the signals data
%   - nMeasures = the number of tests (i.e measures) performed = 5
%   - fs = sampling frequency for the measured data
%   - nfft = the number of points on which compute the fft

% OUTPUT PARAMETERS
%   - frf = the output frf computed
%   - coherence = the coherence of the computed frf
%   - freq = the frequency array of samples, spaced df = fs/ nfft


%% FFT AND PSD
% Preallocate matrixes to store the fft of the input/output signals and
% their power spectral density
fftIn = zeros(nMeasures,nfft);          % fft preallocation (input)
PDSInput = zeros(nMeasures,nfft/2+1);   % PSD preallocation (input)
fftOut = zeros(nMeasures,nfft);         % fft preallocation (input)
PDSOutput = zeros(nMeasures,nfft/2+1);  % PSD preallocation (input)

% Get the spectra and power spectra using the function  computefftandpowerspectra
for i=1:nMeasures
    
    % obs: squeeze transforms a n*1 matrix in a simple array of length n
    % Save the i-th signal (out/input) in the i-th position of the matrixes
    [fftIn(i,:),PDSInput(i,:)] = computefftandpowerspectra(squeeze(measures(i,1,:)), nfft, fs);
    [fftOut(i,:),PDSOutput(i,:)] = computefftandpowerspectra(squeeze(measures(i,2,:)), nfft, fs);
end

% Average over the different measures power spectral densiry
AvgPsdIn = mean(PDSInput,1);   % Average PSD (input)
AvgPsdOut = mean(PDSOutput,1); % Average PSD (output)

% Compute the Cross Spectral Density (CSD) between input and output
% response
crossSpectralDensity = zeros(nMeasures, nfft); % Initialize CSD array
for i = 1:nMeasures
    % Compute CSD for each measurement, output*(conjugate of input), all
    % normalized by N*fs
    crossSpectralDensity(i, :) = (fftOut(i, :) .* conj(fftIn(i, :))) / (length(squeeze(measures(1,1,:))) * fs);
end

% Average over the columns measurements
AvgCSD = mean(crossSpectralDensity, 1);

% Extract single-sided spectrum
singleSideCSD = AvgCSD(1:nfft/2+1);
singleSideCSD(2:end-1) = 2 * singleSideCSD(2:end-1);

% Compute the coherence: cross spectral density^2 / power spectral density
% of input times the one in output
coherence = abs(singleSideCSD).^2 ./ (AvgPsdIn .* AvgPsdOut);

% Noise on the output FRF estimator for mobility
frf = singleSideCSD ./ AvgPsdIn;     % Mobility FRF: output/input

% Integrate acceleration to obtain velocity/force (mobility).
df = fs / nfft;                       % Distance between samples in frequency          
freqArray = (0:nfft/2) * df;          % Frequency axis
frf = frf./(1i*2*pi*freqArray);
