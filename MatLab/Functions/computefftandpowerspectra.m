function [fftSignal,powerDensitySpectrum] = computefftandpowerspectra(signal, nfft, fs)
    
    % Simple fft with nfft freq bins
    fftSignal = fft(signal,nfft);
    
    % Power spectral density
    doubleSidedPSD = (abs(fftSignal).^2) / (length(signal)*fs);

    % consider a single side and molitply by to the freq that are not 0 and
    % nyquest to preserve the energy
    singleSidedPSD = doubleSidedPSD(1:nfft/2+1);
    singleSidedPSD(2:end-1) = singleSidedPSD(2:end-1) * 2;
    powerDensitySpectrum = singleSidedPSD;
end