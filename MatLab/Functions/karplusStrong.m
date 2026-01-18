function [frequencies, magnitude] = karplusStrong()
    fs = 96000;        
    f0 = 440; 
    g = 0.99;
    
    y = ks_lpcomb(f0, g, fs);
    t = 0:1/fs:(length(y)-1)/fs;
    
    % Plot
    fig = figure;
    plot(t, y)
    title('Output Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    exportgraphics(fig, 'Plots\OutputSig.png', 'Resolution', 300);
    
    % Play Output
    sound(y, fs);
    audiowrite("Sounds\KarplusStrong Output.wav", y, fs);


    N = length(y);
    Y = fft(y);
    magnitude = abs(Y) / N; 
    magnitude = magnitude(1:N/2+1); 
    
    % Frequencies
    frequencies = (0:N/2) * (fs / N);
    
    % Plot
    fig = figure;
    plot(frequencies, 20*log10(magnitude)); % Magnitudine in dB
    title('Output Signal Spectrum');
    xlabel('Frequenza (Hz)');
    xlim([0 10000]);
    ylabel('Magnitudine (dB)');
    grid on;
    exportgraphics(fig, 'Plots\OutputMag.png', 'Resolution', 300);
    
    assignin('base', 'y_ks', y );
    assignin('base', 't_ks', t );
end

function y = ks_lpcomb(f0, g, fs)
    % Compute delay length in samples (integer delay)
    % Approx pitch: f0 ≈ fs/m
    m = round(fs/f0);  % legnth of delay line
    disp("m = " + m);
    disp("Actual delay(fs/f0) = " + fs/f0 );

    % Initialize delay line object (struct with buffer + input/output samples)
    d = dline_init(m); % create the delay line object

    % Create excitation: random noise (pluck) (random)
    % rand in [0,1] -> scale to [-1,1] -> then halve amplitude
    x = (rand(1, fs)*2 - 1)/2;

    % low-pass filter of the white noise
    [b, a] = butter(1, 2400/(fs/2));
    x = filter(b, a, x);
    audiowrite("Sounds\KarplusStrong Input.wav", x, fs);

    % Zero-padding
    x = [x zeros(1,fs*1)];
    t1 = 0:1/fs:(length(x)-1)/fs;

    % Plot
    fig = figure;
    plot(t1,x);
    title('Input Signal');
    xlabel('Time (s)');
    ylabel('Amplitude');
    exportgraphics(fig, 'Plots\InputSig.png', 'Resolution', 300);

    Nx = length(x);
    X = fft(x);
    magnitudex = abs(X) / Nx; 
    magnitudex = magnitudex(1:Nx/2+1); 
    
    % frequencies
    frequenciesx = (0:Nx/2) * (fs / Nx);
    
    % Plot
    fig = figure;
    plot(frequenciesx, 20*log10(magnitudex)); % Magnitudine in dB
    title('Input Signal Spectrum');
    xlabel('Frequency (Hz)');
    xlim([0 10000]);
    ylabel('Magnitudine (dB)');
    grid on;
    exportgraphics(fig, 'Plots\InputMag.png', 'Resolution', 300);

    % Allocate output signal
    y = zeros(1, length(x));

    % Store previous loop value for the 2-point averager
    % (implements: y(n) = 0.5*(a(n) + a(n-1)))
    a_past = 0; %initialize auxiliary variable

    for n = 1:length(x)
        a = x(n) + g*d.y;
        y(n) = 0.5*(a + a_past);
        a_past = a;

        d.x = y(n);
        d = dline_compute(d);
    end
end


function f = dline_init(d)
    f.x = 0;
    f.y = 0;
    if(floor(d) == d)
        f.d = d;
    else
        error("Not valid delay")
    end
    f.in = zeros(1,d);
end

function f = dline_compute(f)
f.y = f.in(1);
f.in = [f.in(2:length(f.in)), f.x];
end