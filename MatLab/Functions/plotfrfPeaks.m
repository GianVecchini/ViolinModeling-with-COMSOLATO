function plotfrfPeaks(freqArray,frf,coherence,fMin,fMax,figNum, trim,locs,pks)

figure(figNum);
% 1) FRF magnitude in dB
subplot(4,1,1)
semilogx(freqArray, 20*log10(frf), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('|FRF| [dB]')
title('FRF Magnitude (Mobility)')
xlim([fMin fMax]);
xline(locs, 'r'); hold on;
if trim == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end
hold on;
for i = 1:length(locs)
    semilogx(locs(i), pks(i), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 6); hold on;
end

% 2) FRF phase
subplot(4,1,2)
semilogx(freqArray, angle(10.^(frf / 20))*180/pi, 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [deg]')
title('FRF Phase')
xlim([fMin fMax]);
xline(locs, 'r'); hold on;
if trim  == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end

% 3) Coherence log 10
subplot(4,1,3)
semilogx(freqArray(2:end), log10(coherence(2:end)), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('\gamma^2')
title('Coherence log_{10}')
ylim([-1 0])
xlim([fMin fMax]);
xline(locs, 'r'); hold on;
if trim  == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end

% 4) Coherence
subplot(4,1,4)
semilogx(freqArray(2:end), coherence(2:end), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('\gamma^2')
title('Coherence')
ylim([0 1])
xlim([fMin fMax]);
xline(locs, 'r'); hold on;
if trim  == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end
end