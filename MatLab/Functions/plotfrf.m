function plotfrf(freqArray,frf,coherence,fMin,fMax,figNum, trim)

fig = figure(figNum);
% 1) FRF magnitude in dB
subplot(3,1,1)
semilogx(freqArray, 20*log10(abs(frf)), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('|FRF| [dB]')
title('Mobility Magnitude')
xlim([fMin fMax]);
if trim == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end

% 2) FRF phase
subplot(3,1,2)
semilogx(freqArray, angle(frf), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('Phase [rad]')
title('Mobility Phase')
xlim([fMin fMax]);
yticks([-pi, -pi/2, 0, pi/2, pi]);
yticklabels({'-\pi', '-\pi/2', '0', '\pi/2', '\pi'});
if trim  == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end

% 3) Coherence log 10
% subplot(4,1,3)
% semilogx(freqArray(2:end), log10(coherence(2:end)), 'LineWidth', 1.2)
% grid on
% xlabel('Frequency [Hz]')
% ylabel('\gamma^2')
% title('Coherence log_{10}')
% ylim([-1 0])
% xlim([fMin fMax]);
% if trim  == 0
%     xlim([100 15000]);
%     xline(fMin, '--', 'Color','r','LineWidth',2);
%     xline(fMax, '--', 'Color','r','LineWidth',2);
% end

% 4) Coherence
subplot(3,1,3)
semilogx(freqArray(2:end), coherence(2:end), 'LineWidth', 1.2)
grid on
xlabel('Frequency [Hz]')
ylabel('\gamma^2')
title('Coherence')
ylim([0 1])
 xlim([fMin fMax]);
if trim  == 0
    xlim([100 15000]);
    xline(fMin, '--', 'Color','r','LineWidth',2);
    xline(fMax, '--', 'Color','r','LineWidth',2);
end
exportgraphics(fig, 'Plots\LabMobity.png', 'Resolution', 300);
end