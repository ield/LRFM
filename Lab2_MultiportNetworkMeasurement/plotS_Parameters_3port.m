function [] = plotS_Parameters_3port(sPar, f, save_path, filename)
% Function to plot the S-parameters of a device: magnitude and phase
% Plot the S-parameters of the device (magnitude)
s11 = 20*log10(abs(sPar(1, 1, :)));        s11 = s11(:);
s21 = 20*log10(abs(sPar(2, 1, :)));        s21 = s21(:);
s12 = 20*log10(abs(sPar(1, 2, :)));        s12 = s12(:);
s22 = 20*log10(abs(sPar(2, 2, :)));        s22 = s22(:);
s31 = 20*log10(abs(sPar(3, 1, :)));        s31 = s31(:);
s23 = 20*log10(abs(sPar(2, 3, :)));        s23 = s23(:);
s32 = 20*log10(abs(sPar(3, 2, :)));        s32 = s32(:);
s13 = 20*log10(abs(sPar(1, 3, :)));        s13 = s13(:);
s33 = 20*log10(abs(sPar(3, 3, :)));        s33 = s33(:);

figure('Color',[1 1 1]);

plot(f/1e9, s11); hold on;
plot(f/1e9, s12); hold on;
plot(f/1e9, s13); hold on;
plot(f/1e9, s21); hold on;
plot(f/1e9, s22); hold on;
plot(f/1e9, s23); hold on;
plot(f/1e9, s31); hold on;
plot(f/1e9, s32); hold on;
plot(f/1e9, s33); hold on;

xlim([f(1) f(end)]/1e9);
ylim([-50 10]);
xlabel('Frequency (GHz)');
ylabel('dB');
legend('s_{11}', 's_{12}', 's_{13}', 's_{21}', 's_{22}', 's_{23}', ...
    's_{31}', 's_{32}', 's_{33}', 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the S-parameters of the device (phase)
phase_s11 = angle(sPar(1, 1, :));    phase_s11 = phase_s11(:)*180/pi;
phase_s12 = angle(sPar(1, 2, :));    phase_s12 = phase_s12(:)*180/pi;
phase_s13 = angle(sPar(1, 3, :));    phase_s13 = phase_s13(:)*180/pi;
phase_s21 = angle(sPar(2, 1, :));    phase_s21 = phase_s21(:)*180/pi;
phase_s22 = angle(sPar(2, 2, :));    phase_s22 = phase_s22(:)*180/pi;
phase_s23 = angle(sPar(2, 3, :));    phase_s23 = phase_s23(:)*180/pi;
phase_s31 = angle(sPar(3, 1, :));    phase_s31 = phase_s31(:)*180/pi;
phase_s32 = angle(sPar(3, 2, :));    phase_s32 = phase_s32(:)*180/pi;
phase_s33 = angle(sPar(3, 3, :));    phase_s33 = phase_s33(:)*180/pi;

figure('Color',[1 1 1]);

plot(f/1e9, phase_s11); hold on;
plot(f/1e9, phase_s12); hold on;
plot(f/1e9, phase_s13); hold on;
plot(f/1e9, phase_s21); hold on;
plot(f/1e9, phase_s22); hold on;
plot(f/1e9, phase_s23); hold on;
plot(f/1e9, phase_s31); hold on;
plot(f/1e9, phase_s32); hold on;
plot(f/1e9, phase_s33); hold on;

xlim([f(1) f(end)]/1e9);
xlabel('Frequency (GHz)');
ylabel('(º)');
legend('s_{11}', 's_{12}', 's_{13}', 's_{21}', 's_{22}', 's_{23}', ...
    's_{31}', 's_{32}', 's_{33}', 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');

end

