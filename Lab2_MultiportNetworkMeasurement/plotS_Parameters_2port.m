function [] = plotS_Parameters_2port(sPar, f, save_path, filename)
% Function to plot the S-parameters of a device: magnitude and phase
% Plot the S-parameters of the device (magnitude)
s11 = 20*log10(abs(sPar(1, 1, :)));        s11 = s11(:);
s21 = 20*log10(abs(sPar(2, 1, :)));        s21 = s21(:);
s12 = 20*log10(abs(sPar(1, 2, :)));        s12 = s12(:);
s22 = 20*log10(abs(sPar(2, 2, :)));        s22 = s22(:);

figure('Color',[1 1 1]);

plot(f/1e9, s11); hold on;
plot(f/1e9, s12); hold on;
plot(f/1e9, s21); hold on;
plot(f/1e9, s22); hold on;
xlim([f(1) f(end)]/1e9);
ylim([-50 10]);
xlabel('Frequency (GHz)');
ylabel('dB');
legend('s_{11}', 's_{12}', 's_{21}', 's_{22}', 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the S-parameters of the device (phase)
phase_s11 = angle(sPar(1, 1, :));    phase_s11 = phase_s11(:)*180/pi;
phase_s12 = angle(sPar(1, 2, :));    phase_s12 = phase_s12(:)*180/pi;
phase_s21 = angle(sPar(2, 1, :));    phase_s21 = phase_s21(:)*180/pi;
phase_s22 = angle(sPar(2, 2, :));    phase_s22 = phase_s22(:)*180/pi;

figure('Color',[1 1 1]);

plot(f/1e9, phase_s11); hold on;
plot(f/1e9, phase_s12); hold on;
plot(f/1e9, phase_s21); hold on;
plot(f/1e9, phase_s22); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Phase (º)');
legend('s_{11}', 's_{12}', 's_{21}', 's_{22}', 'Location', 'best');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');
end

