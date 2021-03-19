function [] = plotS_Parameters_2port(f, s11_mod, s11_phase, s21_mod, ...
    s21_phase, s12_mod, s12_phase, s22_mod, s22_phase, ymin, ymax, ...
    save_path, filename)

% Function to plot the S-parameters of a device: magnitude and phase
% Plot the S-parameters of the device (magnitude)
figure('Color',[1 1 1]);

plot(f/1e9, s11_mod); hold on;
plot(f/1e9, s12_mod); hold on;
plot(f/1e9, s21_mod); hold on;
plot(f/1e9, s22_mod); hold on;
xlim([f(1) f(end)]/1e9);
ylim([ymin ymax]);
xlabel('Frequency (GHz)');
ylabel('dB');
legend('s_{11}', 's_{12}', 's_{21}', 's_{22}', 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the S-parameters of the device (phase)
figure('Color',[1 1 1]);

plot(f/1e9, s11_phase); hold on;
plot(f/1e9, s12_phase); hold on;
plot(f/1e9, s21_phase); hold on;
plot(f/1e9, s22_phase); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Phase (º)');
legend('s_{11}', 's_{12}', 's_{21}', 's_{22}', 'Location', 'best');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');
end

