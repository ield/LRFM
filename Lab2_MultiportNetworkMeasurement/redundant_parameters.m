function [] = redundant_parameters(s1, s2, f, legend1, legend2, ...
    save_path, filename)
% Function to plot the redundant S-parameters in magnitude and phase (in
% order to see whether they and different), and to plot the relative error
% after calculateing it.
figure('Color',[1 1 1]);

% Plot the magnitude
plot(f/1e9, 20*log10(abs(s1))); hold on;
plot(f/1e9, 20*log10(abs(s2))); hold on;
xlim([f(1) f(end)]/1e9);
ylim([-50 10]);
xlabel('Frequency (GHz)');
ylabel('dB');
legend(legend1, legend2, 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the S-parameters of the device (phase)

figure('Color',[1 1 1]);

plot(f/1e9, angle(s1)); hold on;
plot(f/1e9, angle(s2)); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Phase (º)');
legend(legend1, legend2, 'Location', 'best');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');

% Plot the error
error = abs((s1 - s2)./s1)*100;

figure('Color',[1 1 1]);

plot(f/1e9, error); hold on;
xlabel('Frequency (GHz)');
xlim([f(1) f(end)]/1e9);
ylabel('Error (%)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_error'],'epsc');
end

