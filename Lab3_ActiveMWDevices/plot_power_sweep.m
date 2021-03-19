function [] = plot_power_sweep(p, spar_mag, spar_phase, ...
    save_path, filename)
% Function to plot the S-parameter of a device specified as input (in dB
% and degrees)
figure('Color',[1 1 1]);

% In the magnitude, the function also finds the P1dB. This is done before
% plotting
p_1dB = spar_mag(1) - 1;
loc_p1db = find(spar_mag <= p_1dB, 1);

% Plot the sweep of the s21
plot(p, spar_mag); hold on;     
% Plot the line of the P1dB
plot(p, p_1dB*ones(size(spar_mag)), '--'); hold on;
% Plot the intersection point
plot(p(loc_p1db), spar_mag(loc_p1db), 'o');

xlim([p(1) p(end)]);
xlabel('Input power (dBm)');
ylabel('s_{21} (dB)');
legend('s_{21}', 's_{21 max} - 1dB', 'P_{1dB}', 'location', 'southwest');
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the S-parameters of the device (phase)
figure('Color',[1 1 1]);

plot(p, spar_phase); hold on;
xlabel('Input power (dBm)');
xlim([p(1) p(end)]);
ylabel('Phase (º)');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');

% Plot the output power according to the gain and the input power, both
% expected and measured
figure('Color',[1 1 1]);

% Plot the measured output power according to the gain and the input power
plot(p, p + spar_mag); hold on;
% Plot the expected output power according to the gain and the input power
plot(p, p + spar_mag(1), '--'); hold on;
% Plot the P1dB
plot(p(loc_p1db), p(loc_p1db) + spar_mag(loc_p1db), 'o');

xlim([p(1) p(end)]);
xlabel('Input power (dBm)');
ylabel('Output power (dBm)');
legend('P_{out} measured', 'P_{out} with G = cte', 'P_{1dB}', ...
    'location', 'southeast');
set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_out'],'epsc');

end

