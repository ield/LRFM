function [] = plot_Transmission_Parameters(sPar, parToPlot, f, ...
    save_path, filename, special)
% This function plots the the M S-parameters of sPar determined in 
% parToPlot, a Mx2 matrix. If the element is special, some extra things are
% done. For instance, a coupler (special = 1)
if nargin == 5
    special = 0;
end

% Plot the selected S-parameters of the device (magnitude)
figure('Color',[1 1 1]);
% It is defined the legend of the plot
legend_plot = cell(size(parToPlot, 1),1);
for ii = 1:size(parToPlot, 1)
    % It is plotted the s-parameter mn (Smn). First, it is selected
    m = parToPlot(ii, 1);
    n = parToPlot(ii, 2);
    phase_s_mn = 20*log10(abs(sPar(m, n, :)));        phase_s_mn = phase_s_mn(:);
    
    plot(f/1e9, phase_s_mn); hold on;
    % Update the legend
    legend_plot{ii}= sprintf('s_%i_%i',m, n);
end

xlim([f(1) f(end)]/1e9);
if special == 0
    ylim([-50 10]);
elseif special == 1 % If it is a -3dB coupler
    ylim([-5 1]);
end

xlabel('Frequency (GHz)');
ylabel('dB');
legend(legend_plot, 'Location', 'southeast');

set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_mag'],'epsc');

% Plot the selected S-parameters of the device (phase)
figure('Color',[1 1 1]);
% The legend is the same as in the former plot
for ii = 1:size(parToPlot, 1)
    % It is plotted the s-parameter mn (Smn). First, it is selected
    m = parToPlot(ii, 1);
    n = parToPlot(ii, 2);
    phase_s_mn = angle(sPar(m, n, :));      
    phase_s_mn = phase_s_mn(:)*180/pi;
    
    plot(f/1e9, phase_s_mn); hold on;
end

xlim([f(1) f(end)]/1e9);
xlabel('Frequency (GHz)');
ylabel('(º)');
legend(legend_plot, 'Location', 'southeast');


set(gcf,'position',[100,100,400,300]);
set(gca,'FontSize',11, 'fontname','Times New Roman');

saveas(gca, [save_path, filename, '_phase'],'epsc');

% Nothing more to do if the device is not special
if special == 0 
    return;
end
% If we are considering a coupler we plot the phase deifference s21 - s31
if special == 1
    figure('Color',[1 1 1]);
    phase_s21 = angle(sPar(2, 1, :));
    phase_s21 = unwrap(phase_s21(:))*180/pi;
    phase_s31 = angle(sPar(3, 1, :));
    phase_s31 = unwrap(phase_s31(:))*180/pi;
    
    plot(f/1e9, phase_s21 - phase_s31); hold on;
    xlim([f(1) f(end)]/1e9);
    xlabel('Frequency (GHz)');
    ylabel('(º)');
    
    set(gcf,'position',[100,100,400,300]);
    set(gca,'FontSize',11, 'fontname','Times New Roman');
    
    saveas(gca, [save_path, filename, '_phaseDiff'],'epsc');

end

end

