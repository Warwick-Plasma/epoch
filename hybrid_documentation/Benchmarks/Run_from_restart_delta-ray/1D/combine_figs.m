% combine_figs.m
% Combines all delta-ray plots

% Open all
fig1 = open('restart_plot1.fig');
ax1 = fig1.CurrentAxes;
fig2 = open('restart_plot2.fig');
ax2 = fig2.CurrentAxes;
fig3 = open('full_plot1.fig');
ax3 = fig3.CurrentAxes;
fig4 = open('full_plot2.fig');
ax4 = fig4.CurrentAxes;

% Take lines from plots
restart1_line = ax1.Children(1);
restart2_line = ax2.Children(1);
full1_line = ax3.Children(1);
full2_line = ax4.Children(1);

% Extract data
full1_x = full1_line.XData;
full1_y = full1_line.YData;
full2_x = full2_line.XData;
full2_y = full2_line.YData;
restart1_x = restart1_line.XData;
restart1_y = restart1_line.YData;
restart2_x = restart2_line.XData;
restart2_y = restart2_line.YData;

% Combine plots together (energy spectrum)
fig5 = figure;
hold on;
plot(full1_x, full1_y, '-', 'LineWidth', 4);
plot(restart1_x, restart1_y, '--', 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('KE [MeV]');
ylabel('dN/dKE [MeV^{-1}]');
legend('w/ res.', 'w/o res.', 'Location', 'northwest');

% Save result
saveas(fig5, 'restart_benchmark1');
saveas(fig5, 'restart_benchmark1','png');

% Combine plots together (angular distribution)
fig6 = figure;
hold on;
plot(full2_x, full2_y, '-', 'LineWidth', 4);
plot(restart2_x, restart2_y, '--', 'LineWidth', 2, ...
    'Color', [0.9290, 0.6940, 0.1250]);
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('\theta [rad]');
ylabel('dN/d\theta [rad^{-1}]');
legend('w/ res.', 'w/o res.', 'Location', 'northeast');

% Save result
saveas(fig6, 'restart_benchmark2');
saveas(fig6, 'restart_benchmark2','png');
