% combine_figs.m
% Combines restart and non-restart figs together for comparison

% Open both figures
fig1 = open('full_plot.fig');
ax1 = fig1.CurrentAxes;
fig2 = open('restart_plot.fig');
ax2 = fig2.CurrentAxes;

% Take the Davies and Hanson lines from Figure 1
restart1_line = ax1.Children(2);
restart2_line = ax1.Children(1);

% Take the Urban line from Figure 2
full1_line = ax2.Children(2);
full2_line = ax2.Children(1);

% Extract data
full1_x = full1_line.XData;
full1_y = full1_line.YData;
full2_x = full2_line.XData;
full2_y = full2_line.YData;
restart1_x = restart1_line.XData;
restart1_y = restart1_line.YData;
restart2_x = restart2_line.XData;
restart2_y = restart2_line.YData;

% Combine plots together
fig3 = figure;
hold on;
plot(full1_x, full1_y, '-', 'LineWidth', 4);
plot(full2_x, full2_y, '-', 'LineWidth', 4);
plot(restart1_x, restart1_y, '--', 'LineWidth', 4);
plot(restart2_x, restart2_y, '--', 'LineWidth', 4, ...
    'Color', [0.4660, 0.6740, 0.1880]);

% Figure formatting
xlabel('Time [s]');
ylabel('Temperature [eV]');
ax = gca;
ax.FontSize = 16;
xlims = ax.XLim;
ylims = ax.YLim;
grid on

legend('Ti w/ res.', 'Te w/ res.', 'Ti w/o res.', 'Te w/o res.', ...
    'Location', 'SouthEast');

% Save result
saveas(fig3, 'restart_benchmark');
saveas(fig3, 'restart_benchmark','png');


