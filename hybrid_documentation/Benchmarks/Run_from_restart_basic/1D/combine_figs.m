% combine_figs.m
% Does final formatting for the Hanson test (combines the Davies_plot.fig
% and the Urban_plot.fig files)

% Open both figures
fig1 = open('full_plot.fig');
ax1 = fig1.CurrentAxes;
fig2 = open('restart_plot.fig');
ax2 = fig2.CurrentAxes;

% Take the Davies and Hanson lines from Figure 1
restart_line = ax1.Children(2);
Hanson_line = ax1.Children(1);

% Take the Urban line from Figure 2
full_line = ax2.Children(2);

% Extract data
full_x = full_line.XData;
full_y = full_line.YData;
restart_x = restart_line.XData;
restart_y = restart_line.YData;
Han_x = Hanson_line.XData;
Han_y = Hanson_line.YData;

% Combine plots together
fig3 = figure;
hold on;
plot(restart_x, restart_y, 'LineWidth', 2);
plot(Han_x, Han_y, ' x', 'LineWidth',2);
plot(full_x, full_y, '--', 'LineWidth', 2);

% Figure formatting
xlabel('\theta [deg]');
ylabel('Frac. scat. / sq. deg. [deg^{-2}]');
ax = gca;
ax.FontSize = 16;
ax.YAxis.Scale = 'Log';
xlims = ax.XLim;
ylims = ax.YLim;
grid on

legend('w/ restart','Hanson (1951)','w/o restart');

% Save result
saveas(fig3, 'restart_benchmark');
saveas(fig3, 'restart_benchmark','png');


