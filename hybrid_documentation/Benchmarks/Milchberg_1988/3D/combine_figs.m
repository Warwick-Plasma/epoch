% combine_figs.m
% Does final formatting for the Milchberg test (combines the with_fit.fig
% and the without_fit.fig files)

% Open both figures
fig1 = open('with_fit.fig');
ax1 = fig1.CurrentAxes;
fig2 = open('without_fit.fig');
ax2 = fig2.CurrentAxes;

% Take the Davies and Hanson lines from Figure 1
with_line = ax1.Children(1);
Milch_data = ax1.Children(2);

% Take the Urban line from Figure 2
without_line = ax2.Children(1);

% Extract data
without_x = without_line.XData;
without_y = without_line.YData;
with_x = with_line.XData;
with_y = with_line.YData;
milch_x = Milch_data.XData;
milch_y = Milch_data.YData;

% Combine plots together
fig3 = figure;
plot(milch_x, milch_y, ' x', 'LineWidth',2);
hold on;
plot(with_x, with_y, 'LineWidth', 2);
plot(without_x, without_y, 'LineWidth', 2);

% Figure formatting
xlabel('Te [eV]');
ylabel('Resistivity [\Omegam]');
ax = gca;
ax.FontSize = 16;
grid on
legend('Milchberg (1988)', '(\lambda_1, \lambda_2)=(7, 3.5)', ...
    '(\lambda_1, \lambda_2)=(1, 1)');

% Save result
saveas(fig3, 'Milchberg_benchmark');
saveas(fig3, 'Milchberg_benchmark','png');


