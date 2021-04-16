% combine_figs.m
% Does final formatting for the Hanson test (combines the Davies_plot.fig
% and the Urban_plot.fig files)

% Open both figures
fig1 = open('Davies_plot.fig');
ax1 = fig1.CurrentAxes;
fig2 = open('Urban_plot.fig');
ax2 = fig2.CurrentAxes;

% Take the Davies and Hanson lines from Figure 1
Davies_line = ax1.Children(2);
Lockwood_line = ax1.Children(1);

% Take the Urban line from Figure 2
Urban_line = ax2.Children(2);

% Extract data
Urb_x = Urban_line.XData;
Urb_y = Urban_line.YData;
Dav_x = Davies_line.XData;
Dav_y = Davies_line.YData;
Loc_x = Lockwood_line.XData;
Loc_y = Lockwood_line.YData;

% Combine plots together
fig3 = figure;
hold on;
plot(Dav_x, Dav_y, 'LineWidth', 2);
plot(Urb_x, Urb_y, 'LineWidth', 2);
plot(Loc_x, Loc_y, ' x', 'LineWidth', 2);

% Figure formatting
ax = gca;
ax.FontSize = 16;
grid on
xlabel('x [CSDA range (195\mum)]');
ylabel('E_{dep} per e^- [MeVcm^{2}/g]');
legend('EPOCH-Davies','EPOCH-Urban','Lockwood (1980)');

% Save result
saveas(fig3, 'Lockwood_benchmark');
saveas(fig3, 'Lockwood_benchmark','png');


