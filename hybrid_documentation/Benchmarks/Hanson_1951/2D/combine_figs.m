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
Hanson_line = ax1.Children(1);

% Take the Urban line from Figure 2
Urban_line = ax2.Children(2);

% Extract data
Urb_x = Urban_line.XData;
Urb_y = Urban_line.YData;
Dav_x = Davies_line.XData;
Dav_y = Davies_line.YData;
Han_x = Hanson_line.XData;
Han_y = Hanson_line.YData;

% Combine plots together
fig3 = figure;
hold on;
plot(Dav_x, Dav_y, 'LineWidth', 2);
plot(Urb_x, Urb_y, 'LineWidth', 2);


% Figure formatting
xlabel('\theta [deg]');
ylabel('Frac. scat. / sq. deg. [deg^{-2}]');
ax = gca;
ax.FontSize = 16;
ax.YAxis.Scale = 'Log';
xlims = ax.XLim;
ylims = ax.YLim;
grid on

% Overlay Geant4 result
load('Geant4_Hanson.mat');
plot(G4x,G4y);
plot(Han_x, Han_y, ' x', 'LineWidth',2);
ax.XLim = xlims;
ax.YLim = ylims;
legend('EPOCH-Davies','EPOCH-Urban','Geant4-Urban','Hanson (1951)');

% Save result
saveas(fig3, 'Hanson_benchmark');
saveas(fig3, 'Hanson_benchmark','png');


