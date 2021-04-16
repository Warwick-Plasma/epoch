% plot_Evans.m
% Performs post processing on EPOCH data to reproduce Fig 4 in Evans
% (2005), and also makes a heatmap to show qualitative similarities to 
% Fig 3. The bulk of the heating in this simulation is Ohmic heating, so
% this provides a benchmark for the hybrid field solver and Ohmic heating
% routines.

% 4 errorbar limits from Evans Fig 4, obtained using WebPlotDigitizer
Evans_raw = [5.089928057553962, 750.8064516129031
             4.982014388489212, 600
             8.974820143884894, 550
             8.974820143884894, 452.41935483870964
             14.046762589928061, 549.1935483870966
             14.046762589928068, 401.61290322580624
             25.971223021582734, 201.61290322580635
             25.917266187050355, 251.61290322580635];

% Process into depth [um], temp [eV], errorbar [eV] data
Evans_depth = [mean(Evans_raw(1:2,1)), mean(Evans_raw(3:4,1)), ...
    mean(Evans_raw(5:6,1)), mean(Evans_raw(7:8,1))];
Evans_temp = [mean(Evans_raw(1:2,2)), mean(Evans_raw(3:4,2)), ...
    mean(Evans_raw(5:6,2)), mean(Evans_raw(7:8,2))];
Evans_etemp = 0.5*(Evans_raw(1:2:7,2)-Evans_raw(2:2:8,2));

%% Add plotting scripts to path

% If MATLAB does not recognise GetDataSDF as a function, add plot epoch to
% the path
if (exist('GetDataSDF') ~=2)
    current_dir = pwd;
    % Keep going up the tree until we find plot_epoch
    for i = 1:10
        cd('..');
        if (exist('plot_epoch') == 7)
            cd('plot_epoch');
            % Add relevant plotting scripts to path
            begin_plotting;
            break
        end
    end 
    cd(current_dir);
end

% Physical constants
c = 299792458;
me = 9.10938356e-31;
q0 = 1.60217662e-19;

%% Get field values

% Ignore the user-input section of get_field
data = GetVariableSDF('0032.sdf','Hybrid.Electron_temperature.data');
x_grid = GetVariableSDF('0032.sdf','Grid.Grid.x');
y_grid = GetVariableSDF('0032.sdf','Grid.Grid.y');
z_grid = GetVariableSDF('0032.sdf','Grid.Grid.z');

% If I forget to dump the grid (which I did the first time I ran this),
% manually enter the grid variables which match the original run
if x_grid(1) == 'C'
    x_grid = 0 : 0.5e-6 : 32.2e-6;
    y_grid = -10.0e-6 : 0.5e-6 : 10.0e-6;
    z_grid = -10.0e-6 : 0.5e-6 : 10.0e-6;
end

%% Average temperature of central bins

% Average central 11 bins - between -0.5um and 0.5um in y, get temp [eV]
kb = 1.38064852e-23;
q0 = 1.60217662e-19;
epoch_temp = mean(mean(data(:,end/2-5:end/2+5,end/2-5:end/2+5),3),2) ...
    * kb / q0; 

% Convert cell edges to cell centres
x_vals = 0.5*(x_grid(2:end) + x_grid(1:end-1))/1.0e-6; % [um]
y_vals = 0.5*(y_grid(2:end) + y_grid(1:end-1))/1.0e-6; % [um]
z_vals = 0.5*(z_grid(2:end) + z_grid(1:end-1))/1.0e-6; % [um]

%% Create temperature-depth plot

% EPOCH data
fig1 = figure;
plot(x_vals, epoch_temp, 'LineWidth', 2);
ax = gca;
grid on;
ax.FontSize = 16;
xlabel('Depth [\mum]');
ylabel('Temperature [eV]');

% Overlay Evans data
hold on;
errorbar(Evans_depth, Evans_temp, Evans_etemp, 'LineWidth', 2, ...
    'LineStyle', 'none');
legend('EPOCH','Evans (2005)','Location','southwest');

%% Create heatmap 

% Plot MATLAB data
fig2 = figure;
[x_mesh, y_mesh] = meshgrid(x_vals, y_vals);
surf(x_mesh', y_mesh', data(:,:,end/2)*kb/q0, 'EdgeColor', 'none', ...
    'FaceColor', 'interp');
col = colorbar;
view(2);
xlim([min(min(x_mesh)), max(max(x_mesh))]);
ylim([min(min(y_mesh)), max(max(y_mesh))]);
ax = gca;
ax.FontSize = 16;
xlabel('x [\mum]');
ylabel('y [\mum]');
col.Label.String = 'Temperature [eV]';

%% Save figures

saveas(fig1, 'EPOCH_vs_Evans');
saveas(fig1, 'EPOCH_vs_Evans', 'png');
saveas(fig2, 'Temperature_heatmap');
saveas(fig2, 'Temperature_heatmap', 'png');
