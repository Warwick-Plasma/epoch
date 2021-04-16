% plot_lockwood.m
% Performs post processing of EPOCH results to compare energy loss as a
% function of depth with the plot in Lockwood 1980 (report pg. 100 & 101,
% PDF pg. 92, 93). This is a depth/dose curve for 0.5 MeV electrons in Ta

% Lockwood data (published in report)
Ta_500kev_r = [0.034, 0.07, 0.1, 0.144, 0.193, 0.253, 0.333, 0.488];
Ta_500kev_E = [3.9, 4.23, 4.02, 3.35, 2.61, 1.61, 0.71, 0.1];


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

%% Get field arrays

temperature = GetVariableSDF('0001.sdf', 'Hybrid.Electron_temperature.data');
x_grid = GetVariableSDF('0001.sdf', 'Grid.Grid.x');
y_grid = GetVariableSDF('0001.sdf', 'Grid.Grid.y');
z_grid = GetVariableSDF('0001.sdf', 'Grid.Grid.z');

dx = x_grid(2) - x_grid(1);
dy = y_grid(2) - y_grid(1);
dz = z_grid(2) - z_grid(1);

ppc = 50;
ny = 8;
nz = 8;
total_weight = ppc*ny*nz;

%% Convert temperature increase back to energy increase

% Temperature increase in each x bin (max 0.0 prevents rounding issues)
dT = max(temperature - 300.0, 0.0);

% Heat capacity at 300 K
% C = 0.3 + 1.2*T'*(2.2 + T')/(1.1 + T')^2, T' = Z^(-4/3)*T_eV
kb = 1.38064852e-23;
q0 = 1.60217662e-19;
T_eV = 300.0*kb/q0;
Tprime = 73^(-4.0/3.0) * T_eV;
C = 0.3 + 1.2 * Tprime * (2.2 + Tprime)/(1.1 + Tprime)^2;

% Get energy deposition
cell_vol = dx * dy * dz;
ne = 4.015e30;
dE = cell_vol * C * ne * dT * kb;
dE = sum(sum(dE, 3), 2);           % In 3d, sum over y and z directions

%% Create plot

fig = figure;
hold on;
mat_rho = 16.65; %g/cm³
dx_cm = dx*100;
lin = loghist(dE/q0/1e6/dx_cm/mat_rho/total_weight, x_grid/1.94e-4);
lin.LineWidth = 2;

% Plot Lockwood data
plot(Ta_500kev_r, Ta_500kev_E, 'x', 'LineWidth', 2);
ax = gca;
ax.FontSize = 16;
grid on
xlabel('x [CSDA range (195\mum)]');
ylabel('E_{dep} per e^- [MeVcm^{2}/g]');

% Save figures
saveas(fig, 'Lockwood_plot');
saveas(fig, 'Lockwood_plot', 'png');
