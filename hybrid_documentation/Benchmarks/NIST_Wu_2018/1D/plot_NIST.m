% plot_NIST.m
% Performs post processing on EPOCH data to reproduce Fig 4 (a) in Wu
% (2018). This benchmarks the absoute stopping powers of collisional and 
% radiative energy loss for e- between 10 keV to 1 GeV in Al targets

load('Nist_dat');

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

%% Get ID arrays

% Ignore the user-input section of id_track
skip_user = true;

% Manually set id_track variables
in_subset = false;
particle_name = 'Electron';
SDF_start = 1;
SDF_end = -1;
SDF_dir = '.';
get_x = true;
get_y = false;
get_z = false;
get_px = true;
get_py = true;
get_pz = true;
get_weight = false;
get_time = false;
use_parallel = false;
replace_NaN_last_val = false;

% Run id_track. The id_var(i,j) returns the value of variable var for
% particle i at timestep j
id_track; 

%% Post process to bin stopping powers

id_KE = sqrt((id_px.^2 + id_py.^2 + id_pz.^2)*c^2 + me^2*c^4) - me*c^2;
clear id_px id_py id_pz

id_dE = id_KE(:,2:end)-id_KE(:,1:end-1);
id_dx = id_x(:,2:end)-id_x(:,1:end-1);
clear id_x id_y id_z
       
stop_pow = id_dE(~isnan(id_dE))./id_dx(~isnan(id_dE))/1.0e6/q0; % [MeV/m]

id_KE = id_KE(:,1:end-1); % Energies at the start of each loss in id_dE
stop_KE = id_KE(~isnan(id_dE))/1.0e6/q0; % [MeV]

%% Bin histogram

NIST_edges = [NIST_KE(1); 0.5*(NIST_KE(1:end-1)+NIST_KE(2:end)); ...
    NIST_KE(end)];
[~, ~, ibin] = histcounts(stop_KE, NIST_edges);

mean_stop_pow = zeros(1,max(ibin));
for i=1:max(ibin)
    mean_stop_pow(i) = mean(stop_pow(ibin(:)==i))*1.0e-6; % [MeV/um]
end

%% Create plot

% Plot MATLAB data
fig = figure;
lin = loghist(-mean_stop_pow,NIST_edges);
lin.LineWidth = 2;
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('KE_e [MeV]');
ylabel('Stopping power [MeV/\mum]');

% Overlay NIST data
hold on;
% NIST quotes stopping power in MeVcm²/g, Al density is 2.7 g/cm³
plot(NIST_KE, NIST_tot*2.7 / (1.0e4), 'x', 'LineWidth', 2);
legend('EPOCH', 'NIST', 'Location', 'northwest');
 
saveas(fig, 'EPOCH_vs_NIST');
saveas(fig, 'EPOCH_vs_NIST', 'png');

