% plot_G4_delta.m
% Performs post processing on EPOCH data to reproduce the electron energy
% and angular spectra from a Geant4 run. The G4 simulation had elastic
% scatter and bremsstrahlung deactivated, and simulated 100,000 electrons
% with KE 50 MeV passing through 100 um of Al. The output is stored in
% G4_100um_Al_50MeV.mat

% Load G4 data [MeV/c]
load G4_100um_Al_50MeV.mat
g4_px = px; 
g4_py = py;
g4_pz = pz;

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
mc2_MeV = me*c^2/1.0e6/q0;

%% Get field values

% Pull out EPOCH data
epoch_px = GetVariableSDF('0001.sdf', 'electron_probe.Px.data');
epoch_py = GetVariableSDF('0001.sdf', 'electron_probe.Py.data');
epoch_pz = GetVariableSDF('0001.sdf', 'electron_probe.Pz.data');

%% Bin data

% Get Geant4 kinetic energies [MeV]
g4_KE = sqrt((g4_px.^2 + g4_py.^2 + g4_pz.^2) + mc2_MeV^2) - mc2_MeV;
[g4_vals_KE, g4_edges_KE] = histcounts(g4_KE);
g4_vals_KE = g4_vals_KE./(g4_edges_KE(2:end)-g4_edges_KE(1:end-1));

% Get EPOCH kinetic energies [MeV]
epoch_KE = (sqrt((epoch_px.^2 + epoch_py.^2 + epoch_pz.^2)*c^2 ...
    + me^2*c^4) - me*c^2)/1.0e6/q0;
[epoch_vals_KE, epoch_edges_KE] = histcounts(epoch_KE);
epoch_vals_KE = epoch_vals_KE./ ...
    (epoch_edges_KE(2:end)-epoch_edges_KE(1:end-1));

% Get Geant4 angles [rad]
g4_ang = atan2(sqrt(g4_py.^2 + g4_pz.^2), g4_px);
[g4_vals_ang, g4_edges_ang] = histcounts(g4_ang);
g4_vals_ang = g4_vals_ang./(g4_edges_ang(2:end)-g4_edges_ang(1:end-1));

% Get Epoch angles [rad]
epoch_ang = atan2(sqrt(epoch_py.^2 + epoch_pz.^2), epoch_px);
[epoch_vals_ang, epoch_edges_ang] = histcounts(epoch_ang);
epoch_vals_ang = epoch_vals_ang./ ...
    (epoch_edges_ang(2:end)-epoch_edges_ang(1:end-1));

%% Create energy spectrum plot

% Epoch data
fig1 = figure;
lin1 = loghist(epoch_vals_KE, epoch_edges_KE);
lin1.LineWidth = 2;
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('KE [MeV]');
ylabel('dN/dE [MeV^{-1}]');
ax.XLim = [0, 50.5];

% Overlay G4 data
hold on;
lin2 = loghist(g4_vals_KE, g4_edges_KE);
lin2.LineWidth=2;
legend('EPOCH','Geant4','Location','northwest');

%% Create angular scatter plot

% Epoch data
fig2 = figure;
lin1 = loghist(epoch_vals_ang, epoch_edges_ang);
lin1.LineWidth = 2;
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('\theta [rad]');
ylabel('dN/d\theta [rad^{-1}]');

% Overlay G4 data
hold on;
lin2 = loghist(g4_vals_ang, g4_edges_ang);
lin2.LineWidth=2;
legend('EPOCH','Geant4','Location','northwest');

 %% Save figures

saveas(fig1, 'Energy_spectrum');
saveas(fig1, 'Energy_spectrum', 'png');
saveas(fig2, 'Angular_dist');
saveas(fig2, 'Angular_dist', 'png');
