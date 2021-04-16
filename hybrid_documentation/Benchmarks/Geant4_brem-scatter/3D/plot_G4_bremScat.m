% plot_G4_bremScat.m
% Performs post processing on EPOCH data to reproduce the electron angular 
% distribution from a Geant4 run. The G4 simulation simulated 100,000 
% electrons with total energy 100 MeV passing through 5 mm of Au. 
%
% The Geant4 simulation has ionisation loss and elastic scatter
% deactivated, and the px, py and pz values of electrons passing a probe at 
% 5 mm have been recorded in the files '../electrons_p_100MeV_5mm_Au.txt'.
%
% All photon processes were deactiavted in the Geant4 physics library (no
% photoelectric affect attenuation). The momentum values of photons after
% were printed after their first step, and are recorded in
% '../photons_p_100MeV_5mm_Au.txt'

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

%% Get particle data

% Pull out EPOCH data
epoch_e_px = GetVariableSDF('0001.sdf', 'electron_probe.Px.data');
epoch_e_py = GetVariableSDF('0001.sdf', 'electron_probe.Py.data');
epoch_e_pz = GetVariableSDF('0001.sdf', 'electron_probe.Pz.data');
epoch_e_w = GetVariableSDF('0001.sdf', 'electron_probe.weight.data');

epoch_p_px = GetVariableSDF('0001.sdf', 'Particles.Px.Photon.data');
epoch_p_py = GetVariableSDF('0001.sdf', 'Particles.Py.Photon.data');
epoch_p_pz = GetVariableSDF('0001.sdf', 'Particles.Pz.Photon.data');
epoch_p_w = GetVariableSDF('0001.sdf', 'Particles.Weight.Photon.data');

%% Get Geant4 data

% Tsai data (documented)

file_e = fopen('../Geant4_data/electrons_p_tsai.txt');
tsai_e_p = fscanf(file_e, '%g %g %g', [3, inf])';
fclose(file_e);

file_p = fopen('../Geant4_data/photons_p_tsai.txt');
tsai_p_p = fscanf(file_p, '%g %g %g', [3, inf])';
fclose(file_p);

% Grichine data (undocumented default)

file_e = fopen('../Geant4_data/electrons_p_grichine.txt');
grich_e_p = fscanf(file_e, '%g %g %g', [3, inf])';
fclose(file_e);

file_p = fopen('../Geant4_data/photons_p_grichine.txt');
grich_p_p = fscanf(file_p, '%g %g %g', [3, inf])';
fclose(file_p);

%% Bin data

% Get Geant4 kinetic energies [MeV]
tsai_e_KE = sqrt((tsai_e_p(:,1).^2 + tsai_e_p(:,2).^2 ...
    + tsai_e_p(:,3).^2) + mc2_MeV^2) - mc2_MeV;
tsai_p_KE = sqrt(tsai_p_p(:,1).^2 + tsai_p_p(:,2).^2 + tsai_p_p(:,3).^2);
grich_e_KE = sqrt((grich_e_p(:,1).^2 + grich_e_p(:,2).^2 ...
    + grich_e_p(:,3).^2) + mc2_MeV^2) - mc2_MeV;
grich_p_KE = sqrt(grich_p_p(:,1).^2 + grich_p_p(:,2).^2 ...
    + grich_p_p(:,3).^2);

% Remove G4 photons below 500 keV (to match EPOCH photons)
tsai_p_p = tsai_p_p(tsai_p_KE > 0.5,:);
tsai_p_KE = tsai_p_KE(tsai_p_KE>0.5);
grich_p_p = grich_p_p(grich_p_KE > 0.5,:);
grich_p_KE = grich_p_KE(grich_p_KE>0.5);

% Bin Geant4 KE (Tsai)
[tsai_vals_e_KE, tsai_edges_e_KE] = histcounts(tsai_e_KE);
[tsai_vals_p_KE, tsai_edges_p_KE] = histcounts(tsai_p_KE);
tsai_vals_e_KE = tsai_vals_e_KE ./ ...
    (tsai_edges_e_KE(2:end)-tsai_edges_e_KE(1:end-1));
tsai_vals_p_KE = tsai_vals_p_KE ./ ...
    (tsai_edges_p_KE(2:end)-tsai_edges_p_KE(1:end-1));

% Bin Geant4 KE (Grichine)
[grich_vals_e_KE, grich_edges_e_KE] = histcounts(grich_e_KE);
[grich_vals_p_KE, grich_edges_p_KE] = histcounts(grich_p_KE);
grich_vals_e_KE = grich_vals_e_KE ./ ...
    (grich_edges_e_KE(2:end)-grich_edges_e_KE(1:end-1));
grich_vals_p_KE = grich_vals_p_KE ./ ...
    (grich_edges_p_KE(2:end)-grich_edges_p_KE(1:end-1));

% Get EPOCH kinetic energies [MeV]
epoch_e_KE = (sqrt((epoch_e_px.^2 + epoch_e_py.^2 + epoch_e_pz.^2)*c^2 ...
    + me^2*c^4) - me*c^2)/1.0e6/q0;
epoch_p_KE = sqrt((epoch_p_px.^2 + epoch_p_py.^2 + epoch_p_pz.^2)*c^2)...
    /1.0e6/q0;
[~, epoch_edges_e_KE, iKE] = histcounts(epoch_e_KE);
epoch_vals_e_KE = zeros(length(epoch_edges_e_KE)-1,1);
for i = 1:length(iKE)
    epoch_vals_e_KE(iKE(i)) = epoch_vals_e_KE(iKE(i)) + epoch_e_w(i); 
end
[~, epoch_edges_p_KE, iKE] = histcounts(epoch_p_KE);
epoch_vals_p_KE = zeros(length(epoch_edges_p_KE)-1,1);
for i = 1:length(iKE)
    epoch_vals_p_KE(iKE(i)) = epoch_vals_p_KE(iKE(i)) + epoch_p_w(i); 
end
epoch_vals_e_KE = epoch_vals_e_KE./ ...
    (epoch_edges_e_KE(2:end)-epoch_edges_e_KE(1:end-1))';
epoch_vals_p_KE = epoch_vals_p_KE./ ...
    (epoch_edges_p_KE(2:end)-epoch_edges_p_KE(1:end-1))';

% Get Geant4 Tsai angles [rad]
tsai_e_ang = ...
    atan2(sqrt(tsai_e_p(:,2).^2 + tsai_e_p(:,3).^2), tsai_e_p(:,1));
tsai_p_ang = ...
    atan2(sqrt(tsai_p_p(:,2).^2 + tsai_p_p(:,3).^2), tsai_p_p(:,1));
[tsai_vals_e_ang, tsai_edges_e_ang] = histcounts(tsai_e_ang);
[tsai_vals_p_ang, tsai_edges_p_ang] = histcounts(tsai_p_ang);
tsai_vals_e_ang = tsai_vals_e_ang./ ...
    (tsai_edges_e_ang(2:end)-tsai_edges_e_ang(1:end-1));
tsai_vals_p_ang = tsai_vals_p_ang./...
    (tsai_edges_p_ang(2:end)-tsai_edges_p_ang(1:end-1));

% Get Geant4 Grich angles [rad]
grich_e_ang = ...
    atan2(sqrt(grich_e_p(:,2).^2 + grich_e_p(:,3).^2), grich_e_p(:,1));
grich_p_ang = ...
    atan2(sqrt(grich_p_p(:,2).^2 + grich_p_p(:,3).^2), grich_p_p(:,1));
[grich_vals_e_ang, grich_edges_e_ang] = histcounts(grich_e_ang);
[grich_vals_p_ang, grich_edges_p_ang] = histcounts(grich_p_ang);
grich_vals_e_ang = grich_vals_e_ang./ ...
    (grich_edges_e_ang(2:end)-grich_edges_e_ang(1:end-1));
grich_vals_p_ang = grich_vals_p_ang./...
    (grich_edges_p_ang(2:end)-grich_edges_p_ang(1:end-1));

% Get Epoch angles [rad]
epoch_e_ang = atan2(sqrt(epoch_e_py.^2 + epoch_e_pz.^2), epoch_e_px);
epoch_p_ang = atan2(sqrt(epoch_p_py.^2 + epoch_p_pz.^2), epoch_p_px);
[~, epoch_edges_e_ang, iang] = histcounts(epoch_e_ang);
epoch_vals_e_ang = zeros(length(epoch_edges_e_ang)-1,1);
for i = 1:length(iang)
    epoch_vals_e_ang(iang(i)) = epoch_vals_e_ang(iang(i)) + epoch_e_w(i); 
end
[~, epoch_edges_p_ang, iang] = histcounts(epoch_p_ang);
epoch_vals_p_ang = zeros(length(epoch_edges_p_ang)-1,1);
for i = 1:length(iang)
    epoch_vals_p_ang(iang(i)) = epoch_vals_p_ang(iang(i)) + epoch_p_w(i); 
end
epoch_vals_e_ang = epoch_vals_e_ang./ ...
    (epoch_edges_e_ang(2:end)-epoch_edges_e_ang(1:end-1))';
epoch_vals_p_ang = epoch_vals_p_ang./ ...
    (epoch_edges_p_ang(2:end)-epoch_edges_p_ang(1:end-1))';

%% Create energy spectrum plot

% Epoch data (electrons)
fig1 = figure;
lin1 = loghist(epoch_vals_e_KE, epoch_edges_e_KE');
lin1.LineWidth = 4;
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('KE [MeV]');
ylabel('dN/dE [MeV^{-1}]');
title('Electrons, 5 mm Au');
% Overlay G4 data
hold on;
lin2 = loghist(tsai_vals_e_KE, tsai_edges_e_KE');
lin2.LineWidth=2;
lin3 = loghist(grich_vals_e_KE, grich_edges_e_KE');
lin3.LineWidth=2;
legend('EPOCH','G4 Tsai','G4 Grichine','Location','northeast');

% Epoch data (photons)
fig2 = figure;
lin1 = loghist(epoch_vals_p_KE, epoch_edges_p_KE');
lin1.LineWidth = 4;
ax = gca;
ax.YScale = 'log';
grid on;
ax.FontSize = 16;
xlabel('KE [MeV]');
ylabel('dN/dE [MeV^{-1}]');
title('Photons, 5 mm Au');
% Overlay G4 data
hold on;
lin2 = loghist(tsai_vals_p_KE, tsai_edges_p_KE');
lin2.LineWidth=2;
lin3 = loghist(grich_vals_p_KE, grich_edges_p_KE');
lin3.LineWidth=2;
legend('EPOCH','G4 Tsai','G4 Grichine','Location','northeast');

%% Create angular scatter plot
 
 % Epoch data (e-)
 fig3 = figure;
 lin1 = loghist(epoch_vals_e_ang, epoch_edges_e_ang');
 lin1.LineWidth = 4;
 ax = gca;
 ax.YScale = 'log';
 grid on;
 ax.FontSize = 16;
 xlabel('\theta [rad]');
 ylabel('dN/d\theta [rad^{-1}]');
 title('Electrons, 5 mm Au');
 % Overlay G4 data
 hold on;
 lin2 = loghist(tsai_vals_e_ang, tsai_edges_e_ang');
 lin2.LineWidth=2;
 lin3 = loghist(grich_vals_e_ang, grich_edges_e_ang');
 lin3.LineWidth=2;
 legend('EPOCH','G4 Tsai','G4 Grich','Location','northeast');
 
 % Epoch data (e-)
 fig4 = figure;
 lin1 = loghist(epoch_vals_p_ang, epoch_edges_p_ang');
 lin1.LineWidth = 4;
 ax = gca;
 ax.YScale = 'log';
 grid on;
 ax.FontSize = 16;
 xlabel('\theta [rad]');
 ylabel('dN/d\theta [rad^{-1}]');
 title('Photons, 5 mm Au');
 % Overlay G4 data
 hold on;
 lin2 = loghist(tsai_vals_p_ang, tsai_edges_p_ang');
 lin2.LineWidth=2;
 lin3 = loghist(grich_vals_p_ang, grich_edges_p_ang');
 lin3.LineWidth=2;
 legend('EPOCH','G4 Tsai','G4 Grichine','Location','northeast');
 
  %% Save figures
 
 saveas(fig1, 'Energy_spectrum_e');
 saveas(fig1, 'Energy_spectrum_e', 'png');
 saveas(fig2, 'Energy_spectrum_p');
 saveas(fig2, 'Energy_spectrum_p', 'png');
 saveas(fig3, 'Angular_dist_e');
 saveas(fig3, 'Angular_dist_e', 'png');
 saveas(fig4, 'Angular_dist_p');
 saveas(fig4, 'Angular_dist_p', 'png');
