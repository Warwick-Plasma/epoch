% plot_g4_rester.m
% Performs post processing on Geant4 data to reproduce Fig 17 in Rester
% (1970). This benchmarks the bremsstrahlung photon production and angular
% distributions of the radiating electrons.

% Data from paper (read using WebPlotDigitizer)
load('rester_data.mat');

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

%% Get photon data

G4_file = '/home/sjm630/bremsstrahlung-G4/Rester_build/photons_p.txt';
fileID = fopen(G4_file);
data = fscanf(fileID, '%g %g %g', [3, Inf]);
fclose(fileID);

% Particle data
px = data(1,:); % [MeV/c]
py = data(2,:);
pz = data(3,:);
weight = 0*px + 1;

% Get photon energies and angles
c = 299792458;
energies = sqrt(px.^2 + py.^2 + pz.^2); % [MeV]
theta = atan2(sqrt(py.^2 + pz.^2), px);

% Total weight of present electrons
el_no = 1.0e5;
theta_width = 5/180*pi;
  
%% Pull out the energy distributions for a given set of angles

% Around 10 degrees
theta_10 = 10/180*pi;
thet_low = theta_10-theta_width;
thet_high = theta_10+theta_width;
solid_10 = 2*pi*(cos(thet_low) - cos(thet_high));
energies_10 = energies(theta > thet_low & theta < thet_high);

% Around 30 degrees
theta_30 = 30/180*pi;
thet_low = theta_30-theta_width;
thet_high = theta_30+theta_width;
solid_30 = 2*pi*(cos(thet_low) - cos(thet_high));
energies_30 = energies(theta > thet_low & theta < thet_high);

% Around 60 degrees
theta_60 = 60/180*pi;
thet_low = theta_60-theta_width;
thet_high = theta_60+theta_width;
solid_60 = 2*pi*(cos(thet_low) - cos(thet_high));
energies_60 = energies(theta > thet_low & theta < thet_high);

%% Bin energies

% Bin photons according to their energies
[vals_10, edges_10, ibin] = histcounts(energies_10);
vals_10 = 0*vals_10;
% Find photon energy contribution in each bin
for i = 1:length(ibin)
    vals_10(ibin(i)) = vals_10(ibin(i)) + weight(i)*energies_10(i);
end

[vals_30, edges_30, ibin] = histcounts(energies_30);
vals_30 = 0*vals_30;
for i = 1:length(ibin)
    vals_30(ibin(i)) = vals_30(ibin(i)) + weight(i)*energies_30(i);
end

[vals_60, edges_60, ibin] = histcounts(energies_60);
vals_60 = 0*vals_60;
for i = 1:length(ibin)
    vals_60(ibin(i)) = vals_60(ibin(i)) + weight(i)*energies_60(i);
end

% Normalise to kdn/dk/dOmega
intens_10 = vals_10./(edges_10(2:end)-edges_10(1:end-1))/solid_10/el_no;
intens_30 = vals_30./(edges_30(2:end)-edges_30(1:end-1))/solid_30/el_no;
intens_60 = vals_60./(edges_60(2:end)-edges_60(1:end-1))/solid_60/el_no;

%% Create plot

fig = figure;
hold on
q0 = 1.60217662e-19;
lin_10 = loghist(intens_10,edges_10);
lin_30 = loghist(intens_30,edges_30);
lin_60 = loghist(intens_60,edges_60);

lin_10.LineWidth = 1;
lin_30.LineWidth = 1;
lin_60.LineWidth = 1;

def_b = [0 0.4470 0.7410];      % Default blue
def_o = [0.8500 0.3250 0.0980]; % Default orange
def_g = [0.9290 0.6940 0.1250]; % Default gold

lin_10.Color = def_b;
lin_30.Color = def_o;
lin_60.Color = def_g;

plot(rester_10(:,1),rester_10(:,2), ' x', 'LineWidth', 1, 'Color', def_b);
plot(rester_30(:,1),rester_30(:,2), ' x', 'LineWidth', 1, 'Color', def_o);
plot(rester_60(:,1),rester_60(:,2), ' x', 'LineWidth', 1, 'Color', def_g);

ax = gca;
ax.YScale = 'log';
ax.FontSize = 16;
grid on;
xlabel('Photon Energy [MeV]');
ylabel('kdn/dkd\Omega per e^- [sr^{-1}]');
