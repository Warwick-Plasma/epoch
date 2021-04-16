% plot_data.m
% Creates a plot of the energy spectrum and angular distribution of a
% particle type at the probe boundary

%% User input

is_electron = true;

%% Pre-amble

c = 299792458;
me = 9.109383560000000e-31;
qe = 1.602176565e-19;

if is_electron
    E_file = 'electrons_energy.txt';
    p_file = 'electrons_p.txt';
end

%% Extract data

% Get particle energies
fileID = fopen(E_file);
E = fscanf(fileID,'%g');
fclose(fileID);

% Get electron momenta
fileID = fopen(p_file);
p_temp = fscanf(fileID,'%g %g %g\n');
fclose(fileID);
mask = mod(floor(linspace(1,length(p_temp),length(p_temp))),3);
px = p_temp(mask==1);
py = p_temp(mask==2);
pz = p_temp(mask==0);
clear p_temp

%% Energy distribution

figure;
[vals, edges] = histcounts(E);
centres = 0.5*(edges(2:end)+edges(1:end-1));
vals = vals./(edges(2:end)-edges(1:end-1));
plot(centres, vals, 'LineWidth', 2);
xlabel('E [MeV]');
ylabel('dN/dE [MeV^{-1}]');
ax = gca;
ax.FontSize = 16;
ax.YAxis.Scale = 'Log';
grid on

%% Angular distribution

figure;
theta = atan2(sqrt(py.^2 + pz.^2),px);
[vals, edges] = histcounts(theta);
centres = 0.5*(edges(2:end)+edges(1:end-1));
vals = vals./(edges(2:end)-edges(1:end-1));
plot(centres, vals, 'LineWidth', 2);
xlabel('\theta [rad]');
ylabel('dN/d\theta [rad^{-1}]');
ax = gca;
ax.FontSize = 16;
ax.YAxis.Scale = 'Log';
grid on