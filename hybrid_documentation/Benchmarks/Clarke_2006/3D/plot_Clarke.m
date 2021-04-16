% plot_Clarke.m
% Performs post processing on EPOCH data to reproduce Fig 13 in Clarke
% (2006), showing the bremsstrahlung energy distribution after a laser-shot

% Clarke Fig. 13, complete with errorbars. Energy x-axis has been converted
% to units of MeV

Clarke_E = [79.74506213781405 
            496.52039296190753
            795.047698632799
            996.9867441557798
            4019.703048156543
            13041.700036817969
            31670.528303095427]*1.0e-3;

Clarke_dEdE = [4175318936560.358
               82888694169.76395
               54930292577.17582
               50372598801.096924
               291038097.8287608
               29314645.016655203
               3435976.277539999];
           
Clarke_eldEdE = [71.53550134185765
                 249.52648267318224
                 496.52039296190804
                 600.4901635817216
                 3415.235905750088
                 8144.7566458059155
                 23069.629995124185]*1.0e-3;
Clarke_eldEdE = Clarke_E - Clarke_eldEdE;
             
Clarke_ehdEdE = [1320.0083274979147
                 2981.5561199135745
                 3807.1753802048756
                 4019.7030481565307
                 103687.73600531375
                 18066.779075568367
                 55017.14300400689]*1.0e-3;
Clarke_ehdEdE = Clarke_ehdEdE - Clarke_E;

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

%% Estimate laser power

% Intensity in W/m²
I0 = 4.0e20 * 1.0e4;

% In 3D, inject double Gaussian up to r = HWHM
y_fwhm = 5.0e-6;
y_sig = 0.5*y_fwhm/sqrt(2*log(2));
y_int = pi*y_sig^2;

% Integral of temporal Gaussian between +/- FW0.1M/2
t_fwhm = 800.0e-15;
t_sig = 0.5*t_fwhm/sqrt(2*log(2));
t_fw01m = 2.0*sqrt(2.0*log(10))*t_sig;
t_int = sqrt(0.5*pi)*t_sig*(erf(0.5*t_fw01m/(sqrt(2)*t_sig)) ...
    - erf(-0.5*t_fw01m/(sqrt(2)*t_sig)));

% Laser energy
E_las = I0 * t_int * y_int;

%% Get X-rays

% Ignore the user-input section of get_field
px = GetVariableSDF('0100.sdf','Particles.Px.Brem.data');
py = GetVariableSDF('0100.sdf','Particles.Py.Brem.data');
pz = GetVariableSDF('0100.sdf','Particles.Pz.Brem.data');
weight = GetVariableSDF('0100.sdf','Particles.Weight.Brem.data'); 

%% X-ray energy spectrum

% Get X-ray energies in MeV
E_xray_MeV = sqrt(px.^2 + py.^2 + pz.^2)*c/(1.0e6*q0);

% Only get X-rays in the "40-degree cone"
theta = atan2(sqrt(py.^2 + pz.^2), px)*180/pi;
in_cone = theta < 40;
E_xray_cone = E_xray_MeV(in_cone);
weight_cone = weight(in_cone);

% Bin X-ray energies into logarithmically spaced bins
[~, log_edge, iE] = histcounts(log10(E_xray_cone));
bin_edges = 10.^log_edge;
bin_centres = 0.5*(bin_edges(2:end) + bin_edges(1:end-1));

% Sort photon count into bins
bin_vals = 0*bin_centres;
for i = 1:length(iE)
    bin_vals(iE(i)) = bin_vals(iE(i)) + weight_cone(i);
end

% Convert to dN/dE
bin_vals = bin_vals./(bin_edges(2:end)-bin_edges(1:end-1));

% Convert to dN/dE . 1/E_las
bin_vals = bin_vals/E_las;

%% Create benchmark plot

fig1 = figure;
errorbar(Clarke_E, Clarke_dEdE,0*Clarke_E,0*Clarke_E,Clarke_eldEdE, ...
    Clarke_ehdEdE, 'x ', 'LineWidth', 2);
ax = gca;
ax.XScale = 'log';
ax.YScale = 'log';
ax.FontSize = 16;
xlabel('E\gamma [MeV]');
ylabel('dN/dE\gamma/E_{las} [MeV^{-1}J^{-1}]');
grid on;
hold on;
lin = loghist(bin_vals,bin_edges);
lin.LineWidth = 2;
legend('Clarke (2006)', 'EPOCH-hybrid', 'Location', 'SouthWest');

% Confirm full emission by plotting temporal rate
fig2 = figure;
skip_user = true;
in_subset = false;
subset_name = '';
particle_name = 'Brem';
SDF_start = 0;
SDF_end = -1;
SDF_dir = '.';
plot_energy = false;
mass = 0;
KE_min = 0;
plot_count = false;
removed_partilces = false;
dim_no = 1;
use_parallel = false;
particle_rate;
xlabel('Time [s]');
ylabel('dN/dt [s^{-1}m^{-1}]');

saveas(fig1, 'Clarke_benchmark');
saveas(fig1, 'Clarke_benchmark','png');
saveas(fig2, 'Xray_emission');
saveas(fig2, 'Xray_emission','png');
