% plot_matlab_test.m
% Performs post processing on EPOCH data to inspect the energy/particle
% loss after two bunches hit a tnsa boundary.
%
% We should find that half the electrons have escaped, and the remaining
% electrons have half their original energy

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
qe = 1.60217662e-19;
me = 9.10938356e-31;
eps0 = 8.85418782e-12;
h_planck = 6.62607004e-34;
hbar = h_planck/(2.0*pi);
c = 299792458;
amu = 1.66054e-27;
kb = 1.38064852e-23;

% Input deck variables
input_p_loss = 5.341e-22;
input_ang = 20; % [degree]

%% Get field values

% Pull out initial EPOCH data
data = GetDataSDF('0001.sdf');
px0_esc = data.Particles.Px.Electron_esc.data;
px0_ref = data.Particles.Px.Electron_ref.data;

% Pull out data after one reflux
data = GetDataSDF('0020.sdf');
px1_ref = data.Particles.Px.Electron_ref.data;
py1_ref = data.Particles.Py.Electron_ref.data;
pz1_ref = data.Particles.Pz.Electron_ref.data;

KE_bef = (px0_esc(1)/qe*c - me*c^2/qe);
if ~exist('data.Particles.Px.Electron_esc.data', 'var')
    fprintf(['\nParticles of initial kinetic energy %g eV have all' ...
        ' escaped\n\n'], KE_bef);
else
    fprintf('\nFAIL: Particles of initial kinetic energy %g eV remain', ...
        KE_bef);
end

p_bef = px0_ref(1);
p_aft = sqrt(px1_ref(1)^2 + py1_ref(1)^2 + pz1_ref(1)^2);
reduction = abs(p_aft) - abs(p_bef);
if ~exist('data.Particles.Px.Electron_esc.data', 'var')
    fprintf(['Particles of initial p = %g kgm/s end with' ...
        ' p = %g kgm/s\nThis is a %g kgm/s change\nExpected change:' ...
        ' %g kgm/s\n\n'], p_bef, p_aft, reduction, -input_p_loss);
end

%% Have angles been sampled correctly?

% Get theta values
theta = atan2(sign(py1_ref).*sqrt(py1_ref.^2 + pz1_ref.^2), -px1_ref) ...
    *180/pi;

% Plot histogram of scatter angles
fig = figure;
histogram(theta, 'Normalization','countdensity');
grid on
xlabel(['\theta [', char(176), ']']);
ylabel(['dN/d\theta [(', char(176), ')^{-1}]']);
ax = gca;
ax.FontSize = 16;

% Overlay expected distribution
hold on
part_no = length(px1_ref);
plot(0.5*input_ang*[-1,-1,1,1], part_no/input_ang*[0,1,1,0], 'r', ...
    'LineWidth', 2);

%% Save figure

saveas(fig, 'tnsa_angle_test');
saveas(fig, 'tnsa_angle_test', 'png');