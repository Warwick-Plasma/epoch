% plot_matlab_test.m
% Performs post processing on EPOCH data to obtain

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

%% Extract particles which trigger the probe

% Make get_probe_species ignore its user-input section
skip_user = true;

% Manually set the values of the get_probe_species variables
probe_name = 'x1_probe';
SDF_start = 0;
SDF_end = -1;
SDF_dir = '.';
get_px = true;
get_py = true;
get_pz = true;
get_weight = true; % All weights are set to 1.0 in the input deck
get_time = true;
get_id = true;
get_x = false;
get_y = false;
get_z = false;

% Run get_probe_species
get_probe_species;

% Save variables
x1probe_t = probe_time;
x1probe_id = id;
x1probe_w = weight;
x1probe_px = px;
x1probe_py = py;
x1probe_pz = pz;

% Repeat for y injector
probe_name = 'x2_probe';
get_probe_species;
x2probe_t = probe_time;
x2probe_id = id;
x2probe_w = weight;
x2probe_px = px;
x2probe_py = py;
x2probe_pz = pz;

%% Confirm injection and probes are accurate

% Particles enter at 1, 2, 3, 4, and 5 fs, and trigger probes after moving 
% 5 um
c = 299792458;
me = 9.10938356e-31;
x1probe_vx = x1probe_px*c^2./ ...
    (sqrt((x1probe_px.^2 + x1probe_py.^2 + x1probe_pz.^2)*c^2 + me^2*c^4));
time_to_probe = 5.0e-6/x1probe_vx(1);

expected_entry_times = [1.0, 2.0, 3.0, 4.0, 5.0]*1.0e-15;
expected_probe_times = expected_entry_times + time_to_probe;

fig1 = figure;
% Plot theory expectation line
plot(expected_entry_times, expected_probe_times, 'LineWidth', 2);
grid on;
ax = gca;
ax.FontSize = 16;
xlabel('Entry time');
ylabel('Probe time');

% Plot EPOCH line
hold on
plot(expected_entry_times, x1probe_t, 'x');
plot(expected_entry_times, x2probe_t, 'o');
legend('Expected', 'xmin-injector', 'xmax-injector', 'Location', 'NorthWest');

saveas(fig1, 'file_injector_probe_time');
saveas(fig1, 'file_injector_probe_time', 'png');

%% Additional ID check

fprintf('\nIn the xmin-probe, we detect ID values: %i, %i, %i, %i, %i\n', ...
    x1probe_id(1), x1probe_id(2), x1probe_id(3), x1probe_id(4), x1probe_id(5));
fprintf('In the xmax-probe, we detect ID values: %i, %i, %i, %i, %i\n\n', ...
    x2probe_id(1), x2probe_id(2), x2probe_id(3), x2probe_id(4), x2probe_id(5));
