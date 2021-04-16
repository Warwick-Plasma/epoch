% plot_hanson.m
% Cycles through the generated SDF Files and pulls out all the probe
% particles. Bin the data into fractional scatter per square degree and
% overlay the data from Hanson (1951) - "Measurement of Multiple Scattering
% of 15.7-Mev Electrons"

% Hanson data (conversion of figure points to data performed on 
% WebPlotDigitizer). Figure 3, 18.66 mg/cm² curve
Hanson = [0.30604444313309287, 0.019275267286819657
0.7950020586700823, 0.018451249194366148
1.2839596742070736, 0.017662457893146567
1.9575136093515875, 0.01456874673954819
2.212842346106025, 0.013375327943278456
2.60762418742883, 0.011890968627170512
3.025581155188852, 0.01045924573851335
4.253728741046124, 0.005882006680565839
4.716426577900207, 0.004459798881416349
5.01703264500883, 0.003682268878578864
6.268623587933377, 0.002092746665421494
6.359446479271802, 0.0017288538815572148
7.1929463391213, 0.0011052829720963738
8.024970930505464, 0.0006288650453052094
9.321570973154138, 0.0003146812882121259
10.363110509678432, 0.00017518939502146638
11.313022462974114, 0.0001107862979055808
12.286914234597187, 0.0000738672025638468
12.332325680266408, 0.00006713862883715454
13.46938214415519, 0.000044745597131401216
15.1173643121829, 0.000025674282458956593
17.04586207312772, 0.000015687056851041443
18.69599008619591, 0.00001066436200632143
21.671526111581258, 0.0000056018473505442885
26.27589263800801, 0.0000023452570878223374
29.649161041646835, 0.0000013828085942321613];

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

%% Extract particles which trigger the probe

% Make get_probe_species ignore its user-input section
skip_user = true;

% Manually set the values of the get_probe_species variables
probe_name = 'electron_probe';
SDF_start = 0;
SDF_end = -1;
SDF_dir = '.';
get_px = true;
get_py = true;
get_pz = true;
get_weight = false; % All weights are set to 1.0 in the input deck
get_time = false;
get_id = false;
get_x = false;
get_y = false;
get_z = false;

% Run get_probe_species
get_probe_species;

%% Process data into fractional scatter per solid angle

fig = figure;
theta = atan2(sqrt(py.^2 + pz.^2), px);
[vals, edges] = histcounts(theta);

% Normalise to fractional scattering per square degree
bin_sq_deg = 2*pi*(cos(edges(1:end-1)) - cos(edges(2:end)))*(180/pi)^2;
vals = (vals/length(px)) ./ bin_sq_deg;

% Convert bin edges to degrees
edges = edges * 180/pi;

% Plot our data
lin = loghist(vals, edges); % loghist has adpative bin sizes
lin.LineWidth = 2;
hold on;
xlabel('\theta [deg]');
ylabel('Frac. scat. per sq. deg. [deg^{-2}]');
ax = gca;
ax.FontSize = 16;
ax.YAxis.Scale = 'Log';
grid on

% Plot Hanson data
plot(Hanson(:,1), Hanson(:,2), ' x', 'LineWidth', 2);

%% Save plot

saveas(fig, 'Hanson_plot.fig');
saveas(fig, 'Hanson_plot', 'png');
