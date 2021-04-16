% plot_G4_delta.m
% Performs post processing on EPOCH data to reproduce the electron and ion 
% temperature evolutions in a MATLAB prototype. This ensures consistency
% between the two codes

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

%% Get field values

% Pull out EPOCH data
epoch_time = zeros(1,21);
epoch_Te = zeros(1,21);
epoch_Ti = zeros(1,21);
for i = 0:20
    if i < 10
        file_name = ['000', int2str(i),'.sdf'];
    elseif i < 100
        file_name = ['00', int2str(i), '.sdf'];
    elseif i < 1000
        file_name = ['0', int2str(i), '.sdf'];
    else
        file_name = [int2str(i), '.sdf'];
    end
    data=GetDataSDF(file_name);
    
    epoch_time(i+1) = data.time;
    epoch_Te(i+1) = data.Hybrid.Electron_temperature.data(1,1);
    epoch_Ti(i+1) = data.Hybrid.Ion_temperature.data(1,1);
end


%% Create plot

% MATLAB prototype
fig1 = figure;
hold on
xlabel('Time [ps]');
ylabel('Temperature [eV]');
grid on;
ax = gca;
ax.FontSize = 16;

% Epoch data
plot(epoch_time*1.0e12,epoch_Te*kb/qe,'--','LineWidth',6);
plot(epoch_time*1.0e12,epoch_Ti*kb/qe,'--','LineWidth',6);
legend('MATLAB Te', 'MATLAB Ti', 'Location', 'SouthEast');

 %% Save figures

saveas(fig1, 'Temperature_evolution');
saveas(fig1, 'Temperature_evolution', 'png');
