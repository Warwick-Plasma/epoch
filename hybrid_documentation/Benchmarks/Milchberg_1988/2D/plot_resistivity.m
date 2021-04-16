% check_resistivity.m
% Reads the temperature in EPOCH, and when we set the valence number to 
% match that in the input deck, we compare the output resistivity to that
% expected from our model

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

q0 = 1.60217662e-19;
kb = 1.38064852e-23;

%% Get EPOCH data

Te = GetVariableSDF('0000.sdf', 'Hybrid.Electron_temperature.data');
eta_epoch = GetVariableSDF('0000.sdf', 'Hybrid.Resistivity.data');

Te = Te(:,1)*kb/q0;
eta_epoch = eta_epoch(:,1);

%% Compare to Milchberg

load Milch_data
fig = figure;
plot(temperature_eta(:,1),temperature_eta(:,2),'x ', 'LineWidth', 2);
grid on
hold on
plot(Te, eta_epoch, 'LineWidth', 2);
xlabel('Te [eV]');
ylabel('Resistivity [\Omegam]');
ax = gca;
ax.FontSize = 16;
legend('Milchberg (1988)', 'EPOCH');

saveas(fig, 'reduced_Lee_More');
saveas(fig, 'reduced_Lee_More', 'png');
