% Equilibration_prototype.m
% Obtains the thermal evolution of an electron and ion species over time,
% using the equations found in Spitzer (1962) - "Physics of Fully Ionized
% Gases", Second Edition, Number 3, Interscience Tracts of Physics and
% Astronomy. Available for free on archive.org (last accessed 5th June 
% 2020)

%% User input

Te0 = 100.0;    % Initial electron temperature [eV]
Ti0 = 50.0;     % Initial ion temperature [eV]

targ_Z = 13;
targ_A = 27;
targ_ni = 6.02e28;

t_end = 1000.0e-15;
dt = 1.0e-15;

%% Constants

% Physical
qe = 1.60217662e-19;
me = 9.10938356e-31;
eps0 = 8.85418782e-12;
h_planck = 6.62607004e-34;
hbar = h_planck/(2.0*pi);
c = 299792458;
amu = 1.66054e-27;
kb = 1.38064852e-23;

% Derived
c_eq = 2.0/3.0*(1.0/(2.0*pi))^1.5 * qe^4/kb^1.5 ...
    * dt*sqrt(me*targ_A*amu)/eps0^2;

%% Pre-amble

rho = targ_ni * targ_A * amu/1000; % [g/cm³]
rho1 = rho / (targ_A * targ_Z);

%Convert temperatures to Kelvin (like EPOCH)
time = 0:dt:t_end;
Te = 0*time + Te0*qe/kb;
Ti = 0*time + Ti0*qe/kb;

%% Algorithm

for i = 2:length(time)
    
    % Thomas-Fermi to get Z*
    T1 = kb*Te(i)/qe/targ_Z^(4.0/3.0);
    Tf = T1./(1.0 + T1);
    A = 0.003323*T1.^0.9718 + 9.26148e-5*T1.^3.10165;
    B = -exp(-1.763 + 1.43175*Tf + 0.31546*Tf.^7);
    C = -0.366667*Tf + 0.983333;
    Q = (rho1.^C + A.^C.*rho1.^(B.*C)).^(1./C);
    x = 14.3139*Q.^0.6624;
    Z_st = targ_Z*x./(1 + x + sqrt(1 + 2*x));
    
    % Lee-More methods to get Coulomb logarithm
    fermi_e = hbar^2/(2*me*qe)*(3*pi^2.*Z_st*targ_ni).^(2/3);
    debye_huk = (Z_st * targ_ni * qe / eps0 * ...
        (1 / (sqrt((kb*Te(i)/qe)^2 + fermi_e^2))))^(-0.5);
    R0 = (4*pi*targ_ni/3)^(-1/3);
    b_max = max([debye_huk,R0]);
    ve = sqrt(3*kb*Te(i)/me);
    uncert = h_planck/(2*me*ve);
    impact = Z_st*qe^2/(4*pi*eps0*me*ve^2);
    b_min = max([uncert, impact]);
    cou_log = max([0.5*log(1 +(b_max/b_min)^2), 2]);
        
    % Update temperatures
    diff = Ti(i-1) - Te(i-1);
    Te(i) = Te(i-1) + diff * c_eq * Z_st^2 * targ_ni * cou_log / ...
        (Te(i-1)*amu*targ_A + Ti(i-1)*me)^(1.5);
    Ti(i) = Ti(i-1) - diff * c_eq * Z_st^3 * targ_ni * cou_log / ...
        (Te(i-1)*amu*targ_A + Ti(i-1)*me)^(1.5);
end
