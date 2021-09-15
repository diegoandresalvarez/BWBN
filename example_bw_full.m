%% Bouc-Wen-Baber-Noori (BWBN) hysteresis model
%
% xdd + 2*xi*w0*xd + alpha*w0^2*x(1) + (1-alpha)*w0^2*z = u;
% 
% nueps  = nu0  + deltanu *e;
% Aeps   = A0   - deltaA  *e;
% etaeps = eta0 + deltaeta*e;
%
% zu  = (1/(nueps*(beta+gamma)))^(1/n);
% vs1 = (1 - exp(-p*e))*vs0;
% vs2 = (psi0 + deltapsi*e)*(lambda + vs1);
% h   = 1 - vs1*exp(-((z*sign(xd) - q*zu)^2)/(vs2^2));
%
%
%      h*(Aeps*xd - nueps*(beta*abs(xd)*(abs(z)^(n-1))*z + gamma*xd*(abs(z)^n)))
% zd = -------------------------------------------------------------------------
%                               etaeps;
%
%
%   References:
%
%   - https://en.wikipedia.org/wiki/Bouc-Wen_model_of_hysteresis
%     (Diego wrote most of the WIKIPEDIA page)
%
%   - FOLIENTE, Greg C. "Hysteresis modeling of wood joints and structural
%     systems". Journal of Structural Engineering. Vol. 121. Nro. 6. June.
%     1995.
%
%   - ASTM E2126-11: "Standard test methods for cyclic (reversed) load test 
%     for shear resistance of vertical elements of the lateral force 
%     resisting systems for buildings".
%
% -------------------------------------------------------
% | Developed by:   Diego Andres Alvarez Marin          |
% |                 diegotorquemada@gmail.com           |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% |                                                     |
% |                 Gilberto Alejandro Ortiz Garcia     |
% |                 gialorga@gmail.com                  |
% |                 Universidad Nacional de Colombia    |
% |                 Manizales, Colombia.                |
% -------------------------------------------------------
%
%   Date: 09 - Sep - 2011

%% Beginning:
clear, clc, close all

tic
DEGRADATION = true; % introduce degradation effect
PINCHING    = true; % introduce pinching effect

%% loading experimental data:
m  = 456;                           % mass (kg)
k  = 6.2684;                        % stiffness (kN/mm)
w0 = sqrt(k*(10^6)/m);              % natural frequency (rad/s)

%% External excitation given as an acceleration (u = mm/s^2)
% type_u = 1; % Force (kN)
% type_u = 2; % Acceleration (mm/s^2)
% type_u = 3; % sinusoidal
type_u = 4;   % Patter ASTM-E2126 (Load (kN))
switch type_u
  case 1
    % Here we multiply by 1000 because of the units. The data at the 
    % laboratory was measured in 'kN', and we need acceleration in 'm/s^2'.
    load_kN = dlmread('data/load_kN.txt');
    t       = 0:0.02:(length(load_kN)*0.02)-0.02;    
    u       = 1000*load_kN/m;
  case 2
    % Here we multiply by 9.81 because of the units. The earthquake has
    % units of 'gravity', and we need acceleration in 'm/s^2'.
    u_m = dlmread('data/Northridge_1994_25_sec.txt');
    t   = u_m(:,1);
    u   = 9.81*u_m(:,2);
  case 3
    % Simulated input signal (sinusoidal case)
    t       = 0:0.02:20;
    load_kN = t.*sin(2*pi*t);   % Load in 'kN'
    u       = 1000*load_kN/m;
  case 4
    % Here we multiply by 1000 because of the units. The data returned by
    % function 'generate_signal_astm' is in 'kN', and we need acceleration
    % in 'm/s^2'.
    [t, load_kN] = generate_signal_astm(20, 0.02);
    u            = 1000*load_kN/m;
  otherwise
    error('Invalid external excitation');
end

N       = length(u);                    % number of observations
dt      = t(2)-t(1);                    % Runge-Kutta time step (sec)

%% Setting Bouc-Wen parameters:
xi        =   0.1798;
alpha     =   0.2183;
beta      =   3.0250;
gamma     =  -1.3156;
n         =   1.5173;
paramBW = [ ...
    w0           % natural frequency (rad/s)
    xi           % damping ratio (Adim.)
    alpha        % ratio of post-yield to pre-yield stiffness
    beta         % Bouc-Wen parameter
    gamma        % Bouc-Wen parameter
    n            % Bouc-Wen parameter
];

%% degradation and pinching effect
param_degradation = zeros(0,1);
param_pinching    = zeros(0,1);
if DEGRADATION
    nu0       =   0.4092;
    deltanu   =   3.2397;
    A0        =   2.4929;
    deltaA    =   0.8155;
    eta0      =   3.4873;
    deltaeta  =  -0.8492;
    param_degradation = [ ...   
            nu0          % strength degradation
            deltanu      % strength degradation parameter
            A0           % hysteresis amplitude
            deltaA       % control parameter of 'A' with respect to the energy
            eta0         % stiffness degradation
            deltaeta     % stiffness degradation parameter
    ];    

    % pinching effect
    if PINCHING
        p         =   9.7076;
        vs0       =   1.0962;
        psi0      =   1.2478;
        deltapsi  =  -3.6339;
        lambda    =  -3.2035;
        q         =   1.5662;
        param_pinching = [ ...
            % If you do not want to introduce Pinching effect, comment the
            % following six (6) lines. If you want to introduce Pinching you
            % must introduce also degradation
            p            % parameter which controls initial pinching
            vs0          % pinching severity
            psi0         % pinching parameter
            deltapsi     % parameter which controls change of pinching
            lambda       % pinching parameter
            q            % pinching parameter
        ];        
    end
end

param = [paramBW; param_degradation; param_pinching ];

%% Definition of the state function
BW_real = @(x,u) diff_eq_real(x,u,param);
Fexact  = @(x_k,u_k) rk_discrete(BW_real,x_k,u_k,dt);

%% Initial condition:
% Initial condition x, xd, z, e (Displacement, velocity. hysteretic
% displacement, hysteretic energy) set to zero
x_0 = zeros(4,1);
x_k = zeros(length(x_0)); % vector that contains the system response

x_k(:,1) = x_0;           % Initial state

%% Computing system response
for i = 2:N
    x_k(:,i) = Fexact(x_k(:,i-1),u(i-1));
end
x_k = x_k';
% x_k(:,1) - [mm]   displacement
% x_k(:,2) - [mm/s] velocity
% x_k(:,3) - [mm]   hysteretic displacement 
% x_k(:,4) - [J/kg] dissipated hysteretic energy

%% Computing "Restoring Force" (Fz):
Fz = alpha*k*x_k(:,1) + (1-alpha)*k*x_k(:,3);
toc

%% Plot the results:

%% Displacement:
figure;
plot(t,x_k(:,1),'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Displacement (mm)', 'FontSize', 16);
title('Time vs. Displacement', 'FontSize', 18);
grid on

%% Hysteresis cycles:
figure;
plot(x_k(:,1),Fz,'b');
xlabel('Displacement (mm)', 'FontSize', 16);
ylabel('Restoring force (kN)', 'FontSize', 16);
title('Displacement vs. Restoring force', 'FontSize', 18);
grid on

%% Total dissipated energy
% Compute 'dissipated elastic energy'
%
%                                            /t_f
%  dissipated_elastic_energy = alpha * w0^2 *|     est_displ * est_vel dt
%                                            /t_0
%
% Rememeber that this is a cummulative measure.
%
diss_elastic_energy = 1e-6*alpha*(w0^2)*cumtrapz(t, x_k(:,1).*x_k(:,2));

% tot_diss_energy = diss_elastic_energy + diss_hysteretic_energy
tot_diss_energy = diss_elastic_energy + x_k(:,4);

figure;
plot(t,tot_diss_energy,'b');
xlabel('Time (s)', 'FontSize', 16);
ylabel('Dissipated energy (J/kg)', 'FontSize', 16);
title('Time vs. Dissipated energy', 'FontSize', 18);
grid on

%save histeresis_4 t u Fz x_k tot_diss_energy

%% END
