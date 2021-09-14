function [t,x_t] = generate_signal_astm(dur, dt)
%% [t,x_t] = generate_signal_astm(dur, dt)
%
% This function returns the pattern load shown in the standard ASTM-E2126
% (Pattern B).
%
% Input data:
%
% - dur: Duration of the load (s)
% - dt:  Sampling time (s)
%
% Output data:
%
% - t:   Time vector
% - x_t: Signal
%
%   Bibliography:
%
%  - ASTM E2126-11: "Standard test methods for cyclic (reversed) load test for
%    shear resistance of vertical elements of the lateral force resisting
%    systems for buildings".
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
%   Date: 20 - Mar - 2012

%% 

% Generate time vector 't'
t = 0:dt:dur;

% Generate signal
n   = length(t);        % Compute number of entries for 'x_t'
x_t = zeros(n,1);       % Allocate space in memory
A   = 1;                % Initial amplitude signal
j   = 1;                % Auxiliary counter
for i = 2:n
  x_t(i) =  A*sin(2*pi*t(i));
  if x_t(i)*x_t(i-1) < 0
    j = j+1;
    if j == 6
      A = 1.5*A;
      j = 0;
    end
  end
end

%% Plot resulting signal
figure
plot(t,x_t);
xlabel('Time (s)',    'FontSize', 16);
ylabel('Load (kN)',   'FontSize', 16);
title('Pattern load', 'FontSize', 18);
grid on;

%% END