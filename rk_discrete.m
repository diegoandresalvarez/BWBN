function x_th = rk_discrete(diff_eq,x_t,u,h)
%%  x_th = rk_discrete(diff_eq,x_t,u,h)
%   Continuous to discrete nonlinear system using fourth order Runge-Kutta
%
%   Input data:
%
%   - diff_eq: Handle function with the differential equation
%   - x_t    : State vector at time 't'
%   - u      : Exogenous input at time 't'
%   - h      : Integration step
%
%   Output data:
%
%   - x_th  : State vector at time 't+h'
%
%   Source:
%
%   - http://mathworld.wolfram.com/Runge-KuttaMethod.html
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
%   Date: 25 - Aug - 2011

%% 4th order Runge-Kutta:
k1 = h*diff_eq(x_t,        u);
k2 = h*diff_eq(x_t + k1/2, u);
k3 = h*diff_eq(x_t + k2/2, u);
k4 = h*diff_eq(x_t + k3,   u);

x_th = x_t + (k1 + 2*k2 + 2*k3 + k4)/6;

%% END
