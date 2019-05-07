% Exponential elastic-plastic solution. Comparison between "cohesive"
% numerical model and exponential asymptotic solution

close all
clear all


%% Model set-up

% Non-dimensional parameters
rho_hat = 1e-4;                           % Density 
omega = 1;                             % Rescaled nominal strain-rate
lambda = 1e-3;                            % Plastic work

% Plastic strain-rate parameters
T_p = 1.01;
s_p = 1.01;
beta1 = 1;
beta2 = 1;
gammadot0 = 1;

% Reaction Parameters
T_R = 1.02;
Omega_hat = 1;
A0_hat = 5e-2;
E_hat = 1e3;

% Small parameter 
eps = T_R^2/E_hat;

% Filename
filename = sprintf('SB_rho1em4_omega1_lambda1em3_A5em2_Tp%.0f_Sp%.0f_TR%.0f',T_p*1000,s_p*1000,T_R*1000);
%filename = 'Approach_I';
%% Centreline heat flux

% Heat flux scaling
delta = 0.1;

% Heat flux inhomogeneity (non-dimensional)
h_flux = @(t) 4*t.*exp(-4*t.^2);
dhdt = @(t)  4*(1 - 8*t.^2).*exp(-4*t.^2);

q_flux = @(t) delta*h_flux(t);

%% Elastic-platic parameters

% Critical plastic timescale
t_p = fzero(@(t) ...
    beta1*(T_p - T_elastic_centreline_fun(delta,t,h_flux)) + beta2*(s_p - s_elastic_fun(omega,t)),1e-6);

% EP parameters
a =  delta*integral(@(t) (pi*(t_p-t)).^(-1/2).*dhdt(t),1e-12,t_p);
b = delta*h_flux(t_p);
beta3 = beta1*a + beta2*omega;
rho_hat_0 = rho_hat/eps;

% Combined parameters
Lambda_p = (beta2*(rho_hat_0*beta3)^(1/2))/(beta1*lambda*(1+omega*t_p));
Lambda_R = (Omega_hat*A0_hat)/(b*beta3^(1/2));
Lambda_t = a/beta3;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Numerical solution "cohesive scheme"

% Non-dimensional elastic wave speed
S = sqrt(1/rho_hat);

% Computational grid
N = 201;                       % Number of gridpoints
y = linspace(0,5,N)';          % Spatial domain (half space) 
dy = (y(end)-y(1))/(N-1);      % Grid size
dt = dy/S;                     % Timestep

% Stopping time 
t_stop = 2*t_p;

% Run solver
[T_centreline_max,t_end] =  cohesive_solver_EPreactive(rho_hat,S,lambda,omega,gammadot0,eps,beta1,beta2,T_p,s_p,t_p,q_flux,N,y,dy,dt,filename,t_stop,Omega_hat,A0_hat,E_hat,T_R);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Centreline temperature -- asymptotic 

% Stopping temperature T ~ 1 + O(1) i.e. when correction grows to order 1
T_stop = 2;

% Run solver
[tau,tau_R,T_1,s_1] = sb_centreline_EPreactive(h_flux,beta1,beta2,beta3,rho_hat_0,lambda,omega,gammadot0,a,b,Lambda_p,Lambda_R,Lambda_t,eps,delta,t_p,T_R,T_stop);

% Reaction time 
t_R = t_p + eps*tau_R;

% Save asymptotic solution
save(sprintf('output/asymptotic_%s',filename))
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


