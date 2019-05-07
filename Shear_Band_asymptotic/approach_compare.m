% Solution of centreline temperature for reactive shear band using an
% elastic-plastic constitutive law

close all
clear all

% Approach I
name = 'numerical_SB_rho1em4_omega1em1_lambda1em3_A5em2_Tp1010_Sp1010_TR1020_N201';
load(sprintf('output/%s.mat',name))
load(sprintf('output/asymptotic_%s',FILENAME))
% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau')

t_R1 = t_R;

figure(1)
hold on
plot(t_p + eps*(tau'),T_elastic + eps*T_1','k-')

figure(2)
hold on
plot(t_p + eps*(tau'),s_elastic + eps*s_1','k-')

%% Approach II Plastic

load Approach_II.mat

t_R2 = t_R;

figure(1)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','k-.')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat2 - tau_hat2(1)), T_R + eps*theta_psi2(:,1),'k-.')
figure(2)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau') + eps*s_1','k-.')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat2 - tau_hat2(1)), s_R + eps*theta_psi2(:,2),'k-.')
figure(3)
hold on
plot(tau_hat2,theta_psi2(:,1),'k-.')
figure(4)
hold on
plot(tau_hat2,theta_psi2(:,2),'k-.')
%% Approach III Reaction stage

% Solve ODEs for inner temperature and stress
options=odeset('Stats','on','RelTol',1e-6);
[tau_hat3,theta_psi3] = ode15s(@ODE_reactive_sb_inner,[-(t_p+eps*tau_R)/(eps^(1/2)) 1/(10*eps)],[0 0],options,...
    beta1,beta2,beta3,gammadot0,0,A0_hat,Omega_hat,tau_R,s_R,T_1(end),s_1(end));

figure(1)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','k:')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat3 - tau_hat3(1)), T_R + eps*theta_psi3(:,1),'k:')
figure(2)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau') + eps*s_1','k:')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat3 - tau_hat3(1)), s_R + eps*theta_psi3(:,2),'k:')
figure(3)
hold on
plot(tau_hat3,theta_psi3(:,1),'k:')
figure(4)
hold on
plot(tau_hat3,theta_psi3(:,2),'k:')





