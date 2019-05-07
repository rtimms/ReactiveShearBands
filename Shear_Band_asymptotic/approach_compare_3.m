close all
clear all

%% Numerical 

name = 'numerical_Approach_I_N201';
load output/numerical_Approach_I_N201.mat
load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

load(sprintf('output/asymptotic_%s',FILENAME))

% Increment for numerical plots
[~,idx] = min(abs(t_R - t));
%tt = [1:4:idx idx+2:2:numel(t)];
tt = [1:8:numel(t)];

figure(1);
hold on
plot(t(tt),T(1,tt),'bo')


figure(2);
hold on
plot(t(tt),stress(1,tt),'bo')

%% Aproach I

load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

t_R1 = t_R;




figure(1)
hold on
plot(t_p + eps*(tau'),T_elastic + eps*T_1','k-')

figure(2)
hold on
plot(t_p + eps*(tau'),s_elastic + eps*s_1','k-')

%% Approach II
load output/Approach_II.mat

t_R2 = t_R;

figure(1)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau') + eps*T_1','k-.')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat - tau_hat(1)), T_R + eps*theta_psi(:,1),'k-.')
figure(2)
hold on
plot(t_p + eps*(tau'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau') + eps*s_1','k-.')
plot(t_p + eps*(tau_R) + eps^(3/2)*(tau_hat - tau_hat(1)), s_R + eps*theta_psi(:,2),'k-.')

%% Approach III
load output/Approach_III.mat

figure(1)
hold on
plot(t_p + eps*(tau(t_idx)'),arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau(t_idx)') + eps*T_1(t_idx)','k--')
plot(time_reaction3,T_reaction3,'k--')
figure(2)
hold on
plot(t_p + eps*(tau(t_idx)'),arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau(t_idx)') + eps*s_1(t_idx)','k--')
plot(time_reaction3,(beta1/beta2)*(T_p - T_reaction3) + s_p + (eps/beta2)*log(((g_star/t_star)*time_reaction3/gammadot0)*eps^(1/2)),'k--')

%%
figure(1)
plot(t_p*ones(1,2),linspace(1,1.025,2),'k:',t_R1*ones(1,2),linspace(1,1.025,2),'r:',...
    t_R2*ones(1,2),linspace(1,1.025,2),'r:',t_star*ones(1,2),linspace(1,1.025,2),'r:')
ylim([1 1.025])
xlim([0 0.18])
legend('Numerical','Approach I','Approach II','Approach III','location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_approch3_compare_T.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(2)
plot(t_p*ones(1,2),linspace(0.995,1.01,2),'k:',t_R1*ones(1,2),linspace(0.995,1.01,2),'r:',...
    t_R2*ones(1,2),linspace(0.995,1.01,2),'r:',t_star*ones(1,2),linspace(0.995,1.01,2),'r:')
ylim([0.995 1.01])
xlim([0 0.18])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_approch3_compare_s.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(1)
legend off
xlim([0.14 0.18])
ylim([1.014 1.022])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_approch3_compare_T_zoom.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(2)
xlim([0.14 0.18])
ylim([0.999 1.0055])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_approch3_compare_s_zoom.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
