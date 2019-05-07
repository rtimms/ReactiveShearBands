%% Centreline Plots

close all
clear all

%% Load numerical data beta4 = 0

name = 'numerical_SB_harden_beta4_0_decoupled_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

% Plots

% Increment for numerical plots
tt = [1:50:numel(t)];

figure(101);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*(tau'), T_elastic + eps*T_1','k-')

figure(102);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*(tau'), T_elastic + eps*T_1','k-')

figure(103);
hold on
plot(t(tt),T(1,tt),'bx',t_p + eps*(tau'), T_elastic + eps*T_1','k-')

figure(201);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*(tau'),s_elastic + eps*s_1','k-')

figure(202);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*(tau'),s_elastic + eps*s_1','k-')

figure(203);
hold on
plot(t(tt),stress(1,tt),'bx',t_p + eps*(tau'),s_elastic + eps*s_1','k-')

%% Load numerical data beta4 = 1

name = 'numerical_SB_harden_beta4_1_decoupled_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

% Plots

% Increment for numerical plots
tt = [1:50:numel(t)];

figure(101);
hold on
plot(t(tt),T(1,tt),'ro',t_p + eps*(tau'), T_elastic + eps*T_1','k--')


figure(201);
hold on
plot(t(tt),stress(1,tt),'ro',t_p + eps*(tau'),s_elastic + eps*s_1','k--')

%% Load numerical data beta4 = 10

name = 'numerical_SB_harden_beta4_10_decoupled_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

% Plots

% Increment for numerical plots
tt = [1:50:numel(t)];

figure(102);
hold on
plot(t(tt),T(1,tt),'ro',t_p + eps*(tau'), T_elastic + eps*T_1','k--')


figure(202);
hold on
plot(t(tt),stress(1,tt),'ro',t_p + eps*(tau'),s_elastic + eps*s_1','k--')

%% Load numerical data beta4 = 0,1

name = 'numerical_SB_harden_beta4_01_decoupled_N201';

load(sprintf('output/%s.mat',name))
t = fread(fopen(sprintf('output/time_%s.bin',name)),'double');
data = fread(fopen(sprintf('output/out_%s.bin',name)),[5 numel(t)*round(N/10)],'double')';
vel = reshape(data(:,1),length(data(:,1))/length(t),length(t));
stress = reshape(data(:,2),length(data(:,2))/length(t),length(t));
gdot = reshape(data(:,3),length(data(:,3))/length(t),length(t));
g = reshape(data(:,4),length(data(:,4))/length(t),length(t));
T = reshape(data(:,5),length(data(:,5))/length(t),length(t));

% load asymptotic data
load(sprintf('output/asymptotic_%s',FILENAME))

% compute elastic temp and stress
T_elastic = arrayfun(@(t) T_elastic_centreline_fun(delta,t,h_flux),t_p+eps*tau');
s_elastic = arrayfun(@(t) s_elastic_fun(omega,t),t_p+eps*tau');

% Plots

% Increment for numerical plots
tt = [1:50:numel(t)];

figure(103);
hold on
plot(t(tt),T(1,tt),'ro',t_p + eps*(tau'), T_elastic + eps*T_1','k--')


figure(203);
hold on
plot(t(tt),stress(1,tt),'ro',t_p + eps*(tau'),s_elastic + eps*s_1','k--')

%% save plots
figure(101)
xlabel('$t$');
ylabel('$T$');
xlim([0 0.25])
ylim([1 1.04])
legend('Numerical $\beta_4 = 0$','Asymptotic $\beta_4 = 0$','Numerical $\beta_4 = 1$','Asymptotic $\beta_4 = 1$','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_1_T.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
figure(102)
xlabel('$t$');
ylabel('$T$');
xlim([0 0.25])
ylim([1 1.04])
legend('Numerical $\beta_4 = 0$','Asymptotic $\beta_4 = 0$','Numerical $\beta_4 = 1$','Asymptotic $\beta_4 = 1$','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_10_T.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
figure(103)
xlabel('$t$');
ylabel('$T$');
xlim([0 0.25])
ylim([1 1.04])
legend('Numerical $\beta_4 = 0$','Asymptotic $\beta_4 = 0$','Numerical $\beta_4 = 0.1$','Asymptotic $\beta_4 = 0.1$','Location','northwest')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_01_T.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(201)
xlabel('$t$');
ylabel('$s$');
xlim([0 0.25])
ylim([0.98 1.01])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_1_s.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');


figure(202)
xlabel('$t$');
ylabel('$s$');
xlim([0 0.25])
ylim([0.98 1.01])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_10_s.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

figure(203)
xlabel('$t$');
ylabel('$s$');
xlim([0 0.25])
ylim([0.98 1.01])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_01_s.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
