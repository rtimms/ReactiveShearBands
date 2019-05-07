% batch plots for reactive shear band results

close all
clear all

% parameter space
Ap_vec = [100];
AR_vec = [1 10 15];

[Ap_plot,AR_plot] =  meshgrid(Ap_vec,AR_vec);

linestyleS = {':','-.','--','-',':','-.','--','-',':','-.','--','-'};
linestyleC = {'g','g','g','g','b','b','b','b','r','r','r','r'};


% Load and plot results

for i = 1:numel(Ap_plot)
    
    % Set parameters
    Ap = Ap_plot(i)
    AR = AR_plot(i)
    At = 0.01
    
    % Load data
    filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR*100,At*100);
    load(filename)
    
    % Plot temperature perturbation
    figure(10)
    hold on
    plot(eta,f,'linestyle',linestyleS{i},'linewidth',1.5,'color',linestyleC{i})
    
    % Plot stress perturbation
    figure(20)
    hold on
    plot(eta,g,'linestyle',linestyleS{i},'linewidth',1.5,'color',linestyleC{i})
    
end

% Plot D&O threshold Ap = 10

% Set parameters
Ap = 10;
AR = 0;
At = 0.1;

% Load data
filename = sprintf('output/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR*100,At*100);
load(filename)

% Plot temperature perturbation
figure(1)
hold on
plot(eta,f,'k--','linewidth',2)

% Plot stress perturbation
figure(2)
hold on
plot(eta,g,'k--','linewidth',2)



% Plot exact result for Ap = 1 AR = At = 0
timep = linspace(eta0,eta_Final,100);
figure(1)
plot(timep,exp(timep),'k-','linewidth',2)
xlabel('$\eta$')
ylabel('$f$')
xlim([-5 15])
ylim([0 15])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure

matlab2tikz('f.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');

figure(2)
plot(timep,exp(timep),'k-','linewidth',2)
xlabel('$\eta$')
ylabel('$g$')
xlim([-5 15])
ylim([0 15])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure

matlab2tikz('g.tikz', 'showInfo', false, ...
     'parseStrings',false,'standalone', false, ...
     'height', '\figureheight', 'width','\figurewidth');

