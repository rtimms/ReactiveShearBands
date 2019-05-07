% Testing of constitutive laws
close all
% Power-law parameters
m = 0.08;
l = -1.28;

% Isothermal
T0 = 300;
Tref = [300 350 400 450 500]/T0;
logg = linspace(1,3,50)';
s0 = 4*1e7;
for i = 1:length(Tref)
%     % Data from power law model (stress,strain-rate)
%     powerlaw_data = [s (1/m)*log(s)-(l/m)*log(Tref(i))];
%     % Least-squares fit
%     fitmat = [powerlaw_data(:,1) ones(size(powerlaw_data(:,1)))];
%     exp_param = (fitmat'*fitmat)\(fitmat'*powerlaw_data(:,2));
%     beta2 = eps*exp_param(1);

    % Data from power law model (stress,strain-rate)
    powerlaw_data = [logg Tref(i)^l*exp(logg).^m];
    % Least-squares fit
    fitmat = [powerlaw_data(:,1) ones(size(powerlaw_data(:,1)))];
    exp_param = (fitmat'*fitmat)\(fitmat'*powerlaw_data(:,2));
    beta2 = eps*(1/exp_param(1));


        
    figure(200)
    hold on
    plot(powerlaw_data(:,1),powerlaw_data(:,2)*s0*1e-6,'kx',logg,(exp_param(1)*logg + exp_param(2))*s0*1e-6,'k-')
    pause
    beta2vec(:,i) = [beta2 powerlaw_data(1,2)];
end
xlabel('log(strain rate) /\si{s^{-1}}')
ylabel('stress /\si{MPa}')
xlim([-1 8])
ylim([20 80])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('isotherm.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

% Fixed strain-rate
gammaref = [0.1 1 10 100 1000];
T = linspace(1,5/3,50)';
s0 = 4*1e7;
for i = 1:length(gammaref)
    % Data from power law model (stress,strain-rate)
    powerlaw_data = [T T.^l*gammaref(i)^m];
    % Least-squares fit
    fitmat = [powerlaw_data(:,1) ones(size(powerlaw_data(:,1)))];
    exp_param = (fitmat'*fitmat)\(fitmat'*powerlaw_data(:,2));
    beta1 = -beta2*exp_param(1);
    
    figure(300)
    hold on
    plot(powerlaw_data(:,1)*T0,powerlaw_data(:,2)*s0*1e-6,'kx',T*T0,(exp_param(1)*T+exp_param(2))*s0*1e-6,'k-')
    pause
    beta1vec(:,i) = [beta1 powerlaw_data(1,2)*s0*1e-6];
end
xlabel('Temperature /K')
ylabel('stress /MPa')
xlim([250 500])
ylim([20 80])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('constant_strainrate.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
beta1vec
beta2vec

%%
% Comparison of strain rate

% Power law stress (no strain-hardening)
m = 0.08;
l = -1.28;

gdot_pow = @(s,T) s.^(1/m).*T.^(-l/m);

% Exponential plastic strain rate (no strain-hardening)
eps = 1e-2;
Tp = 1.01;
sp = 1.01;
gdot0 = eps^(1/2)*sp^(1/m)*Tp^(-l/m);
beta1 = 0.09;
beta2 = 0.198;

gdot_exp = @(s,T) gdot0*eps^(-1/2)*exp(-(1/eps)*(beta1*(Tp - T) + beta2*(sp-s)));

% Isothermal
s = linspace(0.1,1.1,50)';
T = 500/T0;

figure(1)
hold on
plot(gdot_exp(s,T)*1e10,s*s0*1e-6,'ro')
plot(gdot_pow(s,T)*1e10,s*s0*1e-6,'bx')
xlabel('strain rate /s')
ylabel('stress /MPa')
xlim([0 2]*1e10)
ylim([0 max(s)]*s0*1e-6)

% Difference surface
s = linspace(0.5,1.5,100)';
T = linspace(1,500/T0,100);
[smat,Tmat] = meshgrid(s,T);
figure(2)
levels = [-3 -2 -1 0 1 2];
contour(smat*s0*1e-6,Tmat*T0,log10(abs(gdot_exp(smat,Tmat) - gdot_pow(smat,Tmat))),levels)
colormap gray
xlabel('stress /\si{MPa}')
ylabel('temperature /\si{K}')
zlim([0 1])
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('diff_nohardening.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
contour(smat*s0*1e-6,Tmat*T0,log10(abs(gdot_exp(smat,Tmat) - gdot_pow(smat,Tmat))),levels,'Showtext','on')

%%

%%%%%%%%%%%%%%%%%%%%%%%  Strain-hardening  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m = 0.08;
l = -1.28;
n = 0.32;
eps = 1e-2;
Tp = 1.01;
sp = 1.01;
gdot0 = eps^(1/2)*sp^(1/m)*Tp^(-l/m);
beta1 = 0.19;
beta2 = 0.19;


% Isothermal (range) & constant strain rate
T0 = 300;
s0 = 4*1e7;

Tref = [300 350 400 450 500]/T0;
gdot = 1;

g = linspace(0.05,0.3,50)';
for i = 1:length(Tref)
    % Data from power law model (stress,strain-rate)
    powerlaw_data = [g gdot^m*Tref(i)^l*g.^n];
    % Least-squares fit
    fitmat = [powerlaw_data(:,1) ones(size(powerlaw_data(:,1)))];
    exp_param = (fitmat'*fitmat)\(fitmat'*powerlaw_data(:,2));
    beta3 = -beta2*exp_param(1);
    
    figure(400)
    hold on
    plot(powerlaw_data(:,1),powerlaw_data(:,2)*s0*1e-6,'kx',g,(exp_param(1)*g+exp_param(2))*s0*1e-6,'k-')
    pause
    beta3vec1(:,i) = [beta3 powerlaw_data(end,2)*s0*1e-6];
end
xlabel('strain')
ylabel('stress /MPa')
xlim([min(g) max(g)*1.2])
ylim([5 30])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('isotherm_hardening', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

% Isothermal & constant strain rate (range)
T0 = 300;
s0 = 4*1e7;

Tref = 500/T0;
gdot = [0.1 1 10 100 1000];

g = linspace(0.05,0.3,50)';
for i = 1:length(gdot)
    % Data from power law model (stress,strain-rate)
    powerlaw_data = [g gdot(i)^m*Tref^l*g.^n];
    % Least-squares fit
    fitmat = [powerlaw_data(:,1) ones(size(powerlaw_data(:,1)))];
    exp_param = (fitmat'*fitmat)\(fitmat'*powerlaw_data(:,2));
    beta3 = -beta2*exp_param(1);
    
    figure(500)
    hold on
    plot(powerlaw_data(:,1),powerlaw_data(:,2)*s0*1e-6,'kx',g,(exp_param(1)*g+exp_param(2))*s0*1e-6,'k-')
    pause
    beta3vec2(:,i) = [beta3 powerlaw_data(end,2)*s0*1e-6];
end
xlabel('strain')
ylabel('stress /MPa')
xlim([min(g) max(g)*1.2])
ylim([5 30])
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('const_strainrate_hardening', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');

beta3vec1
beta3vec2

%%
% Comparison of strain rate
% Power law stress (no strain-hardening)
m = 0.08;
n = 0.32;
l = -1.28;

s_pow = @(T,g,gdot) gdot.^m.*g.^n.*T.^l;

% Exponential plastic strain rate (no strain-hardening)
eps = 1e-2;
Tp = 1.01;
sp = 1.01;
gp = 1.01;
gdot0 = eps^(1/2)*sp^(1/m)*Tp^(-l/m)*gp^(-n/m);
beta1 = 0.09;
beta2 = 0.198;
beta3 = -0.2;
s_exp = @(T,g,gdot) -(beta3/beta2)*g - (beta1/beta2)*T +...
    (1/beta2)*(beta1*Tp + beta2*sp + beta3*gp) +...
    (eps/beta2)*(log(gdot) - log(gdot0*eps^(-1/2)));


% Fixed strain rate (1) and isotherm range
gdot = 1;
Tref = [1 400/T0 500/T0 600/T0 700/T0];
g = linspace(0.05,0.3,50)';
for i = 1:length(Tref)
    figure(600)
    hold on
    plot(g,s_pow(Tref(i),g,gdot),'k-.',g,s_exp(Tref(i),g,gdot),'k-')
end

    
    
    
    
    
    
