close all
clear all


load SB_harden_fgk_Ap100_Ak0_At.mat

figure(1)
hold on
plot(eta,f,'k:')
figure(2)
hold on
plot(eta,g,'k:')
figure(3)
hold on
plot(eta,k,'k:')

load SB_harden_fgk_Ap100_Ak10_At.mat

figure(1)
plot(eta,f,'k-.')
figure(2)
plot(eta,g,'k-.')
figure(3)
plot(eta,k,'k-.')

load SB_harden_fgk_Ap100_Ak100_At.mat

figure(1)
plot(eta,f,'k--')
figure(2)
plot(eta,g,'k--')
figure(3)
plot(eta,k,'k--')

load SB_harden_fgk_Ap100_Ak1000_At.mat

figure(1)
plot(eta,f,'k-')
figure(2)
plot(eta,g,'k-')
figure(3)
plot(eta,k,'k-')

figure(1)
xlim([-5 15])
ylim([0 10])
xlabel('$\eta$')
ylabel('$f(\eta)$')
legend('$\Lambda_{\text{k}} = 0$','$\Lambda_{\text{k}} = 0.1$','$\Lambda_{\text{k}} = 1$','$\Lambda_{\text{k}} = 10$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_f.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');



figure(2)
xlim([-5 15])
ylim([0 10])
xlabel('$\eta$')
ylabel('$g(\eta)$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_g.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');



figure(3)
xlim([-5 15])
ylim([0 10])
xlabel('$\eta$')
ylabel('$k(\eta)$')
box on
set(0,'defaulttextinterpreter','none');
cleanfigure
matlab2tikz('SB_harden_k.tikz', 'showInfo', false, ...
    'parseStrings',false,'standalone', false, ...
    'height', '\figureheight', 'width','\figurewidth');
