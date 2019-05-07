close all
clear all
clc

% LX-14 parameters (SI units)

rho = 1849;
G = 3.52*1e8;
c = 1130;
k = 0.439;

W = 5.95*1e6;
A = 5*1e19;
E = 2.206*1e5;
R = 8.314;

T0 = 400;
s0 = 4*1e7;
L = 3.47*1e-3;
w = 1e3;
v0 = L*w;

T = linspace(298,798,100);
plot(T,W*A*exp(-E./(R*T)))
ylim([0 10])

% Non-dimensional parameters
l = (k*s0/(rho*c*G*v0));
Gammadot0 = 1e7;% (v0/l);
t_star = s0/(G*Gammadot0);
Gamma0 = t_star*Gammadot0;
w_hat = w/Gammadot0;
rho_hat = k*Gammadot0/(c*s0);
lambda = s0^2/(rho*c*G*T0);
E_hat = E/(R*T0);
W_hat = W/(c*T0);
A_hat = A*t_star;
q0 = k*T0/l;


TR = 700/300;
A_hat0 = W_hat*A_hat/((E_hat/TR)^(1/2)*exp(E_hat/TR))

log(lambda*(1+w_hat*0.1)/(0.1 + w_hat))
