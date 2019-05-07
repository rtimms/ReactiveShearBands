% Forward difference for thermal equation
function TNew = thermal_solver_reactive(stressOld,gdotOld,TOld,stressNew,gdotNew,TPred,lambda,dy,dt,N,time,q_flux,Omega_hat,A0_hat,E_hat,T_R)

% Generate vector
TNew = TOld;

% Set index
i = 2:N-1;

% Update temperature using forward difference
TNew(i) = TOld(i) + dt*((TOld(i+1) - 2*TOld(i) + TOld(i-1))/(dy^2)) ...
    + dt*lambda*((stressOld(i) + stressNew(i))/2).*((gdotOld(i) + gdotNew(i))/2) ...
    + dt*Omega_hat*A0_hat*(E_hat^(1/2)/T_R)*exp(E_hat/T_R - 2*E_hat./(TOld(i) + TPred(i)));

% Heat flux disturbance
TNew(1) = TOld(1) + dt*((2*TOld(2) - 2*TOld(1) +  2*dy*q_flux(time-dt)))/(dy^2) ...
    + dt*lambda*((stressOld(1) + stressNew(1))/2).*((gdotOld(1) + gdotNew(1))/2) ...
    + dt*Omega_hat*A0_hat*(E_hat^(1/2)/T_R)*exp(E_hat/T_R - 2*E_hat./(TOld(1) + TPred(1)));

% Constant temperature end
TNew(end) = 1;

end
