% Forward difference for thermal equation
function TNew = thermal_solver(stressOld,gdotOld,TOld,stressNew,gdotNew,lambda,dy,dt,N,time,q_flux)

% Generate vector
TNew = zeros(size(TOld));

% Set index
i = 2:N-1;

% Update temperature using forward difference
TNew(i) = TOld(i) + dt*((TOld(i+1) - 2*TOld(i) + TOld(i-1))/(dy^2)) ...
    + dt*lambda*((stressOld(i) + stressNew(i))/2).*((gdotOld(i) + gdotNew(i))/2);

% Heat flux disturbance
TNew(1) = TOld(1) + dt*((2*TOld(2) - 2*TOld(1) +  2*dy*q_flux(time-dt)))/(dy^2) ...
    + dt*lambda*((stressOld(1) + stressNew(1))/2).*((gdotOld(1) + gdotNew(1))/2);

% Constant temperature end
TNew(end) = 1;
end
