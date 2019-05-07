% Newton iteration to update mechanical equations
function [velNew,sNew,gdotNew,gNew,converge] = ...
    mechanical_solver_exp(velOld,sOld,gdotOld,gOld,TOld,rho_hat,S,omega,gammadot0,eps,beta1,beta2,T_p,s_p,y,dy,dt,N)

% Plastic strain-rate
plastic_rate = @(s,T) gammadot0*eps^(-1/2)*exp(-eps^(-1)*(beta1*(T_p - T) + beta2*(s_p - s)));
dplastic_rate_ds = @(s,T) gammadot0*beta2*eps^(-3/2)*exp(-eps^(-1)*(beta1*(T_p - T) + beta2*(s_p - s)));

% Generate vectors
velNew = velOld;
sNew = sOld;
gdotNew = gdotOld;
%gNew = gOld;

% Calculate velocities
v_minus = velOld - gdotOld*dy/2;
v_plus = velOld + gdotOld*dy/2;

% Set index
i = 2:N-1;

% Reset tolerance
iterate = 0;
it_max = 5;
tol = 1;
converge = 1;

while tol > 1e-6
    
    % Count iteration
    iterate = iterate + 1;
    if iterate > it_max
        converge = 0;
        break
    end
    
    % Compute plastic rate 
    r = plastic_rate(sNew,TOld);
    drds = dplastic_rate_ds(sNew,TOld);
    
    % Jacobian
    J = [speye(N-2) (rho_hat*S*dy/2)*speye(N-2)
        diag(-drds(i)) speye(N-2)];
    
    % Newton step for stress and strain-rate
    delta_sol = -J\([ sNew(i) - (sOld(i-1) + sOld(i+1))/2 ...
        - (rho_hat*S/2)*(v_minus(i+1) - v_plus(i-1)) ...
        + (rho_hat*S*dy/2)*gdotNew(i);...
        gdotNew(i) - r(i)] );
    
    sNew(i) = sNew(i) + delta_sol(1:N-2);
    gdotNew(i) = gdotNew(i) + delta_sol(N-1:end);
    
    % Intermediate strain
    % gNew = gOld + (dt/2)*(gdotOld + gdotNew);
    
    % Calculate difference
    tol = norm(delta_sol,inf);
    
end

% Boundary values
%sNew(1) = sOld(2) - rho_hat*S*(velOld(1) + (dy/2)*gdotNew(1) - v_minus(2));
%sNew(end) = sOld(end-1) + rho_hat*S*(velOld(end) - (dy/2)*gdotNew(end) - v_plus(end-1));

sNew(1) = sNew(2);
sNew(end) = sNew(end-1);

gdotNew(1) =  r(1);
gdotNew(end) =  r(end);

% Update strain
gNew = gOld + (dt/2)*(gdotOld + gdotNew);

% Update velocity
v_minusNew(i) = v_plus(i-1) + (1/(rho_hat*S))*(sNew(i) - sOld(i-1));
v_plusNew(i) = v_minus(i+1) - (1/(rho_hat*S))*(sNew(i) - sOld(i+1));

velNew(i) = (v_minusNew(i) + v_plusNew(i))/2;
velNew(1) = omega*y(1);
velNew(end) = omega*y(end);
end
