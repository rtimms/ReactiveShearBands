function [eta,f,g] = fg_solver_test

% Initial f,g and timeVec
k = 1;
f(1) = 0;
g(1) = 0;
eta(1) = -15;

% Parameters
Lambda_p = 1;
Lambda_R = 0;
Lambda_t = 0;
eta_R = 10.3393;
fR = 1.4;

% Set tolerance
tol = 1e-3;                         % Iteration tolerance
d_eta = 1e-2;                       % Initial timestep
stepTol = 1e-2;                     % Set relative change for step size control
converge = 1;

while f(end) < 20 && converge == 1
    
    % Increase timestepping index
    k = k+1;
    
    % Loop to repeat timstep if temperature tolerance not met
    proceed = 0;
    
    while proceed == 0
        
        % Update time
        eta(k) = eta(k-1) + d_eta;
        eta(k)
        
        % Reset tolerance
        diff_f = tol*10;
        diff_g = tol*10;
        iterate = 0;
        
        % Use previous result as initial guess
        fNew = f(k-1);
        gNew = g(k-1);
        
        
        while max(diff_f,diff_g)  > tol && converge == 1;
            
            iterate = iterate + 1;
            
            % Break if maximum iterations exceeded
            if iterate > 10
                f(k) = fNew;
                g(k) = gNew;
                sprintf('Maximum number of iterations exceeded in fg_solver_EPreactive !')
                converge = 0;
            end
            
            % Save old values for computation of diff_f and diff_g
            fOld = fNew;
            gOld = gNew;
            
            % Calculate updated f and g
            fNew = f_EPreactive([f(1:k-1) fNew],[g(1:k-1) gNew],Lambda_R,Lambda_t,eta_R,fR,eta);
            gNew = gOld...
                - (gOld - g_EPreactive(fOld,gOld,Lambda_p,eta(end)))...
                /(1 + g_EPreactive(fOld,gOld,Lambda_p,eta(end)));
            
            % Calculate change in iterated solution
            diff_f = norm(fNew-fOld);
            diff_g = norm(gNew-gOld);
            
        end
        
        % Check if time step needs reducing and update results
        fChange = (fNew - f(k-1));
        gChange = (gNew - g(k-1));
        if any(fChange > stepTol) || any(gChange > stepTol)
            % Reduce timestep
            d_eta = d_eta/2;
        else
            % Update results
            proceed = 1;
            f(k) = fNew;
            g(k) = gNew;

        end
        
    end
    
end

figure(1)
plot(eta,f,'rx')
figure(2)
plot(eta,g,'bo')
figure(3)
plot(eta,f,'rx',eta,g,'bo',eta,exp(eta),'g-')

end

%% Auxilliary

function f = f_EPreactive(f,g,AR,At,tR,fR,time)
f = integral(@(t) f_nest(t,f,g,AR,At,tR,fR,time),time(1),time(end));
end

function g = g_EPreactive(f,g,Ap,time)
g = Ap*exp(time + f - g);
end

function expF = f_nest(t,f,g,AR,At,tR,fR,time)
ff = interp1(time,f,t);
gg = interp1(time,g,t);
expF = (pi*(time(end)-t)).^(-1/2)...
    .*(exp(t + ff - gg) + AR*exp(At*(t-tR) + ff - fR));
end
