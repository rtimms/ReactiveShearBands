% solve until f or g become O(1/eps) which indicates a shear band
 close all
 clear all

Ap_vec = [1 5 10:5:100];
AR_vec = 1:30;
[Ap_mesh,AR_mesh] =  meshgrid(Ap_vec,AR_vec);

param_ind = str2num(getenv('LSB_JOBINDEX'));


% % Set parameters
Ap = Ap_mesh(param_ind);
AR = AR_mesh(param_ind);
Ap = 10;
AR = 1;
At = 0.01;

filename = sprintf('output_NEW/fgreactive_Ap%.0f_AR%.0f_At%.0f',Ap*100,AR*100,At*100);

% Initial f,g and eta
f0 = 0;
g0 = 0;
eta0 = -30;
eta_Final = 20;

% Set tolerance
tol = 1e-2;                        % Iteration tolerance
d_eta = 1e-1;                      % Initial timestep
stepTol = 1e-2;                    % Set relative change for step size control
it_Max_1 = 10;                     % Maxmimum iterations in fg_reactive_solver
it_Max_2 = 10;                      % Maxmimum iterations in eta_R solution

% Set epsilon
eps = 1e-2;
fMax = 20;
gMax = 20;
f_R = 5;

% Iterate to find  eta_R and f_R
tol_R = 1e-3;
eta_R = 1e6;
iterate = 0;

% Iterate eta_R if AR non-zero
if AR == 0
    
    % Compute f and g for inert case
    [f,g,eta,converge] = fg_reactive_solver(f0,g0,eta0,eta_Final,0,0,Ap,0,0,fMax,gMax,tol,d_eta,stepTol,it_Max_1);
    
else
    
    % Initial guess for eta_R using inert solution
    [~,~,eta,~] = fg_reactive_solver(f0,g0,eta0,eta_Final,0,0,Ap,0,0,f_R,1e6,tol,d_eta,stepTol,it_Max_1);
    eta_R_New = eta(end)
    
    % Iterate to find eta_R, and compute f and g
    while abs(eta_R_New - eta_R) > tol_R;
        
        % Break if too many iterations
        iterate = iterate + 1;
        if iterate > it_Max_2
            sprintf('Maximum number of iterations exceeded in computing eta_R! Tolerance exceeded by %f',abs( abs(eta_R_New-eta_R) - tol))
            converge = 0;
            break
        end
        
        % Save old eta_R
        eta_R = eta_R_New;
        
        % Compute f and g for current value of eta_R
        [f,~,eta,~] = fg_reactive_solver(f0,g0,eta0,eta_Final,eta_R,f_R,Ap,AR,At,f_R,1e6,tol,d_eta,stepTol,it_Max_1);
        [~,idx] = min(abs(f-f_R));
        eta_R_New = eta(idx)
        
    end
    
    % Run with converged f_R and eta_R
    [f,g,eta,converge] = fg_reactive_solver(f0,g0,eta0,eta_Final,eta_R,f_R,Ap,AR,At,fMax,gMax,tol,d_eta,stepTol,it_Max_1);
   
end

save(filename)

