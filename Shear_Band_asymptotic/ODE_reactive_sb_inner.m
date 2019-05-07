% ODEs for reaction stage inner layer
function ODE = ODE_reactive_sb_inner(~,RHS,beta1,beta2,beta3,gammadot0,lambda,A0,Omega,t_R,s_R,T1tR,s1tR)

ODE = [lambda*s_R*gammadot0*exp(beta3*t_R + beta1*T1tR + beta2*s1tR + beta1*RHS(1) + beta2*RHS(2)) + Omega*A0*exp(RHS(1))
    -gammadot0*exp(beta3*t_R + beta1*T1tR + beta2*s1tR + beta1*RHS(1) + beta2*RHS(2))];
end
    