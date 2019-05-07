% First order correction to centreline temperature 
function feq = fR_eqn(f,g,AR,At,tR,fR,time)
feq = integral(@(t) f_nest(t,f,g,AR,At,tR,fR,time),time(1),time(end));
end

function expF = f_nest(t,f,g,AR,At,tR,fR,time)
ff = interp1(time,f,t);
gg = interp1(time,g,t);
expF = (pi*(time(end)-t)).^(-1/2)...
    .*(exp(t + ff - gg) + AR*exp(At*(t-tR) + ff - fR));
end


