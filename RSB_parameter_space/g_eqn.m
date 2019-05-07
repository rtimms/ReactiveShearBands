% First order correction to centreline stress
function geq = g_eqn(f,g,A,time)

geq = A*exp(time(end)+f-g);
end