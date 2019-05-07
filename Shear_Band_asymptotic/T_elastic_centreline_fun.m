% Centreline elastic temperature
function Te_centreline = T_elastic_centreline_fun(delta,time,h)

Te_centreline = real( 1 + delta*integral(@(t) (pi*(time-t)).^(-1/2).*h(t),0,time) );

end
