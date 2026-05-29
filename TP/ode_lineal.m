% ODE para la hipotesis lineal (a.ii): K fija, fuerzas fijas
% Estado z = [u; v] con u = desplazamientos, v = velocidades (DOFs libres)
function dz = ode_lineal(t, z, K_red, F_red, M_inv, nLib)
    u = z(1:nLib);
    v = z(nLib+1:end);
    a = M_inv * (F_red - K_red * u);
    dz = [v; a];
end
