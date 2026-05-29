% Calcula la tension normal e en la barra 'barra_a' para todos los
% instantes de tiempo hasta tF
function [sigma_v, tau_v] = calcular_tensiones_hist(t_vec, u_hist, tF, barra_a, ...
                             x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0b)
    nDOF    = 2*nNodos;
    idx     = find(t_vec <= tF);
    sigma_v = zeros(length(idx),1);
    tau_v   = zeros(length(idx),1);
    for k = 1:length(idx)
        u_glob = zeros(nDOF,1);
        u_glob(DOF_lib) = u_hist(idx(k),:)';
        [sigma_v(k), tau_v(k)] = tension_barra(barra_a, u_glob, x0, y0, Ni, Nj, E, L0b);
    end
end
