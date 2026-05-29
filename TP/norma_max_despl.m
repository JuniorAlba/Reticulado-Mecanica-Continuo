% Calcula la norma maxima del vector de desplazamiento nodal
% hasta el tiempo tF, y el instante en que ocurre
function [norm_max, t_norm_max] = norma_max_despl(t_vec, u_hist, DOF_lib, nNodos, tF)
    nDOF  = 2*nNodos;
    idx   = find(t_vec <= tF);
    norms = zeros(length(idx),1);
    for k = 1:length(idx)
        u_glob = zeros(nDOF,1);
        u_glob(DOF_lib) = u_hist(idx(k),:)';
        % Norma del vector de desplazamiento para cada nodo
        u_nodos     = reshape(u_glob, 2, nNodos);
        normas_nodo = sqrt(u_nodos(1,:).^2 + u_nodos(2,:).^2);
        norms(k)    = max(normas_nodo);
    end
    [norm_max, idx_max] = max(norms);
    t_norm_max = t_vec(idx(idx_max));
end
