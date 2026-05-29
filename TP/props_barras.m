% Propiedades geometricas de las barras dada configuracion actual de nodos
function [L_b, ct, st, kb] = props_barras(x, y, Ni, Nj, E, A)
    nB = length(Ni);
    L_b = zeros(1,nB); ct = zeros(1,nB); st = zeros(1,nB); kb = zeros(1,nB);
    for e = 1:nB
        dx     = x(Nj(e)) - x(Ni(e));
        dy     = y(Nj(e)) - y(Ni(e));
        L_b(e) = sqrt(dx^2 + dy^2);
        ct(e)  = dx / L_b(e);
        st(e)  = dy / L_b(e);
        kb(e)  = E * A / L_b(e);
    end
end
