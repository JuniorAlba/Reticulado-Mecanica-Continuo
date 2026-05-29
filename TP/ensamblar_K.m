% Ensamblaje de la matriz de rigidez global K
function K = ensamblar_K(Ni, Nj, ct, st, kb, nNodos)
    nDOF = 2*nNodos;
    K = zeros(nDOF, nDOF);
    for e = 1:length(Ni)
        c = ct(e); s = st(e); k = kb(e);
        Ke = k * [ c^2,  c*s, -c^2, -c*s;
                   c*s,  s^2, -c*s, -s^2;
                  -c^2, -c*s,  c^2,  c*s;
                  -c*s, -s^2,  c*s,  s^2 ];
        dofs = [2*Ni(e)-1, 2*Ni(e), 2*Nj(e)-1, 2*Nj(e)];
        K(dofs,dofs) = K(dofs,dofs) + Ke;
    end
end
