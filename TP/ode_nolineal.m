% ODE para la hipotesis no lineal (a.i): K recalculada en cada paso
% Opcionalmente incluye el amortiguador de arrastre en el nodo c
% Estado z = [u; v] con u = desplazamientos, v = velocidades (DOFs libres)
function dz = ode_nolineal(t, z, x0, y0, Ni, Nj, E, A, P_val, ...
                            DOF_lib, M_inv, nLib, nNodos, ...
                            usar_amort, rho_fl, Cd, r_esf, h_SL, nodo_c)
    nDOF = 2*nNodos;
    u = z(1:nLib);
    v = z(nLib+1:end);

    % Reconstruir desplazamientos globales
    u_glob = zeros(nDOF,1);
    u_glob(DOF_lib) = u;

    % Coordenadas actuales
    x_act = x0 + u_glob(1:2:end)';
    y_act = y0 + u_glob(2:2:end)';

    % Recalcular K con geometria deformada (hipotesis no lineal)
    [~, ct, st, kb] = props_barras(x_act, y_act, Ni, Nj, E, A);
    K_act = ensamblar_K(Ni, Nj, ct, st, kb, nNodos);
    K_red = K_act(DOF_lib, DOF_lib);

    % Fuerzas externas: cargas puntuales verticales en nodos 7 y 10
    F_glob = zeros(nDOF,1);
    F_glob(2*7)  = -P_val;
    F_glob(2*10) = -P_val;

    % Amortiguador de arrastre en nodo c (punto c del enunciado)
    if usar_amort
        vy_c = v(DOF_lib == 2*nodo_c);   % velocidad en y del nodo c
        y_c  = y_act(nodo_c);
        if ~isempty(vy_c)
            Fa = drag_force(vy_c, y_c, h_SL, r_esf, rho_fl, Cd);
            F_glob(2*nodo_c) = F_glob(2*nodo_c) + Fa;
        end
    end

    F_red = F_glob(DOF_lib);
    a = M_inv * (F_red - K_red * u);
    dz = [v; a];
end
