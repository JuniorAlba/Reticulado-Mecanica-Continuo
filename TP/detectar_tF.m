% Detecta el tiempo de colapso tF como el primer instante en que el area
% con signo de algun triangulo cambia de signo (inversion topologica)
%
%   t_vec      - vector de tiempos de la solucion ODE
%   u_hist     - matriz de desplazamientos (filas = pasos, cols = DOFs libres)
%   x0, y0     - coordenadas iniciales de los nodos
%   triangulos - matriz Ntri x 3 con indices de nodos de cada triangulo
%   DOF_lib    - indices de DOFs libres
%   nNodos     - numero de nodos
%   tF_max     - tiempo maximo de simulacion (valor por defecto si no hay colapso)
function tF = detectar_tF(t_vec, u_hist, x0, y0, triangulos, DOF_lib, nNodos, tF_max)
    nDOF = 2*nNodos;
    nTri = size(triangulos,1);
    tF   = tF_max;

    % Areas con signo en la configuracion inicial
    A0 = zeros(nTri,1);
    for tr = 1:nTri
        ni=triangulos(tr,1); nj=triangulos(tr,2); nk=triangulos(tr,3);
        A0(tr) = area_triangulo(x0(ni),y0(ni), x0(nj),y0(nj), x0(nk),y0(nk));
    end
    signo0 = sign(A0);

    for step = 1:length(t_vec)
        u_glob = zeros(nDOF,1);
        u_glob(DOF_lib) = u_hist(step,:)';
        x_act = x0 + u_glob(1:2:end)';
        y_act = y0 + u_glob(2:2:end)';

        for tr = 1:nTri
            ni=triangulos(tr,1); nj=triangulos(tr,2); nk=triangulos(tr,3);
            A_act = area_triangulo(x_act(ni),y_act(ni), x_act(nj),y_act(nj), x_act(nk),y_act(nk));
            if sign(A_act) ~= signo0(tr)
                tF = t_vec(step);
                fprintf('Cambio de signo en triangulo %d en t=%.4f s\n', tr, tF);
                return;
            end
        end
    end
end
