% Calcula la tension normal y tangencial en la barra e para un estado u_glob dado
%
%   e      - indice de la barra a analizar
%   u_glob - vector de desplazamientos globales (2*nNodos x 1)
%   x0,y0  - coordenadas iniciales de los nodos
%   Ni,Nj  - conectividad de las barras
%   E      - modulo de elasticidad
%   L0b    - longitudes iniciales de las barras (1 x nBarras)
function [sigma, tau] = tension_barra(e, u_glob, x0, y0, Ni, Nj, E, L0b)
    ni = Ni(e); nj = Nj(e);

    % Coordenadas actuales de los extremos
    x_act_ni = x0(ni) + u_glob(2*ni-1);
    y_act_ni = y0(ni) + u_glob(2*ni);
    x_act_nj = x0(nj) + u_glob(2*nj-1);
    y_act_nj = y0(nj) + u_glob(2*nj);

    dx_act = x_act_nj - x_act_ni;
    dy_act = y_act_nj - y_act_ni;
    L_act  = sqrt(dx_act^2 + dy_act^2);

    % Tension normal (deformacion axial * modulo)
    sigma = E * (L_act - L0b(e)) / L0b(e);

    % Tension tangencial: siempre cero en barra articulada
    tau = 0;
end
