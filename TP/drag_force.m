% Fuerza de arrastre (Drag Force) sobre una esfera en un fluido
% Actua solo en la componente y, opuesta a la velocidad
%
%   vy_c   - velocidad en y del nodo c [m/s]
%   y_c    - coordenada y actual del centro de la esfera [m]
%   h_SL   - cota de la superficie libre del fluido [m]
%   r_esf  - radio de la esfera [m]
%   rho_fl - densidad del fluido [kg/m^3]
%   Cd     - coeficiente de arrastre [-]
function Fa = drag_force(vy_c, y_c, h_SL, r_esf, rho_fl, Cd)
    dist = h_SL - y_c;   % distancia del centro a la superficie (+ = bajo la sup.)

    if dist <= -r_esf
        % Esfera completamente fuera del fluido (por encima)
        AR = 0;
    elseif dist >= r_esf
        % Esfera completamente sumergida
        AR = pi * r_esf^2;
    else
        % Parcialmente sumergida: seccion circular al nivel del corte
        AR = pi * (r_esf^2 - dist^2);
        AR = max(AR, 0);
    end

    % Fuerza de arrastre: opuesta a la velocidad
    Fa = -0.5 * rho_fl * Cd * AR * abs(vy_c) * vy_c;
end
