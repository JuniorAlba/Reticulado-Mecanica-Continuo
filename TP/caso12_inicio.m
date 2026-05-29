% =========================================================
%  CASO 12 - Reticulado plano 2D
%  Trabajo Practico 1
%  Paso 1: Geometria, rigidez elemental y masas nodales
% =========================================================

clear; clc; close all;

% ---------------------------------------------------------
% PARAMETROS DEL PROBLEMA
% ---------------------------------------------------------
E   = 600;      % Modulo de elasticidad longitudinal
A   = 0.25;     % Area de seccion transversal
rho = 1;        % Densidad de las barras
P   = 3.2;      % Carga uniforme en el tiempo

% ---------------------------------------------------------
% PASO 1 y 2: GEOMETRIA Y CONECTIVIDAD DESDE ARCHIVO DXF
% Lectura automatica con f_LectDxf (provisto por la catedra)
% ---------------------------------------------------------
dxf_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'Consigna', 'dxf', 'dxf', 'Estr12.dxf');
[c_Line, ~, c_Cir, ~, ~] = f_LectDxf(dxf_path);

% Nodos: centros de los circulos (en el orden en que aparecen en el DXF)
m_Cir = cell2mat(c_Cir(:,1));   % [x, y, radio] para cada circulo
x     = round(m_Cir(:,1)' * 1e6) / 1e6;  % coordenadas x (1 x nNodos)
y     = round(m_Cir(:,2)' * 1e6) / 1e6;  % coordenadas y (1 x nNodos)
nNodos = length(x);

% Barras: lineas cuyos dos extremos coincidan con nodos conocidos
% (las lineas de cota/dimension del DXF tienen extremos fuera de la red)
m_Lin = cell2mat(c_Line(:,1));  % [Xi Yi Zi Xj Yj Zj]
tol   = 0.5;                    % tolerancia geometrica [m]
Ni = []; Nj = [];
for e_dxf = 1:size(m_Lin, 1)
    xi = m_Lin(e_dxf,1);  yi = m_Lin(e_dxf,2);
    xj = m_Lin(e_dxf,4);  yj = m_Lin(e_dxf,5);
    di = sqrt((x - xi).^2 + (y - yi).^2);
    dj = sqrt((x - xj).^2 + (y - yj).^2);
    [min_di, ni] = min(di);
    [min_dj, nj] = min(dj);
    if min_di < tol && min_dj < tol && ni ~= nj
        Ni(end+1) = ni;  %#ok<AGROW>
        Nj(end+1) = nj;  %#ok<AGROW>
    end
end
nBarras = length(Ni);

fprintf('Geometria leida desde DXF: %d nodos, %d barras\n', nNodos, nBarras);

% Verificacion: imprimir coordenadas
fprintf('\n--- Coordenadas de los nodos ---\n');
fprintf('Nodo  x     y\n');
for i = 1:nNodos
    fprintf('  %2d  %5.1f  %5.1f\n', i, x(i), y(i));
end

% ---------------------------------------------------------
% PASO 3: PROPIEDADES GEOMETRICAS DE CADA BARRA
% ---------------------------------------------------------
fprintf('\n--- Propiedades geometricas de las barras ---\n');
fprintf('Barra  Ni  Nj    dx      dy      L       cos_t   sin_t   k\n');

L_bar  = zeros(1, nBarras);
cos_t  = zeros(1, nBarras);
sin_t  = zeros(1, nBarras);
k_bar  = zeros(1, nBarras);

for e = 1:nBarras
    ni = Ni(e);
    nj = Nj(e);
    dx = x(nj) - x(ni);
    dy = y(nj) - y(ni);
    L_bar(e)  = sqrt(dx^2 + dy^2);
    cos_t(e)  = dx / L_bar(e);
    sin_t(e)  = dy / L_bar(e);
    k_bar(e)  = E * A / L_bar(e);
    fprintf('  %2d    %2d  %2d   %6.2f  %6.2f  %6.3f  %6.3f  %6.3f  %7.3f\n', ...
        e, ni, nj, dx, dy, L_bar(e), cos_t(e), sin_t(e), k_bar(e));
end

% ---------------------------------------------------------
% PASO 4: MASA CONCENTRADA EN CADA NODO
% Cada barra aporta la mitad de su masa a cada extremo
% masa_barra = rho * A * L
% ---------------------------------------------------------
m_nodo = zeros(1, nNodos);

for e = 1:nBarras
    masa_barra = rho * A * L_bar(e);
    m_nodo(Ni(e)) = m_nodo(Ni(e)) + 0.5 * masa_barra;
    m_nodo(Nj(e)) = m_nodo(Nj(e)) + 0.5 * masa_barra;
end

fprintf('\n--- Masa concentrada en cada nodo ---\n');
fprintf('Nodo   Masa\n');
for i = 1:nNodos
    fprintf('  %2d   %.4f\n', i, m_nodo(i));
end
fprintf('Masa total: %.4f\n', sum(m_nodo));

% ---------------------------------------------------------
% PASO 5: ENSAMBLAJE DE LA MATRIZ DE RIGIDEZ GLOBAL K
% Cada nodo tiene 2 DOF: (u_x, u_y)
% DOF del nodo i: (2i-1) en x, (2i) en y
% ---------------------------------------------------------
nDOF = 2 * nNodos;   % = 30
K = zeros(nDOF, nDOF);

for e = 1:nBarras
    ni = Ni(e);
    nj = Nj(e);
    c  = cos_t(e);
    s  = sin_t(e);
    ke = k_bar(e);

    % Matriz de rigidez local 4x4 en coordenadas globales
    Ke = ke * [ c^2,  c*s, -c^2, -c*s;
                c*s,  s^2, -c*s, -s^2;
               -c^2, -c*s,  c^2,  c*s;
               -c*s, -s^2,  c*s,  s^2 ];

    % DOF globales de este elemento
    dofs = [2*ni-1, 2*ni, 2*nj-1, 2*nj];

    % Ensamblaje (suma en posiciones correspondientes)
    K(dofs, dofs) = K(dofs, dofs) + Ke;
end

fprintf('\nMatriz K global ensamblada: %dx%d\n', nDOF, nDOF);

% ---------------------------------------------------------
% PASO 6: CONDICIONES DE BORDE CINEMATICAS
% Nodos 1 y 15 estan empotrados: todos sus DOF = 0
% DOFs restringidos: 1,2 (nodo 1) y 29,30 (nodo 15)
% ---------------------------------------------------------
DOF_restringidos = [1, 2, 29, 30];
DOF_libres = setdiff(1:nDOF, DOF_restringidos);

% Submatriz de rigidez reducida (solo DOFs libres)
K_red = K(DOF_libres, DOF_libres);

% Subvector de masas reducido
% Matriz de masa diagonal (2 DOF por nodo, misma masa en x e y)
M_diag = zeros(nDOF, 1);
for i = 1:nNodos
    M_diag(2*i-1) = m_nodo(i);   % DOF en x
    M_diag(2*i)   = m_nodo(i);   % DOF en y
end
M_red = diag(M_diag(DOF_libres));   % Matriz diagonal reducida

fprintf('\nDOFs libres: %d (de %d totales)\n', length(DOF_libres), nDOF);
fprintf('DOFs restringidos: %s\n', mat2str(DOF_restringidos));

% ---------------------------------------------------------
% PASO 7: VECTOR DE FUERZAS EXTERNAS
% Carga P aplicada en nodo 7 (DOF 14, dir y negativa)
% Carga P aplicada en nodo 10 (DOF 20, dir y negativa)
% ---------------------------------------------------------
F_global = zeros(nDOF, 1);
F_global(2*7)  = -P;   % nodo 7, direccion y
F_global(2*10) = -P;   % nodo 10, direccion y

F_red = F_global(DOF_libres);

% ---------------------------------------------------------
% PASO 8: SOLUCION ESTATICA (verificacion)
% K_red * u_red = F_red
% ---------------------------------------------------------
u_red = K_red \ F_red;

% Reconstruir vector completo de desplazamientos
u_global = zeros(nDOF, 1);
u_global(DOF_libres) = u_red;

fprintf('\n--- Desplazamientos estaticos (verificacion) ---\n');
fprintf('Nodo   ux          uy\n');
for i = 1:nNodos
    fprintf('  %2d   %10.6f  %10.6f\n', i, u_global(2*i-1), u_global(2*i));
end

% ---------------------------------------------------------
% PASO 9: FRECUENCIAS NATURALES DEL SISTEMA
% Problema de autovalores: K*phi = omega^2 * M * phi
% ---------------------------------------------------------
[V, D] = eig(K_red, M_red);
omega2 = diag(D);
omega2 = sort(omega2);                 % ordenar de menor a mayor
omega  = sqrt(abs(omega2));            % frecuencias naturales [rad/s]
freq   = omega / (2*pi);              % frecuencias [Hz]

fprintf('\n--- Primeras 6 frecuencias naturales ---\n');
fprintf('Modo  omega [rad/s]   f [Hz]\n');
for i = 1:min(6, length(omega))
    fprintf('  %2d   %10.4f    %10.4f\n', i, omega(i), freq(i));
end

fprintf('\n=== Fin del Paso 1 y 2. Listo para dinamica. ===\n');

% ---------------------------------------------------------
% GRAFICA: Estructura inicial con numeracion de nodos y barras
% ---------------------------------------------------------
figure('Name','Estructura inicial - Caso 12','NumberTitle','off');
hold on; axis equal; grid on;

% Barras
for e = 1:nBarras
    xm = (x(Ni(e)) + x(Nj(e))) / 2;
    ym = (y(Ni(e)) + y(Nj(e))) / 2;
    plot([x(Ni(e)) x(Nj(e))], [y(Ni(e)) y(Nj(e))], 'b-', 'LineWidth', 1.5);
    text(xm, ym+0.4, sprintf('%d',e), 'FontSize', 7, 'Color',[0 0.5 0], ...
        'HorizontalAlignment','center');
end

% Nodos libres
plot(x, y, 'ko', 'MarkerFaceColor','w', 'MarkerSize', 8);
for i = 1:nNodos
    text(x(i)+0.5, y(i)+1.0, sprintf('%d',i), 'FontSize', 9, ...
        'Color','r', 'FontWeight','bold');
end

% Empotramientos (nodos 1 y 15)
plot(x([1,15]), y([1,15]), 'rs', 'MarkerFaceColor','r', 'MarkerSize', 10);

% Cargas (nodos 7 y 10)
quiver(x(7),  y(7)+3,  0, -2.5, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.8);
quiver(x(10), y(10)+3, 0, -2.5, 0, 'k', 'LineWidth', 2, 'MaxHeadSize', 0.8);
text(x(7)+0.5,  y(7)+3.5,  'P', 'FontSize', 9, 'FontWeight','bold');
text(x(10)+0.5, y(10)+3.5, 'P', 'FontSize', 9, 'FontWeight','bold');

xlabel('x [m]'); ylabel('y [m]');
title('Reticulado Caso 12 - Configuracion inicial');
legend({'Barras','Nodos','Empotramientos','Carga P'}, 'Location','northeast');
