% =========================================================
%  CASO 12 - Reticulado plano 2D
%  Dinamica: hipotesis lineal (a.ii) y no lineal (a.i)
%  Deteccion de tF, tensiones, amortiguador Drag Force
% =========================================================

clear; clc; close all;

% ---------------------------------------------------------
% PARAMETROS DEL PROBLEMA
% ---------------------------------------------------------
E      = 600;     % Modulo de elasticidad longitudinal
A      = 0.25;    % Area de seccion transversal
rho    = 1;       % Densidad de las barras
P_val  = 3.2;     % Carga uniforme
tF_max = 50;      % Tiempo maximo de simulacion [s]

% Amortiguador (punto c)
rho_fl = 3.0;     % Densidad del fluido
Cd     = 0.47;    % Coeficiente de arrastre (esfera, regimen 1e3<Re<2e5)
r_esf  = 1.5;     % Radio de la esfera
h_SL   = 30.0;    % Cota superficie libre del fluido
nodo_c = 8;       % Nodo donde actua el amortiguador
barra_a = 5;      % Barra a analizar
nodo_b  = 4;      % Nodo b a monitorear

% ---------------------------------------------------------
% COORDENADAS DE LOS NODOS
% ---------------------------------------------------------
x0 = [  0, 10,  5, 20, 15, 25, 20, 30, 35, 40, 45, 40, 50, 55, 60 ];
y0 = [  0,  0, 5,  0, 15, 15, 20, 30, 15, 20, 15,  0,  0, 5,  0 ];

nNodos = numel(x0);

if numel(y0) ~= nNodos
    error('x0 e y0 deben tener la misma cantidad de nodos.');
end

% ---------------------------------------------------------
% CONECTIVIDAD
% ---------------------------------------------------------
Ni = [ 1,  1,  2,  4,  4,  5,  2,  5,  6,  3,  5,  7, 15, 15, 14, 11, 10, 13, 12, 12, 11, 13, 11,  9,  6,  9 ];
Nj = [ 2,  3,  4,  6,  5,  2,  3,  6,  7,  5,  7,  8, 13, 14, 11, 10,  8, 12,  9, 11, 13, 14,  9, 10,  8,  8 ];
nBarras = numel(Ni);

if numel(Nj) ~= nBarras
    error('Ni y Nj deben tener la misma cantidad de barras.');
end

if any([Ni, Nj] < 1) || any([Ni, Nj] > nNodos)
    error('La conectividad Ni/Nj referencia nodos fuera de rango 1..%d.', nNodos);
end

% ---------------------------------------------------------
% TRIANGULOS (para deteccion de cruce de barras)
% Cada triangulo definido por 3 nodos en orden antihorario inicial
% Se identifican visualmente de la figura
% ---------------------------------------------------------
triangulos = [
    1,  2,  3;      %1
    2,  4,  5;      %
    2,  5,  3;      %
    4,  6,  5;      %
    5,  6,  7;      %
    7,  6,  8;      %
    9, 10,  8;
    9, 11, 10;
    9, 12, 11;
    12, 13, 11;
    11,13,14;
    13,15,14;


];
% Nota: verificar orientacion de cada triangulo con area_signo > 0 al t=0

if any(triangulos(:) < 1) || any(triangulos(:) > nNodos)
    error('La matriz triangulos referencia nodos fuera de rango 1..%d.', nNodos);
end

% ---------------------------------------------------------
% MASAS NODALES (calculadas con geometria inicial)
% ---------------------------------------------------------
[L0, ~, ~, ~] = props_barras(x0, y0, Ni, Nj, E, A);
m_nodo = zeros(1, nNodos);
for e = 1:nBarras
    mb = rho * A * L0(e);
    m_nodo(Ni(e)) = m_nodo(Ni(e)) + 0.5*mb;
    m_nodo(Nj(e)) = m_nodo(Nj(e)) + 0.5*mb;
end

% Vector de masa para todos los DOF
nDOF = 2*nNodos;
M_vec = zeros(nDOF,1);
for i = 1:nNodos
    M_vec(2*i-1) = m_nodo(i);
    M_vec(2*i)   = m_nodo(i);
end

% ---------------------------------------------------------
% CONDICIONES DE BORDE
% Nodos 1 y 15 empotrados: DOF 1,2,29,30 fijos
% ---------------------------------------------------------
DOF_rest  = [1, 2, 29, 30];
DOF_lib   = setdiff(1:nDOF, DOF_rest);
nLib      = length(DOF_lib);

M_inv = diag(1 ./ M_vec(DOF_lib));   % inversa de M reducida (diagonal)

% Areas iniciales de los triangulos (para detectar cambio de signo)
nTri = size(triangulos,1);
A0_tri = zeros(nTri,1);
for t = 1:nTri
    ni = triangulos(t,1); nj = triangulos(t,2); nk = triangulos(t,3);
    A0_tri(t) = area_triangulo(x0(ni),y0(ni), x0(nj),y0(nj), x0(nk),y0(nk));
end
fprintf('Areas iniciales de los triangulos:\n');
disp(A0_tri');

% ---------------------------------------------------------
% FUNCION ODE - HIPOTESIS a.ii (LINEAL, pequenos desplazamientos)
% K calculada UNA SOLA VEZ con geometria inicial
% Fuerzas mantienen direccion fija
% ---------------------------------------------------------
[L0b, ct0, st0, kb0] = props_barras(x0, y0, Ni, Nj, E, A);
K0     = ensamblar_K(Ni, Nj, ct0, st0, kb0, nNodos);
K0_red = K0(DOF_lib, DOF_lib);

% Vector de fuerzas estaticas (P hacia abajo en nodos 7 y 10)
F_est = zeros(nDOF,1);
F_est(2*7)  = -P_val;
F_est(2*10) = -P_val;
F_red = F_est(DOF_lib);

% =========================================================
%  SIMULACION a.ii: LINEAL
% =========================================================
fprintf('\n=== Simulacion a.ii: Hipotesis lineal ===\n');

z0    = zeros(2*nLib, 1);
tspan = [0, tF_max];
opts  = odeset('RelTol',1e-6, 'AbsTol',1e-8, 'MaxStep',0.05);

[t_lin, Z_lin] = ode45(@(t,z) ode_lineal(t, z, K0_red, F_red, M_inv, nLib), ...
                        tspan, z0, opts);

u_hist_lin = Z_lin(:, 1:nLib);   % desplazamientos en DOFs libres

% Detectar tF lineal
tF_lin = detectar_tF(t_lin, u_hist_lin, x0, y0, triangulos, DOF_lib, nNodos, tF_max);
fprintf('tF lineal = %.4f s\n', tF_lin);

% =========================================================
%  SIMULACION a.i: NO LINEAL (sin amortiguador)
% =========================================================
fprintf('\n=== Simulacion a.i: Hipotesis no lineal (sin amortiguador) ===\n');

[t_nl, Z_nl] = ode45(@(t,z) ode_nolineal(t, z, x0, y0, Ni, Nj, E, A, P_val, ...
                     DOF_lib, M_inv, nLib, nNodos, ...
                     false, rho_fl, Cd, r_esf, h_SL, nodo_c), ...
                     tspan, z0, opts);

u_hist_nl = Z_nl(:, 1:nLib);
tF_nl = detectar_tF(t_nl, u_hist_nl, x0, y0, triangulos, DOF_lib, nNodos, tF_max);
fprintf('tF no lineal = %.4f s\n', tF_nl);

% =========================================================
%  SIMULACION a.i + AMORTIGUADOR (punto c)
% =========================================================
fprintf('\n=== Simulacion con amortiguador (punto c) ===\n');

[t_am, Z_am] = ode45(@(t,z) ode_nolineal(t, z, x0, y0, Ni, Nj, E, A, P_val, ...
                     DOF_lib, M_inv, nLib, nNodos, ...
                     true, rho_fl, Cd, r_esf, h_SL, nodo_c), ...
                     tspan, z0, opts);

u_hist_am = Z_am(:, 1:nLib);
tF_am = detectar_tF(t_am, u_hist_am, x0, y0, triangulos, DOF_lib, nNodos, tF_max);
fprintf('tF con amortiguador = %.4f s\n', tF_am);

% =========================================================
%  POST-PROCESO
% =========================================================

% --- Tension en barra a=5 a lo largo del tiempo ---
tF_min = min(tF_lin, tF_nl);

[sig_lin, tau_lin] = calcular_tensiones_hist(t_lin, u_hist_lin, tF_lin, barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);
[sig_nl,  tau_nl]  = calcular_tensiones_hist(t_nl,  u_hist_nl,  tF_nl,  barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);
[sig_am,  tau_am]  = calcular_tensiones_hist(t_am,  u_hist_am,  tF_am,  barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);

% --- Coordenada actual del nodo b=4 ---
DOF_b_x = find(DOF_lib == 2*nodo_b-1);
DOF_b_y = find(DOF_lib == 2*nodo_b);

idx_min_lin = find(t_lin <= tF_lin);
idx_min_nl  = find(t_nl  <= tF_nl);
idx_min_am  = find(t_am  <= tF_am);

y_nodo_b_lin = y0(nodo_b) + u_hist_lin(idx_min_lin, DOF_b_y);
y_nodo_b_nl  = y0(nodo_b) + u_hist_nl( idx_min_nl,  DOF_b_y);
y_nodo_b_am  = y0(nodo_b) + u_hist_am( idx_min_am,  DOF_b_y);

% --- Norma maxima del vector desplazamiento ---
[nm_lin, t_nm_lin] = norma_max_despl(t_lin, u_hist_lin, DOF_lib, nNodos, tF_lin);
[nm_nl,  t_nm_nl]  = norma_max_despl(t_nl,  u_hist_nl,  DOF_lib, nNodos, tF_nl);
[nm_am,  t_nm_am]  = norma_max_despl(t_am,  u_hist_am,  DOF_lib, nNodos, tF_am);

fprintf('\n--- Resumen de resultados ---\n');
fprintf('                    Lineal (a.ii)   No lineal (a.i)   Con amortiguador\n');
fprintf('tF [s]              %10.4f      %10.4f        %10.4f\n', tF_lin, tF_nl, tF_am);
fprintf('Norma max despl     %10.4f      %10.4f        %10.4f\n', nm_lin, nm_nl, nm_am);
fprintf('Tiempo norma max    %10.4f      %10.4f        %10.4f\n', t_nm_lin, t_nm_nl, t_nm_am);

% =========================================================
%  GRAFICAS
% =========================================================

figure('Name','Tension normal barra a=5','NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), sig_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   sig_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   sig_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
xlabel('Tiempo [s]'); ylabel('\sigma [Pa]');
title('Tension normal en barra a=5');
legend('Location','northeast'); grid on;

figure('Name','Tension tangencial barra a=5','NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), tau_lin, 'b-', 'LineWidth', 1.5);
xlabel('Tiempo [s]'); ylabel('\tau [Pa]');
title('Tension tangencial en barra a=5 (debe ser cero: barra articulada)');
grid on;

figure('Name','Coordenada actual nodo b=4','NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), y_nodo_b_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   y_nodo_b_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   y_nodo_b_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
xlabel('Tiempo [s]'); ylabel('y_{nodo 4} [m]');
title('Coordenada actual del nodo b=4');
legend('Location','northeast'); grid on;

% --- Vista conjunta para verificar tiempos finales ---
figure('Name','Comparacion a.i, a.ii y reticulado amortiguado','NumberTitle','off');

subplot(2,2,1);
hold on;
plot(t_lin(idx_min_lin), sig_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   sig_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   sig_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
yl = ylim;
plot([tF_lin tF_lin], yl, 'b:', 'HandleVisibility','off');
plot([tF_nl  tF_nl ], yl, 'r:', 'HandleVisibility','off');
plot([tF_am  tF_am ], yl, 'g:', 'HandleVisibility','off');
ylim(yl);
xlabel('Tiempo [s]'); ylabel('\sigma [Pa]');
title('Tension normal barra a=5');
legend('Location','northeast'); grid on;

subplot(2,2,2);
hold on;
plot(t_lin(idx_min_lin), y_nodo_b_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   y_nodo_b_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   y_nodo_b_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
yl = ylim;
plot([tF_lin tF_lin], yl, 'b:', 'HandleVisibility','off');
plot([tF_nl  tF_nl ], yl, 'r:', 'HandleVisibility','off');
plot([tF_am  tF_am ], yl, 'g:', 'HandleVisibility','off');
ylim(yl);
xlabel('Tiempo [s]'); ylabel('y_{nodo 4} [m]');
title('Coordenada actual nodo b=4');
legend('Location','northeast'); grid on;

subplot(2,2,3);
bar([tF_nl, tF_lin, tF_am]);
set(gca,'XTickLabel',{'a.i','a.ii','Amort.'});
ylabel('t_F [s]');
title('Tiempos finales');
grid on;

subplot(2,2,4);
hold on; axis equal; grid on;
idx_final_am = idx_min_am(end);
u_glob_am = zeros(nDOF,1);
u_glob_am(DOF_lib) = u_hist_am(idx_final_am,:)';
x_am = x0 + u_glob_am(1:2:end)';
y_am = y0 + u_glob_am(2:2:end)';
for e = 1:nBarras
    plot([x0(Ni(e)) x0(Nj(e))], [y0(Ni(e)) y0(Nj(e))], 'Color',[0.7 0.7 0.7], 'LineWidth',0.8);
end
for e = 1:nBarras
    plot([x_am(Ni(e)) x_am(Nj(e))], [y_am(Ni(e)) y_am(Nj(e))], 'g-', 'LineWidth',1.5);
end
plot(x_am, y_am, 'ko', 'MarkerFaceColor','g', 'MarkerSize',5);
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Reticulado con amortiguacion, t = %.2f s', t_am(idx_final_am)));

% =========================================================
%  ANIMACION - Configuracion deformada en el tiempo (a.i)
% =========================================================
fprintf('\nGenerando animacion...\n');

figure('Name','Animacion reticulado (no lineal)','NumberTitle','off');
idx_anim = find(t_nl <= tF_nl);
paso_anim = max(1, floor(length(idx_anim)/200));   % max 200 frames

for k = 1:paso_anim:length(idx_anim)
    clf;
    u_glob = zeros(nDOF,1);
    u_glob(DOF_lib) = u_hist_nl(idx_anim(k),:)';
    x_act = x0 + u_glob(1:2:end)';
    y_act = y0 + u_glob(2:2:end)';

    hold on; axis equal; grid on;
    % Configuracion inicial (gris)
    for e = 1:nBarras
        plot([x0(Ni(e)) x0(Nj(e))], [y0(Ni(e)) y0(Nj(e))], 'Color',[0.7 0.7 0.7], 'LineWidth',0.8);
    end
    % Configuracion deformada (azul)
    for e = 1:nBarras
        plot([x_act(Ni(e)) x_act(Nj(e))], [y_act(Ni(e)) y_act(Nj(e))], 'b-', 'LineWidth',1.5);
    end
    plot(x_act, y_act, 'ko', 'MarkerFaceColor','b', 'MarkerSize',5);
    title(sprintf('t = %.2f s  (no lineal)', t_nl(idx_anim(k))));
    xlabel('x [m]'); ylabel('y [m]');
    drawnow;
    pause(0.02);
end

fprintf('\n=== Simulacion completa. ===\n');

% =========================================================
%  FUNCIONES LOCALES
% =========================================================

% ---------------------------------------------------------
% FUNCION: propiedades geometricas de las barras
% dada configuracion actual de nodos (x, y)
% ---------------------------------------------------------
function [L_b, ct, st, kb] = props_barras(x, y, Ni, Nj, E, A)
    nB = length(Ni);
    L_b = zeros(1,nB); ct = zeros(1,nB); st = zeros(1,nB); kb = zeros(1,nB);
    for e = 1:nB
        dx    = x(Nj(e)) - x(Ni(e));
        dy    = y(Nj(e)) - y(Ni(e));
        L_b(e) = sqrt(dx^2 + dy^2);
        ct(e)  = dx / L_b(e);
        st(e)  = dy / L_b(e);
        kb(e)  = E * A / L_b(e);
    end
end

% ---------------------------------------------------------
% FUNCION: ensamblaje de K global
% ---------------------------------------------------------
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

% ---------------------------------------------------------
% FUNCION: area con signo de un triangulo
% ---------------------------------------------------------
function A_sig = area_triangulo(xi, yi, xj, yj, xk, yk)
    A_sig = 0.5 * ((xj-xi)*(yk-yi) - (xk-xi)*(yj-yi));
end

% ---------------------------------------------------------
% FUNCION: fuerza de arrastre en nodo c (solo direccion y)
% ---------------------------------------------------------
function Fa = drag_force(vy_c, y_c, h_SL, r_esf, rho_fl, Cd)
    % Posicion relativa de la esfera respecto a la superficie libre
    % y_c = coordenada actual del nodo c
    % La esfera esta centrada en y_c
    dist = h_SL - y_c;   % positivo si el centro esta bajo la superficie

    if dist <= -r_esf
        % Esfera completamente fuera del fluido (centro muy por encima)
        AR = 0;
    elseif dist >= r_esf
        % Esfera completamente sumergida
        AR = pi * r_esf^2;
    else
        % Parcialmente sumergida: AR = area del circulo a la altura del corte
        % h_inm = profundidad de inmersion del centro = dist
        AR = pi * (r_esf^2 - dist^2);
        AR = max(AR, 0);
    end

    % Fuerza de arrastre (solo componente y, opuesta a la velocidad)
    Fa = -0.5 * rho_fl * Cd * AR * abs(vy_c) * vy_c;
end

function dz = ode_lineal(~, z, K_red, F_red, M_inv, nLib)
    u = z(1:nLib);
    v = z(nLib+1:end);
    a = M_inv * (F_red - K_red * u);
    dz = [v; a];
end

function dz = ode_nolineal(~, z, x0, y0, Ni, Nj, E, A, P_val, ...
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

    % Recalcular K con geometria deformada
    [~, ct, st, kb] = props_barras(x_act, y_act, Ni, Nj, E, A);
    K_act = ensamblar_K(Ni, Nj, ct, st, kb, nNodos);
    K_red = K_act(DOF_lib, DOF_lib);

    % Fuerzas externas: cargas puntuales verticales.
    F_glob = zeros(nDOF,1);
    F_glob(2*7)  = -P_val;
    F_glob(2*10) = -P_val;

    % Amortiguador en nodo c (si corresponde)
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

function tF = detectar_tF(t_vec, u_hist, x0, y0, triangulos, DOF_lib, nNodos, tF_max)
    nDOF = 2*nNodos;
    nTri = size(triangulos,1);
    tF   = tF_max;

    % Areas iniciales
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

function [sigma, tau] = tension_barra(e, u_glob, x0, y0, Ni, Nj, E, L0b)
    ni = Ni(e); nj = Nj(e);

    % Desplazamientos actuales
    x_act_ni = x0(ni) + u_glob(2*ni-1);
    y_act_ni = y0(ni) + u_glob(2*ni);
    x_act_nj = x0(nj) + u_glob(2*nj-1);
    y_act_nj = y0(nj) + u_glob(2*nj);

    dx_act = x_act_nj - x_act_ni;
    dy_act = y_act_nj - y_act_ni;
    L_act  = sqrt(dx_act^2 + dy_act^2);

    % Tension normal (axial)
    sigma = E * (L_act - L0b(e)) / L0b(e);

    % Tension tangencial: CERO en barra articulada (solo transmite axial)
    tau = 0;
end

function [sigma_v, tau_v] = calcular_tensiones_hist(t_vec, u_hist, tF, barra_a, ...
                             x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0b)
    nDOF   = 2*nNodos;
    idx    = find(t_vec <= tF);
    sigma_v = zeros(length(idx),1);
    tau_v   = zeros(length(idx),1);
    for k = 1:length(idx)
        u_glob = zeros(nDOF,1);
        u_glob(DOF_lib) = u_hist(idx(k),:)';
        [sigma_v(k), tau_v(k)] = tension_barra(barra_a, u_glob, x0, y0, Ni, Nj, E, L0b);
    end
end

function [norm_max, t_norm_max] = norma_max_despl(t_vec, u_hist, DOF_lib, nNodos, tF)
    nDOF = 2*nNodos;
    idx  = find(t_vec <= tF);
    norms = zeros(length(idx),1);
    for k = 1:length(idx)
        u_glob = zeros(nDOF,1);
        u_glob(DOF_lib) = u_hist(idx(k),:)';
        % Norma del vector desplazamiento para cada nodo
        u_nodos = reshape(u_glob, 2, nNodos);
        normas_nodo = sqrt(u_nodos(1,:).^2 + u_nodos(2,:).^2);
        norms(k) = max(normas_nodo);
    end
    [norm_max, idx_max] = max(norms);
    t_norm_max = t_vec(idx(idx_max));
end
