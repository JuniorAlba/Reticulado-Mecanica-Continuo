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
% GEOMETRIA Y CONECTIVIDAD DESDE ARCHIVO DXF
% Lectura automatica con f_LectDxf (provisto por la catedra)
% ---------------------------------------------------------
dxf_path = fullfile(fileparts(mfilename('fullpath')), ...
    '..', 'Consigna', 'dxf', 'dxf', 'Estr12.dxf');
[c_Line, ~, c_Cir, ~, ~] = f_LectDxf(dxf_path);

% Nodos: centros de los circulos (en el orden en que aparecen en el DXF)
m_Cir = cell2mat(c_Cir(:,1));    % [x, y, radio] para cada circulo
x0    = round(m_Cir(:,1)' * 1e6) / 1e6;   % coordenadas x iniciales (1 x nNodos)
y0    = round(m_Cir(:,2)' * 1e6) / 1e6;   % coordenadas y iniciales (1 x nNodos)

nNodos = numel(x0);

% Barras: lineas cuyos dos extremos coincidan con nodos conocidos
% (las lineas de cota/dimension del DXF tienen extremos fuera de la red)
m_Lin = cell2mat(c_Line(:,1));   % [Xi Yi Zi Xj Yj Zj]
tol   = 0.5;                     % tolerancia geometrica [m]
Ni_v = []; Nj_v = [];
for e_dxf = 1:size(m_Lin, 1)
    xi = m_Lin(e_dxf,1);  yi = m_Lin(e_dxf,2);
    xj = m_Lin(e_dxf,4);  yj = m_Lin(e_dxf,5);
    di = sqrt((x0 - xi).^2 + (y0 - yi).^2);
    dj = sqrt((x0 - xj).^2 + (y0 - yj).^2);
    [min_di, ni] = min(di);
    [min_dj, nj] = min(dj);
    if min_di < tol && min_dj < tol && ni ~= nj
        Ni_v(end+1) = ni;  %#ok<AGROW>
        Nj_v(end+1) = nj;  %#ok<AGROW>
    end
end
Ni = Ni_v;
Nj = Nj_v;
nBarras = numel(Ni);

fprintf('Geometria leida desde DXF: %d nodos, %d barras\n', nNodos, nBarras);

if numel(y0) ~= nNodos
    error('x0 e y0 deben tener la misma cantidad de nodos.');
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
tF_min = min(tF_lin, tF_nl);   % intervalo comun b.iii

[sig_lin, tau_lin] = calcular_tensiones_hist(t_lin, u_hist_lin, tF_lin, barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);
[sig_nl,  tau_nl]  = calcular_tensiones_hist(t_nl,  u_hist_nl,  tF_nl,  barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);
[sig_am,  tau_am]  = calcular_tensiones_hist(t_am,  u_hist_am,  tF_am,  barra_a, ...
                      x0, y0, Ni, Nj, DOF_lib, nNodos, E, L0);

% Indices hasta cada tF propio (para graficas b.ii)
idx_min_lin = find(t_lin <= tF_lin);
idx_min_nl  = find(t_nl  <= tF_nl);
idx_min_am  = find(t_am  <= tF_am);

% Indices recortados al intervalo comun [0, tF_min] (para comparacion b.iii)
idx_com_lin = find(t_lin <= tF_min);
idx_com_nl  = find(t_nl  <= tF_min);

% --- Coordenada actual del nodo b ---
DOF_b_x = find(DOF_lib == 2*nodo_b-1);  %#ok<NASGU>
DOF_b_y = find(DOF_lib == 2*nodo_b);

y_nodo_b_lin = y0(nodo_b) + u_hist_lin(idx_min_lin, DOF_b_y);
y_nodo_b_nl  = y0(nodo_b) + u_hist_nl( idx_min_nl,  DOF_b_y);
y_nodo_b_am  = y0(nodo_b) + u_hist_am( idx_min_am,  DOF_b_y);

% Coordenada nodo b en intervalo comun (para b.iii)
y_nodo_b_lin_com = y0(nodo_b) + u_hist_lin(idx_com_lin, DOF_b_y);
y_nodo_b_nl_com  = y0(nodo_b) + u_hist_nl( idx_com_nl,  DOF_b_y);

% --- Norma del vector desplazamiento en el tiempo (para b.iii) ---
[nm_lin, t_nm_lin] = norma_max_despl(t_lin, u_hist_lin, DOF_lib, nNodos, tF_lin);
[nm_nl,  t_nm_nl]  = norma_max_despl(t_nl,  u_hist_nl,  DOF_lib, nNodos, tF_nl);
[nm_am,  t_nm_am]  = norma_max_despl(t_am,  u_hist_am,  DOF_lib, nNodos, tF_am);

% Evolucion de la norma maxima en el intervalo comun (para grafica b.iii)
nDOF = 2*nNodos;
norm_ev_lin = zeros(length(idx_com_lin),1);
for k = 1:length(idx_com_lin)
    u_g = zeros(nDOF,1);
    u_g(DOF_lib) = u_hist_lin(idx_com_lin(k),:)';
    un = reshape(u_g,2,nNodos);
    norm_ev_lin(k) = max(sqrt(un(1,:).^2 + un(2,:).^2));
end
norm_ev_nl = zeros(length(idx_com_nl),1);
for k = 1:length(idx_com_nl)
    u_g = zeros(nDOF,1);
    u_g(DOF_lib) = u_hist_nl(idx_com_nl(k),:)';
    un = reshape(u_g,2,nNodos);
    norm_ev_nl(k) = max(sqrt(un(1,:).^2 + un(2,:).^2));
end

fprintf('\n--- Resumen de resultados ---\n');
fprintf('                    Lineal (a.ii)   No lineal (a.i)   Con amortiguador\n');
fprintf('tF [s]              %10.4f      %10.4f        %10.4f\n', tF_lin, tF_nl, tF_am);
fprintf('min(tFa.i, tFa.ii)  %10.4f\n', tF_min);
fprintf('Norma max despl     %10.4f      %10.4f        %10.4f\n', nm_lin, nm_nl, nm_am);
fprintf('Tiempo norma max    %10.4f      %10.4f        %10.4f\n', t_nm_lin, t_nm_nl, t_nm_am);

% =========================================================
%  GRAFICAS
% =========================================================

% =========================================================
%  FIGURA 1: Tension normal barra a (b.ii) - cada caso hasta su propio tF
% =========================================================
figure('Name','Tension normal barra a=5','NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), sig_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   sig_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   sig_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
xlabel('Tiempo [s]'); ylabel('\sigma [Pa]');
title(sprintf('Tension normal en barra a=%d', barra_a));
legend('Location','northeast'); grid on;

% =========================================================
%  FIGURA 2: Tension tangencial barra a (b.ii)
%  tau = 0 siempre en reticulado articulado (sin momento flector)
% =========================================================
figure('Name','Tension tangencial barra a=5','NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), tau_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   tau_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   tau_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
xlabel('Tiempo [s]'); ylabel('\tau [Pa]');
title(sprintf('Tension tangencial en barra a=%d  (\tau = 0: barra articulada)', barra_a));
legend('Location','northeast'); grid on;
ytext = sprintf(['\tau = 0 en todo momento porque cada barra solo puede transmitir\n' ...
    'fuerza axial (union articulada, sin momento flector ni corte).']);
text(0.02, 0.5, ytext, 'Units','normalized', 'FontSize', 9, ...
    'VerticalAlignment','middle', 'BackgroundColor',[1 1 0.85], 'EdgeColor','k');

% =========================================================
%  FIGURA 3: Coordenada actual nodo b (b.ii)
% =========================================================
figure('Name',sprintf('Coordenada actual nodo b=%d',nodo_b),'NumberTitle','off');
hold on;
plot(t_lin(idx_min_lin), y_nodo_b_lin, 'b-',  'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_min_nl),   y_nodo_b_nl,  'r-',  'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
plot(t_am(idx_min_am),   y_nodo_b_am,  'g--', 'LineWidth', 1.5, 'DisplayName','Con amortiguador');
xlabel('Tiempo [s]'); ylabel(sprintf('y_{nodo %d} [m]', nodo_b));
title(sprintf('Coordenada actual del nodo b=%d', nodo_b));
legend('Location','northeast'); grid on;

% =========================================================
%  FIGURA 4 (b.iii): Comparacion a.i vs a.ii en intervalo comun [0, tF_min]
% =========================================================
figure('Name', sprintf('b.iii - Comparacion en [0, min(tFa.i, tFa.ii)] = [0, %.2f s]', tF_min), ...
    'NumberTitle','off');

subplot(2,2,1);   % --- Tension normal en intervalo comun ---
hold on;
plot(t_lin(idx_com_lin), sig_lin(1:length(idx_com_lin)), 'b-', 'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_com_nl),   sig_nl( 1:length(idx_com_nl)),  'r-', 'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
xlabel('Tiempo [s]'); ylabel('\sigma [Pa]');
title(sprintf('\\sigma barra a=%d  en [0, tF_{min}=%.2f s]', barra_a, tF_min));
legend('Location','northeast'); grid on;

subplot(2,2,2);   % --- Coordenada actual nodo b en intervalo comun ---
hold on;
plot(t_lin(idx_com_lin), y_nodo_b_lin_com, 'b-', 'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_com_nl),   y_nodo_b_nl_com,  'r-', 'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
xlabel('Tiempo [s]'); ylabel(sprintf('y_{nodo %d} [m]', nodo_b));
title(sprintf('Coord. actual nodo b=%d  en [0, tF_{min}=%.2f s]', nodo_b, tF_min));
legend('Location','northeast'); grid on;

subplot(2,2,3);   % --- Norma maxima de desplazamiento en intervalo comun ---
hold on;
plot(t_lin(idx_com_lin), norm_ev_lin, 'b-', 'LineWidth', 1.5, 'DisplayName','Lineal (a.ii)');
plot(t_nl(idx_com_nl),   norm_ev_nl,  'r-', 'LineWidth', 1.5, 'DisplayName','No lineal (a.i)');
% Marcar el maximo de cada caso dentro del intervalo comun
[~, ki] = max(norm_ev_lin); plot(t_lin(idx_com_lin(ki)), norm_ev_lin(ki), 'bs', 'MarkerSize',8, 'MarkerFaceColor','b', 'HandleVisibility','off');
[~, ki] = max(norm_ev_nl);  plot(t_nl( idx_com_nl(ki)),  norm_ev_nl(ki),  'r^', 'MarkerSize',8, 'MarkerFaceColor','r', 'HandleVisibility','off');
xlabel('Tiempo [s]'); ylabel('max_I ||u_I|| [m]');
title(sprintf('Norma max. despl. nodal  en [0, tF_{min}=%.2f s]', tF_min));
legend('Location','northwest'); grid on;

subplot(2,2,4);   % --- Tiempos finales y deformada al tF_min ---
bar([tF_nl, tF_lin, tF_am], 0.5);
set(gca,'XTickLabel',{'a.i (NL)','a.ii (Lin)','Amort.'});
ylabel('t_F [s]');
title(sprintf('Tiempos de colapso  |  t_{F,min} = %.2f s', tF_min));
text(1, tF_nl  * 0.5, sprintf('%.2f s',tF_nl),  'HorizontalAlignment','center','Color','w','FontWeight','bold');
text(2, tF_lin * 0.5, sprintf('%.2f s',tF_lin), 'HorizontalAlignment','center','Color','w','FontWeight','bold');
text(3, tF_am  * 0.5, sprintf('%.2f s',tF_am),  'HorizontalAlignment','center','Color','w','FontWeight','bold');
grid on;

% sgtitle no existe en Octave; se usa un axes oculto centrado como titulo de figura
ax_title = axes('Position',[0 0.96 1 0.04], 'Visible','off');
text(ax_title, 0.5, 0.5, ...
    sprintf('b.iii - Comparacion a.i vs a.ii en [0, min(tF)] = [0, %.2f s]  |  Marcadores: norma max. despl.', tF_min), ...
    'HorizontalAlignment','center', 'VerticalAlignment','middle', ...
    'FontSize', 10, 'Units','normalized');

% =========================================================
%  FIGURA 5: Vista conjunta con deformada (referencia geometrica)
% =========================================================
figure('Name','Comparacion a.i, a.ii y reticulado amortiguado - deformada','NumberTitle','off');

subplot(1,2,1);
hold on; axis equal; grid on;
idx_final_nl = idx_min_nl(end);
u_glob_nl_f = zeros(nDOF,1);
u_glob_nl_f(DOF_lib) = u_hist_nl(idx_final_nl,:)';
x_nl_f = x0 + u_glob_nl_f(1:2:end)';
y_nl_f = y0 + u_glob_nl_f(2:2:end)';
for e = 1:nBarras
    plot([x0(Ni(e)) x0(Nj(e))], [y0(Ni(e)) y0(Nj(e))], 'Color',[0.75 0.75 0.75], 'LineWidth',0.8);
end
for e = 1:nBarras
    plot([x_nl_f(Ni(e)) x_nl_f(Nj(e))], [y_nl_f(Ni(e)) y_nl_f(Nj(e))], 'r-', 'LineWidth',1.5);
end
plot(x_nl_f, y_nl_f, 'ko', 'MarkerFaceColor','r', 'MarkerSize',5);
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('No lineal (a.i), t = %.2f s', t_nl(idx_final_nl)));

subplot(1,2,2);
hold on; axis equal; grid on;
idx_final_am = idx_min_am(end);
u_glob_am = zeros(nDOF,1);
u_glob_am(DOF_lib) = u_hist_am(idx_final_am,:)';
x_am = x0 + u_glob_am(1:2:end)';
y_am = y0 + u_glob_am(2:2:end)';
for e = 1:nBarras
    plot([x0(Ni(e)) x0(Nj(e))], [y0(Ni(e)) y0(Nj(e))], 'Color',[0.75 0.75 0.75], 'LineWidth',0.8);
end
for e = 1:nBarras
    plot([x_am(Ni(e)) x_am(Nj(e))], [y_am(Ni(e)) y_am(Nj(e))], 'g-', 'LineWidth',1.5);
end
plot(x_am, y_am, 'ko', 'MarkerFaceColor','g', 'MarkerSize',5);
xlabel('x [m]'); ylabel('y [m]');
title(sprintf('Con amortiguador, t = %.2f s', t_am(idx_final_am)));

% =========================================================
%  ANIMACION - Configuracion deformada en el tiempo (a.i)
%  Se guarda como GIF usando gif.m (provisto por la catedra)
% =========================================================
fprintf('\nGenerando animacion GIF (hipotesis no lineal)...\n');

gif_path = fullfile(fileparts(mfilename('fullpath')), 'animacion_caso12.gif');

fig_anim = figure('Name','Animacion reticulado (no lineal)','NumberTitle','off');
idx_anim = find(t_nl <= tF_nl);
paso_anim = max(1, floor(length(idx_anim)/120));   % max 120 frames en el GIF

% Limites fijos del grafico para toda la animacion
u_glob_all = zeros(nDOF, length(idx_anim));
for k = 1:length(idx_anim)
    u_glob_all(DOF_lib, k) = u_hist_nl(idx_anim(k), :)';
end
x_all = x0' + u_glob_all(1:2:end, :);  % 15 x nFrames
y_all = y0' + u_glob_all(2:2:end, :);
xlims = [min(x_all(:))-2, max(x_all(:))+2];
ylims = [min(y_all(:))-2, max(y_all(:))+2];

first_gif = true;
for k = 1:paso_anim:length(idx_anim)
    u_glob = zeros(nDOF,1);
    u_glob(DOF_lib) = u_hist_nl(idx_anim(k),:)';
    x_act = x0 + u_glob(1:2:end)';
    y_act = y0 + u_glob(2:2:end)';

    clf;
    hold on; axis equal; grid on;
    xlim(xlims); ylim(ylims);

    % Configuracion inicial (gris)
    for e = 1:nBarras
        plot([x0(Ni(e)) x0(Nj(e))], [y0(Ni(e)) y0(Nj(e))], ...
            'Color',[0.7 0.7 0.7], 'LineWidth', 0.8);
    end
    % Configuracion deformada (azul)
    for e = 1:nBarras
        plot([x_act(Ni(e)) x_act(Nj(e))], [y_act(Ni(e)) y_act(Nj(e))], ...
            'b-', 'LineWidth', 1.5);
    end
    plot(x_act, y_act, 'ko', 'MarkerFaceColor','b', 'MarkerSize', 5);
    title(sprintf('t = %.2f s  |  No lineal (a.i)  |  Caso 12', t_nl(idx_anim(k))), ...
        'FontSize', 10);
    xlabel('x [m]'); ylabel('y [m]');
    drawnow;

    % Guardar frame en GIF
    if first_gif
        gif(gif_path, 'overwrite', true, 'frame', gca, ...
            'DelayTime', 1/15, 'LoopCount', 1);
        first_gif = false;
    else
        gif;
    end
end
gif('clear');
fprintf('GIF guardado en: %s\n', gif_path);

fprintf('\n=== Simulacion completa. ===\n');

