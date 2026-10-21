clear; clc; close all;

% --- 1. PARÁMETROS DEL SISTEMA ---
m1 = 2.0;       % Masa del carrito (kg)
m2 = 1.0;       % Masa del péndulo (kg)
k  = 15.0;      % Constante del muelle (N/m)
l  = 1.5;       % Longitud del péndulo (m)
g  = 9.81;      % Gravedad (m/s^2)

% Condiciones iniciales
x0 = 1.0;       % Estiramiento inicial del muelle (m)
xp0 = 0.0;      % Velocidad inicial masa 1 (m/s)
theta0 = deg2rad(45); % Ángulo inicial (rad)
thetap0 = 0.0;  % Velocidad angular inicial (rad/s)

t0 = 0.0;
tf = 15.0;      % Tiempo total de simulación


% --- 2. INTEGRACIÓN NUMÉRICA CON ODE45 ---
tspan = [t0 tf];
y0 = [x0; xp0; theta0; thetap0];
options = odeset('RelTol', 1e-7, 'AbsTol', 1e-9);

% Llamada a ode45.
[t, Y] = ode45(@(t, y) ecuaciones(t, y, m1, m2, k, l, g), tspan, y0, options);

x_sol = Y(:,1);
theta_sol = Y(:,3);


% --- 3. GRÁFICAS DEL MOVIMIENTO ---
figure(1)
plot(t, x_sol, 'LineWidth', 1.5); hold on;
plot(t, theta_sol, 'LineWidth', 1.5);
grid on;
legend('Posición x (m)', 'Ángulo \theta (rad)');
xlabel('Tiempo (s)');
ylabel('Amplitudes');
title('Evolución temporal de las coordenadas');


% --- 4. ANIMACIÓN ---
figure(2)
set(gcf, 'Position', [100, 100, 800, 400]);
ndib = 300;
tdib = linspace(t0, tf, ndib);
xdib = interp1(t, x_sol, tdib);
thdib = interp1(t, theta_sol, tdib);

dt = tdib(2) - tdib(1);

for i = 1:ndib
    clf; hold on; grid on; axis equal;
    axis([-3 5 -2.5 1.5]); % Ajustar según parámetros
    
    % Coordenadas en este instante
    xc = xdib(i);
    yc = 0; % El carrito no se mueve en Y
    xp = xc + l * sin(thdib(i));
    yp = yc - l * cos(thdib(i));
    
    % Dibujar suelo y pared
    plot([-3 5], [-0.2 -0.2], 'k', 'LineWidth', 2); % Suelo
    plot([-3 -3], [-0.2 1.5], 'k', 'LineWidth', 2); % Pared
    
    % Dibujar muelle
    plot([-3 xc-0.5], [0 0], 'k--', 'LineWidth', 1.5);
    
    % Dibujar carrito (m1)
    rectangle('Position', [xc-0.5, yc-0.2, 1, 0.4], 'FaceColor', [0.8 0.8 0.8], 'Curvature', 0.2);
    
    % Dibujar péndulo (hilo y m2)
    plot([xc xp], [yc yp], 'b', 'LineWidth', 2);
    plot(xp, yp, 'ro', 'MarkerSize', 10, 'MarkerFaceColor', 'r');
    
    title(sprintf('Tiempo: %.2f s', tdib(i)));
    drawnow;
    
end

% --- FUNCIÓN DE LAS ECUACIONES DIFERENCIALES ---
function dydt = ecuaciones(~, y, m1, m2, k, l, g)
   
    dydt = zeros(4,1);
    
    x  = y(1);
    v  = y(2);
    th = y(3);
    w  = y(4);
    
    dydt(1) = v;
    dydt(3) = w;
    
    % Determinante de la matriz de masas para despejar aceleraciones
    delta = l * (m1 + m2 * sin(th)^2);
    
    % Aceleraciones despejadas
    termino_x = -k * x + m2 * l * w^2 * sin(th);
    termino_th = -g * sin(th);
    
    dydt(2) = (l * termino_x - m2 * l * cos(th) * termino_th) / delta;
    dydt(4) = ((m1 + m2) * termino_th - cos(th) * termino_x) / delta;
end