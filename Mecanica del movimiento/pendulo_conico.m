% Definicion de las constantes y de las condiciones iniciales 

m = 1;               % masa del pendulo 
L = 1;               % longitud del pendulo 
g = 9.81;            % aceleracion debido a la gravedad 
theta = pi/4;        % este sera mi angulo inicial de oscilacion 
phi = 0;             % este sera mi angulo inicial de rotacion 
omega_theta = 0;     % esta sera mi velocidad angular inicial de oscilacion
omega_phi = 1;       % esta sera mi velocidad angular inicial de rotacion
T = 10;              % este sera el tiempo total de mi simulacion
dt = 0.1;            % este sera el paso del tiempo 
constante = 1;       % esta sera la constante de proporcionalidad para el radio

% Inicializacion de vectores (para almacenar los resultados) 

t = 0:dt:T;                  % este sera nuestro vector tiempo 
Theta = zeros(1, length(t)); % este vector almacenara el angulo de oscilacion
Phi = zeros(1, length(t));   % este vector almacenara el angulo de rotacion
Z1 = zeros(1, length(t));    % este vector almacenara la derivada de theta 
Z2 = zeros(1, length(t));    % este vector almacenar la derivada de phi
R = zeros(1, length(t));     % este vector almacenara el radio

% Ahora asignare las variables iniciales 

Theta(1) = theta; 
Phi(1) = phi;
Z1(1) = omega_theta;
Z2(1) = omega_phi;
R(1) = constante * omega_phi;

% Ahora creare un bucle para la simulacion 
% para esto utilizare el metodo de Euler para la integracion numerica 
% tomando en cuenta que la formula de euler es: y(t + dt) = y(t) + dt * f(t, y(t))

for i = 1:(length(t) - 1)
    Z1(i+1) = Z1(i) - (g/L)*sin(Theta(i))*dt;
    Z2(i+1) = Z2(i) - (g/L)*Theta(i)*sin(Phi(i))*dt;
    Theta(i+1) = Theta(i) + Z1(i+1)*dt;
    Phi(i+1) = Phi(i) + Z2(i+1)*dt;
    R(i+1) = constante * Z2(i+1); 
end

% aqui esta la grafiacion muy basica del pendulo 
figure;                                % esto creara una figura donde se dibujaran los subgraficos
subplot(3,1,1);                        % esto dividira la figura en 3 filas y una columna 
plot(t, Theta);                        % esto dibuja el angulo de oscilacion theta en funcion al tiempo
xlabel('Tiempo (s)');                  % esto establece la etiqueta del eje x del primer subgrafico
ylabel('Ángulo de oscilación (rad)');  % esto establece las etiqueta del eje y del primer subgrafico
title('Péndulo cónico');               % esto solo pone el titulo del primer subagrafico

subplot(3,1,2);                        % aqui estoy seleccionando el segundo sugrafico 
plot(t, Phi);                          % esto dibuja el angulo de rotacion phi en funcion al tiempo
xlabel('Tiempo (s)');                  % esto es la etiqueta del eje x 
ylabel('Ángulo de rotación (rad)');    % esto es la etiqueta del eje y 

subplot(3,1,3);                        % aqui estoy seleccionando el tercer sugrafico 
plot(t, R);                            % aqui se dibuja eñ radio el funcion al tiempo 
xlabel('Tiempo (s)');                  % esto es la etiqueta del eje x 
ylabel('Radio (m)');                   % esto es la etiqueta del eje y