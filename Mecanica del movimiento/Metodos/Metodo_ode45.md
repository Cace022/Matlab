# Metodo Ode45

El metodo ode45 es un solver de MATLAB que es utilizado para resolver ecuaciones diferenciales ordinarias (EDO). Este metodo se basa en el metodo de Runge-Kutta de 4 y 5 orden, que cabe destacar es un metodo avanzado para resolver ecuaciones diferenciales ordinarias. 

## ¿Que es el metodo de Runge-Kutta?

Este metodo no es mas que una generalización del método de Euler para resolver EDO. Mientras que el método de Euler aproxima la solución de la EDO avanzando un pequeño paso en la dirección de la derivada, el método de Runge-Kutta toma varios pasos intermedios para obtener una mejor aproximación de la derivada.

### Ahora como lo podriamos explicar en terminos muy simples 

Bueno imagina que tienes una bola rodando cuesta abajo y quieres saber dónde estará después de cierto tiempo. La ecuación diferencial te dice cómo cambia la posición de la bola con el tiempo. 

### Ejemplo

Imaginemos un problema bastante sencillo: tienes una planta que crece a una tasa muy dependiente de su estado actual, la ecuacion diferencial podria ser: 

$$\frac{dy}{dt} = k \cdot y$$

Aquí, **𝑦** es el tamaño de la planta, **𝑡** es el tiempo, y **𝑘** es una constante que describe la rapidez con la que crece la planta.

```matlab
% Definimos la función que representa la EDO

k = 0.1;
f = @(t, y) k * y;

% Establecemos los parametros y condiciones iniciales

t0 = 0; % Tiempo inicial
tf = 10; % Tiempo final
h = 1; % Tamaño del paso
y0 = 1; % Condición inicial


% Inicialización
t = t0:h:tf; % Vector de tiempo
y = zeros(size(t)); % Vector de solución
y(1) = y0; % Condición inicial

% Método de Runge-Kutta de cuarto orden
for n = 1:(length(t)-1)
    k1 = f(t(n), y(n));
    k2 = f(t(n) + h/2, y(n) + h/2 * k1);
    k3 = f(t(n) + h/2, y(n) + h/2 * k2);
    k4 = f(t(n) + h, y(n) + h * k3);
    
    y(n+1) = y(n) + h/6 * (k1 + 2*k2 + 2*k3 + k4);
end

% Graficar los resultados
plot(t, y);
xlabel('Tiempo t');
ylabel('Tamaño de la planta y');
title('Crecimiento de la planta usando el método de Runge-Kutta');
grid on;
```

Este metodo tiende a ser un poco mas complicado que el de Euler por lo cual se recomienda una practica mas extensa y constante.
