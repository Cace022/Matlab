# Métodos Numéricos para la Simulación de Movimiento

En la simulación de movimientos, a menudo es necesario resolver ecuaciones diferenciales ordinarias (EDOs) que describen el comportamiento dinámico de un sistema. Sin embargo, en muchos casos, no es posible encontrar soluciones analíticas exactas, y se requieren métodos numéricos para obtener aproximaciones. En esta sección, exploraremos dos métodos numéricos comúnmente utilizados: el Método de Euler y el Método de Lagrange.

## Método de Euler

El Método de Euler es un enfoque numérico para resolver ecuaciones diferenciales ordinarias (EDOs) de primer orden mediante aproximaciones por pasos. Es simple de implementar pero puede tener precisión limitada en comparación con métodos más avanzados.

### Concepto

En el Método de Euler, la derivada de una función en un punto se aproxima mediante la tangente a la curva en ese punto. Esta aproximación se utiliza para predecir el valor de la función en el siguiente punto de tiempo.

### Formula

<p align="center">
  <img src="https://i.ytimg.com/vi/Px7vI36WKqg/maxresdefault.jpg" alt="imagen">
</p>

### Ejemplo en MATLAB

```matlab
% Parámetros del Método de Euler
h = 0.1;         % Tamaño del paso de tiempo
t = 0:h:10;      % Vector de tiempo de 0 a 10 con paso h
y0 = 1;          % Condición inicial

% Implementación del Método de Euler
y = zeros(size(t));  % Vector para almacenar los valores de y
y(1) = y0;           % Asignar la condición inicial
for i = 1:length(t)-1
   y(i+1) = y(i) + h * f(t(i), y(i));  % Método de Euler
end

% Función que define la EDO dy/dt = f(t, y)
function dydt = f(t, y)
   dydt = -0.2 * y;  % Ejemplo de ecuación diferencial: y' = -0.2*y
end

% Graficar resultados
figure;
plot(t, y, '-o');
xlabel('Tiempo');
ylabel('y(t)');
title('Método de Euler para EDO');
grid on;