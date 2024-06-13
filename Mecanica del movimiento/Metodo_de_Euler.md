# Métodos Numéricos para la Simulación de Movimiento

En la simulación de movimientos, a menudo es necesario resolver ecuaciones diferenciales ordinarias (EDOs) que describen el comportamiento dinámico de un sistema. Estas ecuaciones pueden modelar una amplia variedad de fenómenos físicos, desde el movimiento de planetas hasta la dinámica de fluidos. Sin embargo, en muchos casos, no es posible encontrar soluciones analíticas exactas para estas ecuaciones, y se requieren métodos numéricos para obtener aproximaciones.

Los métodos numéricos son técnicas que permiten aproximar la solución de EDOs mediante cálculos discretos y repetitivos. A lo largo de los años, se han desarrollado numerosos métodos para este propósito, cada uno con sus propias características, ventajas y limitaciones. Algunos de los métodos más conocidos incluyen el Método de Euler, el Método de Runge-Kutta, y las Series de Taylor, entre otros.

En este primer contacto entenderemos el metodo mas sencillo para resolver un EDOs. ya sea por su simplicidad, facilidad de implementacion y por ser la base para metodos mas avanzados.

## Método de Euler

El Método de Euler es un enfoque numérico para resolver EDOs de primer orden mediante aproximaciones por pasos. A pesar de su simplicidad, este metodo puede ser inexacto y sufre de problemas de estabilidad, especialmente con pasos de tiempo grandes. 

### Concepto Teorico

En el Método de Euler, la derivada de una función en un punto se aproxima mediante la tangente a la curva en ese punto. Esta aproximación se utiliza para predecir el valor de la función en el siguiente punto de tiempo.

Para entenderlo de una manera menos técnica, imagina que estás en un campo de montañas en una noche muy oscura con solo una pequeña linterna. Tu objetivo es llegar al punto más bajo del campo (el valle). La linterna solo te permite ver un pequeño espacio a tu alrededor.

El método de Euler es como decidir en qué dirección caminar basándote solo en la pendiente del terreno que puedes ver con tu linterna. Si el terreno desciende hacia la derecha, decides dar un paso a la derecha. Si desciende hacia la izquierda, das un paso a la izquierda.

No puedes ver el valle desde donde estás, y este método no garantiza que siempre tomes la ruta más rápida o eficiente. Pero si sigues tomando pasos basados en la pendiente inmediata, eventualmente llegarás al valle.

En esta analogía, el valle es la solución de la ecuación diferencial, la pendiente del terreno es la derivada, y la dirección que eliges para caminar es la aproximación que hace el método de Euler.

### Ejemplo Sencillo: Ilustrando el enfriamiento de un Objeto

El enfriamiento de un objeto según la ley de Newton. La ecuación diferencial que describe este proceso es:

$$
\frac{dT}{dt} = -k(T - T_{\text{amb}})
$$

donde:

- $T$ es la temperatura del objeto en el tiempo $t$,
- $T_{\text{amb}}$ es la temperatura ambiente,
- $k$ es una constante de proporcionalidad.

Supongamos que queremos predecir cómo cambia la temperatura del objeto con el tiempo. Utilizando el Método de Euler, podemos aproximar la temperatura en pasos discretos de tiempo $\Delta t$. La fórmula de Euler para este problema sería:

$$
T_{n+1} = T_n + \Delta t(-k(T_n - T_{\text{amb}}))
$$

donde $T_n$ es la temperatura en el tiempo $t_n$ y $T_{n+1}$ es la temperatura en el siguiente paso de tiempo $t_{n+1} = t_n + \Delta t$.

###  ¿Como podemos usar este metodo en MATLAB?

Para demostrar cómo se puede usar el Método de Euler en este software, te mostare paso a paso cómo se realiza la simulación de un péndulo cónico. Este ejemplo es sencillo pero preciso para explicar qué está pasando, incluyendo la declaración de variables y cómo se afecta el radio del péndulo según su velocidad.

#### Para empezar definimos las constantes y condiciones iniciales para la simulación: 

```matlab
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
```

#### Luego inicializamos los vectores que almacenaran todos los resultados:

```matlab
t = 0:dt:T;                  % este sera nuestro vector tiempo 
Theta = zeros(1, length(t)); % este vector almacenara el angulo de oscilacion
Phi = zeros(1, length(t));   % este vector almacenara el angulo de rotacion
Z1 = zeros(1, length(t));    % este vector almacenara la derivada de theta 
Z2 = zeros(1, length(t));    % este vector almacenar la derivada de phi
R = zeros(1, length(t));     % este vector almacenara el radio
```
La función **zeros** en MATLAB se utiliza para crear una matriz de ceros. En esta parte del código, se utiliza para inicializar los vectores Theta, Phi, Z1 y Z2 antes de entrar en el bucle de simulación. para explicarlo mejor podriamos decir que zeros(1, length(t)) crea una matriz de ceros con 1 fila y un número de columnas igual a la longitud del vector de tiempo t. Esto se hace para preasignar memoria para estos vectores, lo que puede hacer que el código se ejecute más rápido.

#### Luego simplemente asignaremos las condiciones iniciales:

```matlab
Theta(1) = theta; 
Phi(1) = phi;
Z1(1) = omega_theta;
Z2(1) = omega_phi;
R(1) = constante * omega_phi;
```

#### Ahora se creara un bucle y se utilizara la ecuacion de euler tomando claramente en cuenta su formula: y(t + dt) = y(t) + dt * f(t, y(t))

```matlab
for i = 1:(length(t) - 1)
    Z1(i+1) = Z1(i) - (g/L)*sin(Theta(i))*dt;
    Z2(i+1) = Z2(i) - (g/L)*Theta(i)*sin(Phi(i))*dt;
    Theta(i+1) = Theta(i) + Z1(i+1)*dt;
    Phi(i+1) = Phi(i) + Z2(i+1)*dt;
    R(i+1) = constante * Z2(i+1); 
end
```

#### Y por ultimo haremos la graficación del mismo 

```matlab
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
```

En esta demostración, hemos utilizado el Método de Euler para simular el movimiento de un péndulo cónico. Este enfoque nos ha permitido aproximar la solución de las ecuaciones diferenciales que describen el sistema de manera sencilla y efectiva. A través de la declaración de variables y la implementación de los cálculos paso a paso, hemos observado cómo la velocidad angular del péndulo afecta el radio de su trayectoria. Aunque el Método de Euler es una técnica básica, su simplicidad lo convierte en una herramienta poderosa para comprender los fundamentos de la simulación numérica.

En futuras secciones, exploraremos métodos numéricos más avanzados que ofrecen mayor precisión y estabilidad, como el Método de Runge-Kutta. Sin embargo, el conocimiento adquirido a través del Método de Euler proporciona una base sólida que facilitará el entendimiento de técnicas más complejas. Esta experiencia práctica con el péndulo cónico demuestra cómo los métodos numéricos pueden ser aplicados a problemas reales, permitiendo el análisis y la predicción del comportamiento dinámico de sistemas físicos.




