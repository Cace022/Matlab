# Cinematica

## Antes de comenzar
Antes de entrar en los conceptos específicos del movimiento, es importante comprender su relevancia fundamental en nuestras vidas. El movimiento es una parte integral de nuestro entorno y experiencia diaria, aunque a veces pase desapercibido.

Desde el simple acto de caminar hasta los complejos movimientos de vehículos y maquinaria, todo a nuestro alrededor se encuentra en constante movimiento. Incluso los objetos aparentemente estáticos, como un edificio o una roca, experimentan movimientos imperceptibles a simple vista, como vibraciones o desplazamientos microscópicos.

**El movimiento es una manifestación fundamental de la naturaleza**, y su comprensión nos permite descifrar los patrones y leyes que rigen el universo. A través de la física y las matemáticas, se han desarrollado modelos y herramientas para representar, predecir y controlar el movimiento de objetos en diversas situaciones.

Comprender el movimiento no solo tiene implicaciones teóricas, sino también aplicaciones prácticas en campos como la ingeniería, la robótica, la astronomía y muchos más. Desde el diseño de vehículos eficientes hasta el lanzamiento de naves espaciales, el conocimiento del movimiento es esencial para el avance tecnológico y científico.

## ¿Qué es el movimiento?
Cuando piensas en movimiento, seguramente crees que solamente es la acción de moverse, pero científicamente podríamos definirlo como el cambio de posición de un cuerpo u objeto con respecto a un sistema de referencia y un intervalo de tiempo determinado. Este se describe utilizando conceptos como posición, velocidad, aceleración y trayectoria, expresados mediante ecuaciones matemáticas derivadas de las leyes de la física.

### Conceptos clave
- **Posición**: La ubicación de un objeto en relación con un sistema de coordenadas específico.
- **Velocidad**: La rapidez con la que un objeto cambia su posición con respecto al tiempo.
- **Aceleración**: El cambio de velocidad de un objeto por unidad de tiempo.
- **Trayectoria**: La línea o camino que describe el movimiento de un objeto.
- **Sistema de referencia**: Un marco de referencia desde el cual se mide el movimiento de un objeto.

### Tipos de movimiento
Como conceptos básicos, tenemos el Movimiento Rectilíneo Uniforme (MRU) y el Movimiento Rectilíneo Uniformemente Acelerado (MRA). Además, existen otros tipos de movimiento, como el movimiento circular, el movimiento parabólico, entre otros, pero hoy nos enfocaremos en los 2 mas conocidos.

#### Movimiento Rectilíneo Uniforme (MRU)
El MRU se refiere a un movimiento en línea recta con una velocidad constante. En otras palabras, un objeto que se mueve en MRU no experimenta cambios en su velocidad ni en su dirección. Algunas características importantes del MRU son:

    - La velocidad es constante en todo momento.
    - La aceleración es cero.
    - La trayectoria es una línea recta.

### Ejemplo en MATLAB
```matlab
% Parametros del MRU
x0 = 0;       % Esta es nuestra posicion inicial (m).
v = 5;        % Esta es nuestra velocidad constante (m/s).
t = 0:0.1:10; % Este es un vector tiempo que va desde 0 a 10 con incrementos de 0.1 segundos.

% Establecemos la ecuación de movimiento 
x = x0 + v * t; 

% Para hacer su grafica utilizamos la siguiente sintaxis
figure; % Esto crea una figura nueva para graficar.
plot(t, x);   % Grafica la posicion 'x' en funcion al tiempo 't'.
xlabel('Tiempo (s)');   % Se etiqueta un eje horizontal como tiempo. 
ylabel('posicion (m)'); % Se etiqueta un eje vertical como posicion.
title('Movimiento Rectilinio Uniforme'); 
grid on;                % muestra una cuadricula en la grafica.
```

### Movimiento Rectilíneo Uniformemente Acelerado (MRA)

Por otro lado, el MRA describe un movimiento en línea recta donde la velocidad cambia constantemente debido a una aceleración constante. En este caso, el objeto experimenta una aceleración constante, lo que significa que su velocidad aumenta o disminuye a un ritmo constante. Las características clave del MRA son:

    - La aceleración es constante y distinta de cero.
    - La velocidad cambia de manera uniforme (aumenta o disminuye constantemente).
    - La trayectoria es una línea recta.

### Ejemplo en MATLAB

```matlab
% Parametros del MRA
x0 = 0; 
v0 = 0; 
a = 2; 
t = 0:0.1:10;

% Establecemos la ecuacion de movimiento 
x = x0 + v0 * t + 0.5 * a * t.^2;

% Para hacer su grafica utilizamos la siguiente sintaxis
figure; 
plot (t, x);   
xlabel('Tiempo (s)');   
ylabel('posicion (m)'); 
title('Movimiento Rectilinio Uniformemente Acelerado'); 
grid on;
```


### Si deseas ampliar lo antes dicho puedes ver este video:

https://youtu.be/aDHVXyFXxCE?si=vv09IqGTXAV2S-eB

