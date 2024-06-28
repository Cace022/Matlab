# Modelo de Hodgkin y Huxley

## Preliminares del concepto 
---

La membrana celular es una bicapa lipídica que aísla el interior de la célula del exterior y regula el paso de sustancias. Está compuesta principalmente de fosfolípidos, que son hidrofóbicos, junto con otros lípidos, proteínas y carbohidratos. Esta estructura actúa como una barrera selectiva, permitiendo el paso de algunas moléculas pequeñas y no cargadas, pero bloqueando las moléculas hidrosolubles y cargadas (iones). 

Para permitir el paso de iones, existen canales y bombas iónicas. Los canales iónicos son proteínas transmembrana que pueden ser de dos tipos: canales de fuga, que están siempre abiertos, y canales regulados, que se abren o cierran en respuesta a estímulos químicos (ligandos) o cambios de voltaje. Los canales son selectivos para ciertos iones, como sodio o potasio, y permiten el flujo de iones según su gradiente electroquímico, sin usar energía. El gradiente electroquímico de un ion es la combinación de dos fuerzas: la fuerza electrostática, que atrae iones hacia áreas con carga opuesta, y el gradiente de concentración, que mueve iones de áreas de alta a baja concentración.

Las bombas iónicas, como la bomba de sodio-potasio, usan energía (ATP) para transportar iones en contra de su gradiente electroquímico. Esta bomba mueve tres iones de sodio al exterior por cada dos iones de potasio al interior, manteniendo el exterior más rico en sodio y contribuyendo a la diferencia de potencial entre el interior y el exterior de la célula, conocido como potencial de membrana.

El potencial de membrana en reposo es el estado estable de una célula cuando no está transmitiendo señales. En este estado, el interior de la célula es más negativo en comparación con el exterior, debido principalmente a la acción de las bombas de sodio-potasio y la permeabilidad diferencial de la membrana a distintos iones. Este potencial de reposo es esencial para la excitabilidad de la célula y su capacidad de generar potenciales de acción.

Un potencial de acción es una señal eléctrica rápida y transitoria que se genera cuando el potencial de membrana supera un umbral específico, provocando la apertura de los canales de sodio regulados por voltaje y permitiendo una rápida entrada de sodio. Esto causa una rápida despolarización de la membrana. Posteriormente, los canales de sodio se inactivan y se abren los canales de potasio, permitiendo la salida de potasio y repolarizando la membrana. Finalmente, un breve periodo de hiperpolarización puede ocurrir antes de que la célula vuelva a su potencial de reposo.


El estado de equilibrio para un ion se alcanza cuando el flujo neto de ese ion a través de la membrana es cero, debido a que las fuerzas del gradiente de difusión y la fuerza electrostática se equilibran. Este valor se conoce como potencial de equilibrio del ion y puede calcularse mediante la ecuación de Nernst:

$$
E_{\text{ion}} = \frac{RT}{zF} \ln \left( \frac{[\text{ion}]_{\text{extracelular}}}{[\text{ion}]_{\text{intracelular}}} \right)
$$

Donde:

- `E_ion` es el potencial de equilibrio del ion
- `R` es la constante de los gases
- `T` es la temperatura en Kelvin
- `z` es la valencia del ion
- `F` es la constante de Faraday
- `[ion]_extracelular` es la concentración del ion fuera de la célula
- `[ion]_intracelular` es la concentración del ion dentro de la célula

---
## Explicacion del modelo 

En el modelo HH, existe un modelo de circuito eléctrico equivalente, que se denomina modelo de conductancia paralela de la membrana, en el que la membrana celular de la neurona se reemplaza como un condensador y el canal iónico enterrado en la membrana celular se reemplaza como una resistencia variable (una resistencia que cambia dinámicamente).

Un canal iónico es un tipo de transportador de membrana que pasa selectivamente a través de iones específicos (por ejemplo, iones de sodio, iones de potasio, etc.), y hay diferentes canales iónicos para cada tipo de ion (hay varios tipos de canales iónicos para el mismo ion)Hay conductancia (el recíproco de la resistencia, que significa la "facilidad de flujo" de la corriente) y potencial de equilibrio. En el modelo HH, se asumen los canales de sodio (Na), los canales de potasio (K) y los canales iónicos de corriente de fuga. En el proximo apartado realizaremos el modelo usando MATLAB y un apartado de como se haria en julia.

---


**Ahora con esto claro empezemos a trabajar:**

1. definimos los parametros de las propiedades eléctricas de la neurona modelada. Se agrupan en una estructura para facilitar su manejo y paso entre funciones.

```matlab
params = struct(...
    'Cm', 1.0, ...    % Capacitancia de la membrana
    'gNa', 120.0, ... % Conductancia máxima del sodio
    'gK', 36.0, ...   % Conductancia máxima del potasio
    'gL', 0.3, ...    % Conductancia de fuga
    'ENa', 50.0, ...  % Potencial de equilibrio del sodio
    'EK', -77.0, ...  % Potencial de equilibrio del potasio
    'EL', -54.387, ... % Potencial de equilibrio de fuga
    'v0', -65.0 ...   % mV, potencial inicial
);
```

Ahora biologicamente podriamos explicarlo de la siguiente manera: 

- Cm: Representa la capacidad de la membrana celular para almacenar carga eléctrica.
- gNa, gK, gL: Conductancias máximas de los canales iónicos de sodio, potasio y fuga.
- ENa, EK, EL: Potenciales de equilibrio para cada ion, determinados por sus concentraciones intra y extracelulares.
- v0: Potencial de reposo de la membrana.

2. Esta es una función que define un sistema de ecuaciones diferenciales ordinarias (ODEs). Cada ecuación describe la tasa de cambio de una variable de estado (V, m, h, n) en función de su estado actual y las entradas.

```matlab
% Función Hodgkin-Huxley
function dudt = hodgkin_huxley(t, u, params, I)
    V = u(1); m = u(2); h = u(3); n = u(4);
    
    am = @(V) 0.1 * (V + 40) / (1 - exp(-0.1 * (V + 40)));
    bm = @(V) 4.0 * exp(-0.0556 * (V + 65));
    ah = @(V) 0.07 * exp(-0.05 * (V + 65));
    bh = @(V) 1.0 / (1 + exp(-0.1 * (V + 35)));
    an = @(V) 0.01 * (V + 55) / (1 - exp(-0.1 * (V + 55)));
    bn = @(V) 0.125 * exp(-0.0125 * (V + 65));
    
    INa = params.gNa * m^3 * h * (V - params.ENa);
    IK = params.gK * n^4 * (V - params.EK);
    IL = params.gL * (V - params.EL);
    
    dVdt = (I - INa - IK - IL) / params.Cm;
    dmdt = am(V) * (1 - m) - bm(V) * m;
    dhdt = ah(V) * (1 - h) - bh(V) * h;
    dndt = an(V) * (1 - n) - bn(V) * n;
    
    dudt = [dVdt; dmdt; dhdt; dndt];
end
```

3. Usa un solver de ODEs (ode45) para calcular la evolución temporal del sistema. I_func simula pulsos de corriente en tiempos específicos.

```matlab
tspan = [0 450];
u0 = [params.v0, 0.05, 0.6, 0.32];
I_func = @(t) 10 * ((t > 50) & (t < 200)) + 35 * ((t > 250) & (t < 400));
[t, sol] = ode45(@(t, u) hodgkin_huxley(t, u, params, I_func(t)), tspan, u0);

figure(1);
plot(t, sol(:, 1));
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Simulación Estándar del Modelo Hodgkin-Huxley');
```

4. Similar a la simulación estándar, pero con un pulso de corriente negativa.

```matlab
tspan_rebound = [0 300];
I_rebound = @(t) -1 * ((t > 50) & (t < 150));
[t_rebound, sol_rebound] = ode45(@(t, u) hodgkin_huxley(t, u, params, I_rebound(t)), tspan_rebound, u0);

figure(2);
plot(t_rebound, sol_rebound(:, 1));
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Rebote Postinhibitorio');
```

5. Realiza múltiples simulaciones con diferentes corrientes constantes, contando los picos para calcular la frecuencia de disparo.

```matlab
tspan_FI = [0 1000];
I_range = 0:0.5:20;
frequencies = zeros(size(I_range));

for i = 1:length(I_range)
    I = I_range(i);
    [~, sol_FI] = ode45(@(t, u) hodgkin_huxley(t, u, params, I), tspan_FI, u0);
    num_spikes = count_spikes(sol_FI(:, 1));
    frequencies(i) = num_spikes / (tspan_FI(2) - tspan_FI(1)) * 1000;
end
```




























