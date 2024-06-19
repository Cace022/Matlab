# Modelo de Hodgkin y Huxley
---

Trataré de explicar de manera muy sencilla este modelo para aquellos que no tienen conocimientos en neurociencia, compactando la parte computacional, biológica y matemática. Para esto, exploraremos a su vez un lenguaje open source como **Julia**.

El modelo de Hodgkin y Huxley describe cómo las neuronas generan y propagan señales eléctricas conocidas como potenciales de acción. Estas señales son fundamentales para que el cerebro y el sistema nervioso funcionen correctamente, permitiendo la comunicación entre células nerviosas. El modelo utiliza ecuaciones matemáticas para explicar estos procesos y ha sido crucial para nuestro entendimiento de la neurociencia.

Para la realizacion de este modelo tomaremos en cuenta la revista **Spiking Neural Networks**: 
https://compneuro-julia.github.io/_static/pdf/SNN_from_scratch_with_python_ver2_1.pdf

---
## Preliminares del concepto 

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

## Modelo en Julia y MATLAB

Julia y MATLAB son ambos lenguajes de programación que se utilizan principalmente para cálculos numericos y cientificos:

- Julia es un lenguaje de programación de alto nivel y alto rendimiento para cálculos técnicos. Se diseñó para abordar la necesidad de velocidad, facilidad de uso y capacidad para programación de alto rendimiento. Julia es popular en la investigación científica y en áreas que requieren cálculos matemáticos, como la física, y la ingeniería.

- MATLAB es un lenguaje de programación y entorno de desarrollo que se utiliza ampliamente en ciencia e ingeniería. MATLAB es conocido por su capacidad para realizar cálculos con matrices y vectores, lo que lo hace muy útil para la matemática lineal, el análisis de datos y la modelización y simulación numérica.

**Ahora con esto claro empezemos a trabajar:**

En julia es necesario importar paquetes adicionales y configurar algunas propiedades de los gráficos. En MATLAB por lo general ya vienen funcionalidades equivalentes.

``` julia
using Base: @kwdef  # Importa @kwdef para definir estructuras con valores predeterminados
using Parameters: @unpack  # Importa @unpack para extraer valores de los campos de una estructura
using PyPlot  # Importa PyPlot para hacer gráficos
using DifferentialEquations #importa metodos de ecuaciones diferenciales 
rc("axes.spines", top=false, right=false)  # Elimina los bordes superior y derecho de los gráficos
```

Ahora definiremos las contantes que se utilizaran en el modelo. Algunas de estas pueden ser fijas mientras que otras púeden variar. En el modelo de HH son cuatro las que pueden variar (v, m, h, n) y estas representan la dinámica presináptica (es el conjunto de procesos que ocurren en la terminal nerviosa antes de liberar los neurotransmisores hacia la siguiente neurona).

#### Julia

``` julia
@kwdef struct parametrosHH{FT}
    Cm::FT = 1 # Capacitancia de la menbrana (uF/cm^2)
    gNa::FT = 120; gk::FT = 36; gL::FT = 0.3 # Conductancias máximas para Na+, K+, fuga (mS/cm^2)
    ENa::FT = 50; Ek::FT = -77; EL::FT = -54 # Potenciales de equilibrio para Na+, K+, fuga (mV)
    tr::FT = 0.5; td::FT = 8 # Constantes de tiempo para la recuperación y la decadencia (ms)
    invtr::FT = 1/tr; invtd::FT = 1/td # Estos valores inversos de se utilizan para ahorrar tiempo de computación
    v0::FT = -20 # Potencial de membrana inicial (mV)
end

@kwdef mutable struct HH{FT}
    param::parametrosHH = parametrosHH{FT}() #Parametros para el modelo HH
    N::UInt16 # Numero de neoruonas o puntos de tiempo (UInt16 representa un numero entero sin signos de 16 bits)
    v::Vector{FT} = fill(-0.65, N) # Potencial de membrana para cada neurona/punto de tiempo
    m::Vector{FT} = fill(0.05, N); h::Vector{FT} = fill(0.6, N) # Probabilidades para los canales de Na+
    n::Vector{FT} = fill(0.32, N); r::Vector{FT} = zeros(N) # Probabilidades para los canales de K+ y variable adicional
end
```

#### MATLAB

``` matlab
```

ahora utilizaremos el metodo de Runge-Kutta (ode45 y solve) 

#### Julia 

``` julia
function hodgkin_huxley!(du, u, p, t)
    @unpack N, v, m, h, n, r = u
    @unpack Cm, gNa, gK, gL, ENa, EK, EL, invtr, invtd, v0 = p

    @inbounds for i = 1:N
        du[i] = ((0.1(v[i]+40)/(1 - exp(-0.1(v[i]+40))))*(1 - m[i]) - 4exp(-(v[i]+65) / 18)*m[i])
        du[i] = ((0.07exp(-0.05(v[i]+65)))*(1.0 - h[i]) - 1/(1 + exp(-0.1(v[i]+35)))*h[i])
        du[i] = ((0.01(v[i]+55)/(1 - exp(-0.1(v[i]+55))))*(1 - n[i]) - (0.125exp(-0.0125(v[i]+65)))*n[i])
        du[i] = 1 / Cm * (Ie[i] - gNa * m[i]^3 * h[i] * (v[i] - ENa) - gK * n[i]^4 * (v[i] - EK) - gL * (v[i] - EL))
        du[i] = ((invtr - invtd) * (1 - r[i])/(1 + exp(-v[i] + v0)) - r[i] * invtd)
    end
end

# Define las condiciones iniciales y los parámetros
u0 = [v0, m0, h0, n0, r0]
p = [Cm, gNa, gK, gL, ENa, EK, EL, invtr, invtd, v0]

# Define el intervalo de tiempo para la solución
tspan = (0.0, 1.0)

# Resuelve la ecuación diferencial
prob = ODEProblem(hodgkin_huxley!, u0, tspan, p)
sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8) 
```

#### MATLAB 

``` matlab
``` 

Ahora ejecutaremos la simulación 

#### Julia 

``` julia
T = 450 
dt = 0.05 
nt = Int(T/dt) 
N = 1 

t = (1:nt)*dt
Ie = repeat(10f0 * ((t .> 50) - (t .> 200)) + 35f0 * ((t .> 250) - (t .> 400)), 1, N)  # injection current

varr, gatearr = zeros(nt, N), zeros(nt, 3, N)

neurons = HH{Float32}(N=N)

@time for i = 1:nt
    update!(neurons, neurons.param, Ie[i, :], dt)
    varr[i, :] = neurons.v
    gatearr[i, :, :] .= [neurons.m; neurons.h; neurons.n]
end
```

#### MATLAB

``` matlab
```

ahora dibujaremos el potencial de menbrana de la neurona 

#### Julia 

``` julia 
figure(figsize=(6, 5))
subplot(3,1,1); plot(t, varr[:, 1], color="black"); ylabel("V (mV)")
subplot(3,1,2); labellist=["m" "h" "n"] 
for i in 1:3
    plot(t, gatearr[:, i, 1], label=labellist[i])
end; 
ylabel("Gating Value"); legend()
subplot(3,1,3); plot(t, Ie[:, 1], color="black"); ylabel(L"Current($\mu$A/cm$^2$)"); xlabel("Times (ms)")
tight_layout()
```

#### MATLAB 

``` matlab
```

































