using DifferentialEquations
using Plots 
pyplot()
using Parameters: @unpack

# Definición de parámetros del modelo
@kwdef struct HHParametros{FT}
    Cm::FT = 1.0       # Capacitancia de la membrana 
    gNa::FT = 120.0    # Conductancia máxima del sodio 
    gK::FT = 36.0      # Conductancia máxima del potasio 
    gL::FT = 0.3       # Conductancia de fuga 
    ENa::FT = 50.0     # Potencial de equilibrio del sodio 
    EK::FT = -77.0     # Potencial de equilibrio del potasio 
    EL::FT = -54.387   # Potencial de equilibrio de fuga 
    tr::FT = 0.5       # ms, tiempo de subida
    td::FT = 8         # ms, tiempo de bajada
    invtr::FT = 1 / tr # Es la inversa del tiempo de subida
    invtd::FT = 1 / td # Es la inversa del tiempo de bajada
    v0::FT = -20       # mV, que sera el potencial inicial
end

@kwdef mutable struct HH{FT}
    p::HHParametros{FT} = HHParametros{FT}() # Parámetros del modelo
    v::FT = -65.0       # Potencial de membrana inicial
    m::FT = 0.05        # Probabilidad de apertura de Na+
    h::FT = 0.6         # Probabilidad de cierre de Na+
    n::FT = 1.32        # Probabilidad de apertura de K+
end

# Funciones alfa y beta para las compuertas del sodio y potasio
function alfa_beta(V)
    am(V) = 0.1 * (V + 40) / (1 - exp(-0.1 * (V + 40)))
    bm(V) = 4.0 * exp(-0.0556 * (V + 65))
    ah(V) = 0.07 * exp(-0.05 * (V + 65))
    bh(V) = 1.0 / (1 + exp(-0.1 * (V + 35)))
    an(V) = 0.01 * (V + 55) / (1 - exp(-0.1 * (V + 55)))
    bn(V) = 0.125 * exp(-0.0125 * (V + 65))
    
    return (am, bm, ah, bh, an, bn)
end

# Uso de las funciones devueltas por alfa_beta
V = -65 # Ejemplo de potencial de membrana
(am, bm, ah, bh, an, bn) = alfa_beta(V)

# Mostrar los valores de las funciones
println("am = $(am(V))")
println("bm = $(bm(V))")
println("ah = $(ah(V))")
println("bh = $(bh(V))")
println("an = $(an(V))")
println("bn = $(bn(V))")

# Definición del sistema de ecuaciones diferenciales
function hodgkin_huxley!(du, u, p, t)
    V, m, h, n = u
    @unpack Cm, gNa, gK, gL, ENa, EK, EL, tr, td, invtr, invtd, v0 = p
    I = 10 * ((t > 50) && (t < 200)) + 35 * ((t > 250) && (t < 400)) # Corriente inyectada

    # Corrientes iónicas
    INa = gNa * m^3 * h * (V - ENa)
    IK = gK * n^4 * (V - EK)
    IL = gL * (V - EL)

    # Ecuaciones de Hodgkin-Huxley
    du[1] = (I - INa - IK - IL) / Cm
    du[2] = am(V) * (1 - m) - bm(V) * m
    du[3] = ah(V) * (1 - h) - bh(V) * h
    du[4] = an(V) * (1 - n) - bn(V) * n
end

# Inicializar modelo HH
hh = HH{Float64}()

# Condiciones iniciales
u0 = [hh.v, hh.m, hh.h, hh.n]

# Intervalo de tiempo
tspan = (0.0, 450.0)

# Definir y resolver el problema
prob = ODEProblem(hodgkin_huxley!, u0, tspan, hh.p)
sol = solve(prob, Tsit5())

# Graficar los resultados
p = plot(sol, idxs=(0, 1), xlabel="Tiempo (ms)", ylabel="Potencial de Membrana (mV)", label="V (mV)", title="Modelo de Hodgkin-Huxley")

# Cambiar la etiqueta del eje y para la primera serie de datos
ylabel!("V(mv)")

# Mostrar la gráfica
display(p)

nt = length(sol.t)
varr = sol[1,:]

spike = (varr[1:nt-1, :].<0) .& (varr[2:nt, :].>0)
num_spikes = sum(spike, dims=1)
println("Número de potenciales: ", num_spikes[1] )

