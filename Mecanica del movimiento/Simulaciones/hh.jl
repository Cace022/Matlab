#Este solo muestra los picos de acción, no la dinámica de las compuertas
using DifferentialEquations
using Plots
using Parameters: @unpack

# Definición de parámetros del modelo
@kwdef struct HHParameter{FT}
    Cm::FT = 1.0      # Capacitancia de la membrana (uF/cm^2)
    gNa::FT = 120.0   # Conductancia máxima del sodio (mS/cm^2)
    gK::FT = 36.0     # Conductancia máxima del potasio (mS/cm^2)
    gL::FT = 0.3      # Conductancia de fuga (mS/cm^2)
    ENa::FT = 50.0    # Potencial de equilibrio del sodio (mV)
    EK::FT = -77.0    # Potencial de equilibrio del potasio (mV)
    EL::FT = -54.387  # Potencial de equilibrio de fuga (mV)
end

# Funciones alfa y beta para las compuertas del sodio y potasio
function αm(V)
    return 0.1 * (V + 40) / (1 - exp(-0.1 * (V + 40)))
end

function βm(V)
    return 4.0 * exp(-0.0556 * (V + 65))
end

function αh(V)
    return 0.07 * exp(-0.05 * (V + 65))
end

function βh(V)
    return 1.0 / (1 + exp(-0.1 * (V + 35)))
end

function αn(V)
    return 0.01 * (V + 55) / (1 - exp(-0.1 * (V + 55)))
end

function βn(V)
    return 0.125 * exp(-0.0125 * (V + 65))
end

# Definición del sistema de ecuaciones diferenciales
function hodgkin_huxley!(du, u, p, t)
    V, m, h, n = u
    @unpack Cm, gNa, gK, gL, ENa, EK, EL = p
    I = 10 * ((t > 50) && (t < 200)) + 35 * ((t > 250) && (t < 400)) # Corriente inyectada

    # Corrientes iónicas
    INa = gNa * m^3 * h * (V - ENa)
    IK = gK * n^4 * (V - EK)
    IL = gL * (V - EL)

    # Ecuaciones de Hodgkin-Huxley
    du[1] = (I - INa - IK - IL) / Cm
    du[2] = αm(V) * (1 - m) - βm(V) * m
    du[3] = αh(V) * (1 - h) - βh(V) * h
    du[4] = αn(V) * (1 - n) - βn(V) * n
end

# Condiciones iniciales
u0 = [-65.0, 0.05, 0.6, 0.32]

# Parámetros del modelo
p = HHParameter{Float64}()

# Intervalo de tiempo
tspan = (0.0, 450.0)

# Definir y resolver el problema
prob = ODEProblem(hodgkin_huxley!, u0, tspan, p)
sol = solve(prob, Tsit5())

# Graficar los resultados
plot(sol, vars=(0, 1), xlabel="Tiempo (ms)", ylabel="Potencial de Membrana (mV)", label="V (mV)", title="Modelo de Hodgkin-Huxley")
plot!(sol, vars=(0, 2), label="m")
plot!(sol, vars=(0, 3), label="h")
plot!(sol, vars=(0, 4), label="n")
