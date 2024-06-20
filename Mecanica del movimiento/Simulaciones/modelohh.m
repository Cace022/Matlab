% Definición de parámetros del modelo
Cm = 1.0;      % Capacitancia de la membrana (uF/cm^2)
gNa = 120.0;   % Conductancia máxima del sodio (mS/cm^2)
gK = 36.0;     % Conductancia máxima del potasio (mS/cm^2)
gL = 0.3;      % Conductancia de fuga (mS/cm^2)
ENa = 50.0;    % Potencial de equilibrio del sodio (mV)
EK = -77.0;    % Potencial de equilibrio del potasio (mV)
EL = -54.387;  % Potencial de equilibrio de fuga (mV)

% Condiciones iniciales
u0 = [-65.0, 0.05, 0.6, 0.32];

% Intervalo de tiempo
tspan = [0.0, 450.0];

% Definir y resolver el problema
[t, u] = ode45(@(t, u) hodgkin_huxley(t, u, Cm, gNa, gK, gL, ENa, EK, EL), tspan, u0);

% Graficar los resultados
figure
plot(t, u(:,1))
xlabel('Tiempo (ms)')
ylabel('Potencial de Membrana (mV)')
title('Modelo de Hodgkin-Huxley')

function du = hodgkin_huxley(t, u, Cm, gNa, gK, gL, ENa, EK, EL)
    V = u(1); m = u(2); h = u(3); n = u(4);
    I = 10 * ((t > 50) && (t < 200)) + 35 * ((t > 250) && (t < 400)); % Corriente inyectada

    % Corrientes iónicas
    INa = gNa * m^3 * h * (V - ENa);
    IK = gK * n^4 * (V - EK);
    IL = gL * (V - EL);

    % Ecuaciones de Hodgkin-Huxley
    du = zeros(4,1);
    du(1) = (I - INa - IK - IL) / Cm;
    du(2) = am(V) * (1 - m) - bm(V) * m;
    du(3) = ah(V) * (1 - h) - bh(V) * h;
    du(4) = an(V) * (1 - n) - bn(V) * n;
end

function val = am(V)
    val = 0.1 * (V + 40) / (1 - exp(-0.1 * (V + 40)));
end

function val = bm(V)
    val = 4.0 * exp(-0.0556 * (V + 65));
end

function val = ah(V)
    val = 0.07 * exp(-0.05 * (V + 65));
end

function val = bh(V)
    val = 1.0 / (1 + exp(-0.1 * (V + 35)));
end

function val = an(V)
    val = 0.01 * (V + 55) / (1 - exp(-0.1 * (V + 55)));
end

function val = bn(V)
    val = 0.125 * exp(-0.0125 * (V + 65));
end