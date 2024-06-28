% Definición de parámetros
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

% Función para contar picos
function num_spikes = count_spikes(V)
    threshold = 0;
    spike_times = find(V(1:end-1) < threshold & V(2:end) >= threshold);
    num_spikes = length(spike_times);
end

% 1. Simulación estándar
tspan = [0 450];
u0 = [params.v0, 0.05, 0.6, 0.32];
I_func = @(t) 10 * ((t > 50) & (t < 200)) + 35 * ((t > 250) & (t < 400));
[t, sol] = ode45(@(t, u) hodgkin_huxley(t, u, params, I_func(t)), tspan, u0);

figure(1);
plot(t, sol(:, 1));
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Simulación Estándar del Modelo Hodgkin-Huxley');

% 2. Rebote postinhibitorio
tspan_rebound = [0 300];
I_rebound = @(t) -1 * ((t > 50) & (t < 150));
[t_rebound, sol_rebound] = ode45(@(t, u) hodgkin_huxley(t, u, params, I_rebound(t)), tspan_rebound, u0);

figure(2);
plot(t_rebound, sol_rebound(:, 1));
xlabel('Tiempo (ms)');
ylabel('Potencial de Membrana (mV)');
title('Rebote Postinhibitorio');

% 3. Curva de frecuencia-corriente (F-I)
tspan_FI = [0 1000];
I_range = 0:0.5:20;
frequencies = zeros(size(I_range));

for i = 1:length(I_range)
    I = I_range(i);
    [~, sol_FI] = ode45(@(t, u) hodgkin_huxley(t, u, params, I), tspan_FI, u0);
    num_spikes = count_spikes(sol_FI(:, 1));
    frequencies(i) = num_spikes / (tspan_FI(2) - tspan_FI(1)) * 1000;
end

figure(3);
plot(I_range, frequencies, '-o');
xlabel('Corriente Inyectada (µA/cm^2)');
ylabel('Frecuencia de Disparo (Hz)');
title('Curva de Frecuencia-Corriente (F-I)');
grid on;