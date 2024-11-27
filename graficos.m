% Problelma 1

clear

% parametros

T_con = 4e-6; % [seg] periodo de conmutacion. no cambia la simulacion

% plots
t0 = T_con * 100; % multiplo del periodo de conmutacion para que se vea mas lindo
duracion_ciclos = 3;

tf = t0 + T_con * duracion_ciclos;


% simulacion
duty_cycle =  26.6667; % [%] original


%%
sim("Problema4_1")

% calculo valor medio
pasos_estacionario = 3000; % pasos para entrar en estado estacionario
% plot(I_inductor1.Data(3000:end))
tiempo_integracion = I_inductor1.Time(end) - I_inductor1.Time(pasos_estacionario);


%% corriente de entrada

% transitorio
figure
plot(I_input)
grid on;
xlim([0,  1e-4])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title("Transitorio de corriente de entrada")


% regimen estacionario
figure
plot(I_input)
grid on;
xlim([t0,  tf]) % Musestreo 6 ciclos de conmutacion
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title("Corriente de entrada")


%% pulsos

figure
hold on
grid on;
plot(pulses, 'LineWidth', 1.5)
xlim([t0, tf])
ylim([0, 1.05])
title("Señal PWM");
ylabel("Señal de control")
xlabel("Tiempo [seg]")

%% inductores

figure

% tension
title('Tensión y corriente por los inductores')
subplot(211)
hold on
grid on;
plot(V_inductor1)
plot(V_inductor2)
xlim([t0, tf])
ylim([-10, 23])
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")
legend('L1', 'L2')
title('Tensión en el inductor')

% corriente
% valor medio

% calculo valor medio como 1/T * int_0^T I(t) dt
% si uso la funcion mean da un resultado incorrecto porque el paso de
% integracion es variable
media_I_inductor_1 = (1/tiempo_integracion)*trapz(I_inductor1.Time(pasos_estacionario:end), I_inductor1.Data(pasos_estacionario:end));
media_I_inductor_2 = (1/tiempo_integracion)*trapz(I_inductor2.Time(pasos_estacionario:end), I_inductor2.Data(pasos_estacionario:end));

media_I_inductor_1_plot = media_I_inductor_1 * ones(1,2);
media_I_inductor_2_plot = media_I_inductor_2 * ones(1,2);

subplot(212)
hold on
grid on;
plot(I_inductor1)
plot(I_inductor2)
plot([t0, tf], media_I_inductor_1_plot, 'b--')
plot([t0, tf], media_I_inductor_2_plot, 'r--')
legend('L1', 'L2', 'Val Medio L1', 'Val Medio L2')
xlim([t0, tf])
ylim([7.4, 8.6])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title('Corriente en el inductor')

% corriente en el transitorio
% explica porque una es mas grande que la otra

figure
subplot(211)
hold on
grid on;
plot(V_inductor1,'LineWidth', 1.5)
plot(V_inductor2,'LineWidth', 1.5)
xlim([0, 10*T_con])
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")
legend('L1', 'L2')
title('Tensión en el transitorio')

subplot(212)
hold on
grid on;
plot(I_inductor1,'LineWidth', 1.5)
plot(I_inductor2,'LineWidth', 1.5)
xlim([0, 10*T_con])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
legend('L1', 'L2')
title('Corriente en el transitorio')

%% salida


figure


% calculo valor medio como 1/T * int_0^T I(t) dt
% si uso la funcion mean da un resultado incorrecto porque el paso de
% integracion es variable
media_V_salida = (1/tiempo_integracion)*trapz(V_output.Time(pasos_estacionario:end), V_output.Data(pasos_estacionario:end));
media_I_salida = (1/tiempo_integracion)*trapz(I_output.Time(pasos_estacionario:end), I_output.Data(pasos_estacionario:end));

media_V_salida_plot = media_V_salida * ones(1,2);
media_I_salida_plot = media_I_salida * ones(1,2);


subplot(211)
hold on
grid on
plot(V_output)
plot([t0, tf], media_V_salida_plot, 'b--')
xlim([t0, tf])
ylim([7.991, 7.997])
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")
title("Tensión de salida")

subplot(212)
hold on
grid on
plot(I_output)
plot([t0, tf], media_I_salida_plot, 'b--')
xlim([t0, tf])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title("Corriente de salida")

%% calculo ripple

ripple_V = max(V_output.Data(3000:end)) - min(V_output.Data(3000:end)) % pico a pico [V]
ripple_I = max(I_output.Data(3000:end)) - min(I_output.Data(3000:end)) % pico a pico [A]

%% Corriente en los capacitores

figure
hold on
grid on
plot(I_cap1)
plot(I_cap2, '.')
xlim([t0, tf])
legend('C1', 'C2')
