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
pasos_estacionario = 6000; % pasos para entrar en estado estacionario
% plot(I_inductor1.Data(5000:end))
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

ripple_V = max(V_output.Data(pasos_estacionario:end)) - min(V_output.Data(pasos_estacionario:end)) % pico a pico [V]
ripple_I = max(I_output.Data(pasos_estacionario:end)) - min(I_output.Data(pasos_estacionario:end)) % pico a pico [A]

%% senoide para comparar con salida
t = t0:(tf-t0)/1000:tf;
amplitud = max(V_output.Data(pasos_estacionario:end)) - media_V_salida;
senoide = media_V_salida + (amplitud) * sin(t*(1/2e-6)*2*pi+pi);

figure
hold on
grid on
plot(V_output)
plot(t, senoide)
xlim([t0, tf])
ylim([7.991, 7.997])
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")

%% Corriente en los capacitores

figure
hold on
grid on
plot(I_cap1, 'LineWidth', 1)
plot(I_cap2, '.', 'LineWidth', 2)
xlim([t0, tf])
legend('C1', 'C2')
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title("Corriente por los capacitores")


%% Graficos del ripple para distintos duty cycle
% tarda mucho en ejecutarse

iteraciones = 21;

datos_V_ripple = zeros(1, iteraciones);
datos_I_ripple = zeros(1, iteraciones);

datos_V_normalizado = zeros(1, iteraciones);

datos_ripple_buck_normal = zeros(1, iteraciones);

datos_V_output = zeros(1, iteraciones);


iter_duty_cycle = 1:(98)/(iteraciones-1):99;

indx=1:iteraciones;

for i = indx
    duty_cycle = iter_duty_cycle(i);
    sim('Problema4_1')
    
    % como tengo paso de simulacion variable, el paso de simulucion cambia
    % con el duty cycle
    pasos_estacionario = find(V_output.Time > 5e-4, 1); 
    tiempo_integracion = I_inductor1.Time(end) - I_inductor1.Time(pasos_estacionario);
                                       
    datos_V_ripple(i) = max(V_output.Data(pasos_estacionario:end)) - min(V_output.Data(pasos_estacionario:end));
    datos_I_ripple(i) = max(I_output.Data(pasos_estacionario:end)) - min(I_output.Data(pasos_estacionario:end)); % pico a pico [A]

    V_salida = (1/tiempo_integracion)*trapz(V_output.Time(pasos_estacionario:end), V_output.Data(pasos_estacionario:end));
    
    datos_V_output(i) = V_salida;
    
    datos_V_normalizado(i) = datos_V_ripple(i)/V_salida;
   
    
end


figure
hold on
grid on
yyaxis left
plot(iter_duty_cycle, datos_V_ripple)
yyaxis right
plot(iter_duty_cycle, datos_I_ripple)


% ripple normalizado calculado para un buck normal
% para que sea igual hay q dividir por 4, que es lo mimso que decir que
% tengo el doble de capacidad y el doble de inductacia. Tiene sentido tener
% el doble de capaciada pq tengo C1 y C2 en paralelo, pero el doble de
% inductancia?
datos_ripple_buck_normal = (1/8) * ((T_con/2)^2*(1-(iter_duty_cycle/100)))/(11e-6*38.8e-6);

figure
grid on
hold on
plot(iter_duty_cycle, datos_V_normalizado, 'LineWidth', 1.5)
plot(iter_duty_cycle, datos_ripple_buck_normal, 'LineWidth', 1.5)
xlabel('Duty Cycle \delta [%]')
ylabel('Ripple normalizado {\Delta}V_{out}/V_{out} [V]')
legend('Multifase', 'Convencional')


%% que pasa con el ripple cuando el I es 50%?
duty_cycle = 50;

sim('Problema4_1')


pasos_estacionario = find(V_output.Time > 8e-4, 1); 
tiempo_integracion = I_inductor1.Time(end) - I_inductor1.Time(pasos_estacionario);

media_V_salida = (1/tiempo_integracion)*trapz(V_output.Time(pasos_estacionario:end), V_output.Data(pasos_estacionario:end));
media_I_salida = (1/tiempo_integracion)*trapz(I_output.Time(pasos_estacionario:end), I_output.Data(pasos_estacionario:end));

media_V_salida_plot = media_V_salida * ones(1,2);
media_I_salida_plot = media_I_salida * ones(1,2);

% Nuevos parametros
t0 = tout(pasos_estacionario);
tf = t0 + T_con*duracion_ciclos;


figure

subplot(211)
hold on
grid on
plot(V_output)
plot([t0, tf], media_V_salida_plot, 'b--')
xlim([t0, tf])
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


% corrientes
figure

hold on
grid on

plot(I_inductor1, 'LineWidth', 1.5)
plot(I_inductor2, 'LineWidth', 1.5)
xlim([t0, tf])
legend('L1', 'L2')
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Problema 2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


clear

% parametros
T_con = 4e-6; % [seg] periodo de conmutacion. no cambia la simulacion

% plots
t0 = T_con * 1000; % multiplo del periodo de conmutacion para que se vea mas lindo
duracion_ciclos = 3;

tf = t0 + T_con * duracion_ciclos;


% simulacion
duty_cycle =  53.84; % [%] original


%% 

sim('Problema4_2')

pasos_estacionario = find(tout > 0.002, 1);

tiempo_integracion = tout(end) - tout(pasos_estacionario);


%% corriente de entrada

% transitorio
figure
plot(I_input,'LineWidth',1.5)
grid on;
xlim([0,  2e-4])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title("Transitorio de corriente de entrada")


% regimen estacionario
figure
plot(I_input, 'LineWidth',1.5)
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
ylabel("Señal de control")
xlabel("Tiempo [seg]")

%% Inductores

figure

% tension
subplot(211)
hold on
grid on;
plot(V_L1 , 'LineWidth',1.2)
plot(V_L2 , 'LineWidth',1.2)
xlim([t0, tf])
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")
legend('L1', 'L2')
title('Tensión en el inductor')

% corriente
% valor medio

% calculo valor medio como 1/T * int_0^T I(t) dt
% si uso la funcion mean da un resultado incorrecto porque el paso de
% integracion es variable
media_I_inductor_1 = (1/tiempo_integracion)*trapz(I_L1.Time(pasos_estacionario:end), I_L1.Data(pasos_estacionario:end));
media_I_inductor_2 = (1/tiempo_integracion)*trapz(I_L2.Time(pasos_estacionario:end), I_L2.Data(pasos_estacionario:end));

media_I_inductor_1_plot = media_I_inductor_1 * ones(1,2);
media_I_inductor_2_plot = media_I_inductor_2 * ones(1,2);

subplot(212)
hold on
grid on;
plot(I_L1 , 'LineWidth',1.2)
plot(I_L2 , 'LineWidth',1.2)
% plot([t0, tf], media_I_inductor_1_plot, 'b--', 'LineWidth', 1.5)
% plot([t0, tf], media_I_inductor_2_plot, 'r--', 'LineWidth', 0.8)
% legend('L1', 'L2', 'Val Medio L1', 'Val Medio L2')
legend('L1','L2')
xlim([t0, tf])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")
title('Corriente en el inductor')

%% Capacitores


% corriente por los capacitores

figure
subplot(211)

plot(pulses, 'LineWidth', 1.5)
hold on
grid on;
xlim([t0, tf])
ylim([0, 1.05])
ylabel("Señal de control")
xlabel("Tiempo [seg]")

subplot(212)
plot(I_C1, 'LineWidth', 1)
hold on
grid on
plot(I_C2, '.', 'LineWidth', 3)
xlim([t0, tf])
legend('C1', 'C2')
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")

%% Salida

figure

media_V_salida = (1/tiempo_integracion)*trapz(tout(pasos_estacionario:end), V_output.Data(pasos_estacionario:end));
media_I_salida = (1/tiempo_integracion)*trapz(tout(pasos_estacionario:end), I_output.Data(pasos_estacionario:end));

media_V_salida_plot = media_V_salida * ones(1,2);
media_I_salida_plot = media_I_salida * ones(1,2);


subplot(211)
hold on
grid on
plot(pulses, 'LineWidth', 1)
plot(V_output, 'LineWidth', 1.5, 'Color', 'r')
plot([t0, tf], media_V_salida_plot, 'b--', 'Color', 'r')
ylabel("Tensión [V]")
xlabel("Tiempo [seg]")
xlim([t0, tf])

ylim([25.94, 26.01])


subplot(212)
hold on
grid on
plot(I_output, 'LineWidth', 1.5)
plot([t0, tf], media_I_salida_plot, 'b--')
xlim([t0, tf])
ylabel("Corriente [A]")
xlabel("Tiempo [seg]")

%% Calculo de ripple

ripple_V = max(V_output.Data(pasos_estacionario:end)) - min(V_output.Data(pasos_estacionario:end))
ripple_I = max(I_output.Data(pasos_estacionario:end)) - min(I_output.Data(pasos_estacionario:end))


%% Duty cicle


iteraciones = 21;

datos_V_ripple = zeros(1, iteraciones);
datos_I_ripple = zeros(1, iteraciones);

datos_V_normalizado = zeros(1, iteraciones);

datos_V_output = zeros(1, iteraciones);


iter_duty_cycle = 1:(89)/(iteraciones-1):90;

indx=1:iteraciones;

figure
hold on

for i = indx
    duty_cycle = iter_duty_cycle(i);
    sim('Problema4_2')
    
    % como tengo paso de simulacion variable, el paso de simulucion cambia
    % con el duty cycle
    pasos_estacionario = find(tout > 2e-3, 1); 
    tiempo_integracion = tout(end) - tout(pasos_estacionario);
                                       
    datos_V_ripple(i) = max(V_output.Data(pasos_estacionario:end)) - min(V_output.Data(pasos_estacionario:end));
    datos_I_ripple(i) = max(I_output.Data(pasos_estacionario:end)) - min(I_output.Data(pasos_estacionario:end)); % pico a pico [A]

    V_salida = (1/tiempo_integracion)*trapz(V_output.Time(pasos_estacionario:end), V_output.Data(pasos_estacionario:end));
    
    datos_V_output(i) = V_salida;
    
    datos_V_normalizado(i) = datos_V_ripple(i)/V_salida;
    
    plot(V_output - V_salida)
   
    
end

xlim([t0, tf])
legend(string(iter_duty_cycle))

figure
hold on
grid on
yyaxis left
plot(iter_duty_cycle, datos_V_ripple)
yyaxis right
plot(iter_duty_cycle, datos_I_ripple)

% rippleV/Vout = ( duty_cycle*T_con )/( RL * C )
datos_ripple_boost_normal = ( (iter_duty_cycle/100) * T_con/2 )/( 1.625 * 22e-6 );

figure
grid on
hold on
plot(iter_duty_cycle, datos_V_normalizado, 'LineWidth', 1.5)
plot(iter_duty_cycle, datos_ripple_boost_normal, 'LineWidth', 1.5)
xlabel('Duty Cycle \delta [%]')
ylabel('Ripple normalizado {\Delta}V_{out}/V_{out} [V]')
legend('Multifase', 'Convencional')