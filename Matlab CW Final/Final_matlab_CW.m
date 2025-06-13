clc; clear; close all;

format long
s = settings;
s.matlab.fonts.editor.code.Name.TemporaryValue = 'Calibri';
set(groot,'defaultLineLineWidth',2)  %sets graph line width as 2
set(groot,'defaultAxesFontSize',24)  %sets graph axes font size as 18
set(groot,'defaulttextfontsize',24)  %sets graph text font size as 18
set(groot,'defaultLineMarkerSize',8) %sets line marker size as 8
set(groot,'defaultAxesXGrid','on')   %sets X axis grid on 
set(groot,'defaultAxesYGrid','on')   %sets Y axis grid on
set(groot,'DefaultAxesBox', 'on')   %sets Axes boxes on

picturewidth = 20; % set this parameter and keep it forever
hw_ratio = 0.65; % feel free to play with this ratio
    
%% Q1

% Parameters
g = 9.81;
L = 1;
density = 500;
density_L = 1000;   % Density at 

% Initial Values
y_0 = 0.1;          % Displacement
dydt_0 = 0;         % Velocity

% Analytical Parameters
omega_2 = (density_L * g) / (density * L);
omega = sqrt(omega_2);

t = 0:0.01:10;
analytical_vals = y_0 * cos(omega .* t) + (dydt_0 / omega) * sin(omega .* t); 

% plot Analytical solution
Analytical_solution = figure;
%set(Analytical_solution,"WindowState","maximized");
set(findall(Analytical_solution,'-property','FontSize'),'FontSize',24);
set(findall(Analytical_solution,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Analytical_solution,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Analytical_solution,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Analytical_solution,'Position');
set(Analytical_solution,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(t,analytical_vals,'k')
hold on
%title("Figure 1: Analytical Solution")
xlabel("Time, (s)")
ylabel("Displacement, (m)")

grid minor
hold off

saveas(Analytical_solution,'Analytical_solution.svg');

%% Q1c:
% Convert 2nd order ode to 1st order system

A = [0, 1; -omega_2, 0];

% Eigen Values
eigen_vals = eig(A);
disp("Eigenvalues of matrix A:");
disp(eigen_vals);

%% Q2a:

% Initialising Values
Y_0 = [y_0; dydt_0];
delta_t = 0.1;
time_range = 0:delta_t:10;

% Calling Explicit Euler Function
Y_Euler = Explicit_Euler(Y_0, A, delta_t, time_range);

% plot Euler solution
Euler_Sol = figure;
%pbaspect([2 1 1])
%set(Euler_Sol,"WindowState","maximized");
set(findall(Euler_Sol,'-property','FontSize'),'FontSize',24);
set(findall(Euler_Sol,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Euler_Sol,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Euler_Sol,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Euler_Sol,'Position');
set(Euler_Sol,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
plot(time_range,Y_Euler(1, :),'r--');
hold on
%title("Figure 2: Explicit Euler: Displacement vs Time")
xlabel("Time, (s)")
ylabel("Displacement, (m)")

grid minor
hold off

saveas(Euler_Sol,'Euler_Sol.svg');

%% Q2b: Runge-Kutta 4th Order Method

% calling Runge-Kutta-4 Function
Y_RK4 = RK4(Y_0, A, delta_t, time_range);

% plot RK4 solution
Runge_Kutta_4 = figure;
%pbaspect([2 1 1])
%set(Runge_Kutta_4,"WindowState","maximized");
set(findall(Runge_Kutta_4,'-property','FontSize'),'FontSize',24);
set(findall(Runge_Kutta_4,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Runge_Kutta_4,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Runge_Kutta_4,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Runge_Kutta_4,'Position');
set(Runge_Kutta_4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
plot(time_range,Y_RK4(1, :),'b-.');
hold on
%title("Figure 2: Runge-Kutta 4th Order: Displacement vs Time")
xlabel("Time, (s)")
ylabel("Displacement, (m)")

grid minor
hold off

saveas(Runge_Kutta_4,'Runge_Kutta_4.svg');

%% Q2c: MATLAB ode45 function

ode_func = @(t, Y_0) [Y_0(2); -omega_2 * Y_0(1)];
[t_ode45, Y_ode45] = ode45(ode_func, [0, 10], Y_0);

ODE45 = figure;
%pbaspect([2 1 1])
%set(ODE45,"WindowState","maximized");
set(findall(ODE45,'-property','FontSize'),'FontSize',24);
set(findall(ODE45,'-property','Interpreter'),'Interpreter','latex') 
set(findall(ODE45,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(ODE45,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(ODE45,'Position');
set(ODE45,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
plot(t_ode45,Y_ode45(:, 1),'c:');
hold on
%title("Figure 2: MATLAB ode45 Solution: Displacement vs Time")
xlabel("Time, (s)")
ylabel("Displacement, (m)");
grid minor
hold off

saveas(ODE45,'ODE45.svg');

%% Q2d: Comparison of Methods

Comparison_of_methods = figure;
%pbaspect([2 1 1])
%set(Comparison_of_methods,"WindowState","maximized");
set(findall(Comparison_of_methods,'-property','FontSize'),'FontSize',24);
set(findall(Comparison_of_methods,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Comparison_of_methods,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Comparison_of_methods,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Comparison_of_methods,'Position');
set(Comparison_of_methods,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(t, analytical_vals, 'k');
hold on
plot(time_range, Y_Euler(1, :), 'r--');
plot(time_range, Y_RK4(1, :), 'b-.');
plot(t_ode45, Y_ode45(:,1), 'g:');
legend('Analytical', 'Eular', 'Runge-Kutta', 'ode45',location = 'best');
xlabel('Time (s)');
ylabel('Displacement y (m)');
%title('Comparison of Solutions');
grid on;

saveas(Comparison_of_methods,'Comparison_of_methods.svg');

%%

Comparison_of_methods_2 = figure;
%pbaspect([2 1 1])
%set(Comparison_of_methods,"WindowState","maximized");
set(findall(Comparison_of_methods_2,'-property','FontSize'),'FontSize',24);
set(findall(Comparison_of_methods_2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Comparison_of_methods_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Comparison_of_methods_2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Comparison_of_methods_2,'Position');
set(Comparison_of_methods_2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(t, analytical_vals, 'k');
hold on
plot(time_range, Y_Euler(1, :), 'r--');
plot(time_range, Y_RK4(1, :), 'b-.');
plot(t_ode45, Y_ode45(:,1), 'g:');
legend('Analytical', 'Eular', 'Runge-Kutta', 'ode45',location = 'north');
xlabel('Time (s)');
ylabel('Displacement y (m)');
axis([0,1.5,-0.2,0.4])
%title('Comparison of Solutions');
grid on;

saveas(Comparison_of_methods_2,'Comparison_of_methods_2.svg');

%% Q2e: Velocity vs Position
Velocity_vs_Position = figure;
%set(Velocity_vs_Position,"WindowState","maximized");
set(findall(Velocity_vs_Position,'-property','FontSize'),'FontSize',24);
set(findall(Velocity_vs_Position,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Velocity_vs_Position,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Velocity_vs_Position,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Velocity_vs_Position,'Position');
set(Velocity_vs_Position,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(Y_Euler(1, :), Y_Euler(2, :), 'r--');
hold on;
plot(Y_RK4(1, :), Y_RK4(2, :), 'b-.');
plot(Y_ode45(:, 1), Y_ode45(:, 2), 'g:');
legend('Euler','Runge-Kutta', 'ode45',location = "best");
xlabel('Position, y (m)');
ylabel('Velocity, v (m)');
%title('Velocity vs Position');
grid on;

saveas(Velocity_vs_Position, 'Velocity_vs_Position.svg');

%% without dy/dt vs y Euler

Velocity_vs_Position2 = figure;
%set(Velocity_vs_Position2,"WindowState","maximized");
set(findall(Velocity_vs_Position2,'-property','FontSize'),'FontSize',24);
set(findall(Velocity_vs_Position2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Velocity_vs_Position2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Velocity_vs_Position2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Velocity_vs_Position2,'Position');
set(Velocity_vs_Position2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

%plot(Y_Euler(1, :), Y_Euler(2, :), 'r--');
hold on;
plot(Y_RK4(1, :), Y_RK4(2, :), 'b-.');
plot(Y_ode45(:, 1), Y_ode45(:, 2), 'g:');
legend('Runge-Kutta', 'ode45',location="best");
xlabel('Position y (m)');
ylabel('Velocity v (m)');
%title('Velocity vs Position');
grid on;
hold off

saveas(Velocity_vs_Position2, 'Velocity_vs_Position2.svg');


%% Q3: Numerical Stability
%% Q3a: Error plots

t__ = time_range;
analytical_values_for_error = y_0 * cos(omega .* t__) + (dydt_0 / omega) * sin(omega .* t__); 
analytical_values_for_error_ode45 = y_0 * cos(omega .* t_ode45) + (dydt_0 / omega) * sin(omega .* t_ode45); 

eular_error = abs(analytical_values_for_error-Y_Euler(1,:));
rk4_error = abs(analytical_values_for_error-Y_RK4(1,:));
ode45_error = abs(analytical_values_for_error_ode45-Y_ode45(:,1));

error_plots = figure;
%set(error_plots, "windowstate", "maximized");
set(findall(error_plots,'-property','FontSize'),'FontSize',24);
set(findall(error_plots,'-property','Interpreter'),'Interpreter','latex') 
set(findall(error_plots,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(error_plots,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(error_plots,'Position');
set(error_plots,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

semilogy(time_range,eular_error,'r--');
hold on
semilogy(time_range,rk4_error,'b-.');
semilogy(t_ode45,ode45_error,'g:');
xlabel("Time (s)");
ylabel('Displacement error (m)');
legend('Explicit Euler','Runge-Kutta 4th Order','ode45',Location='best');
hold off

saveas(error_plots,'error_plots.svg')

%% 3b: step size study

stepsizes = [0.5, 0.05, 0.005];
range = 0:stepsizes(3):10;
RK4_step_study = zeros(3,length(range));
Euler_step_study = zeros(3,length(range));
time__step = zeros(3,length(range));


for i = [1,2,3]

    time_range_step_size = 0:stepsizes(i):10;
    time__step(i,1:length(time_range_step_size)) = time_range_step_size;

    rk_temp = RK4(Y_0, A, stepsizes(i), time_range_step_size);

    analytical_values_for_error_1 = y_0 * cos(omega .* time_range_step_size) + (dydt_0 / omega) * sin(omega .* time_range_step_size); 

    RK4_step_study(i,1:length(time_range_step_size)) = abs(rk_temp(1,:) - analytical_values_for_error_1);

    %RK4_error(i,:) = abs(rk_temp(1,:)-analytical_values_for_error);


    Euler_temp = Explicit_Euler(Y_0, A, stepsizes(i), time_range_step_size);

    %Euler_error(i,:) = abs(Euler_temp(1,:)-analytical_values_for_error);

    Euler_step_study(i,1:length(time_range_step_size)) = abs(Euler_temp(1,:) - analytical_values_for_error_1);
end

time__step(time__step(:,2:end) == 0) = "Nan";
RK4_step_study(RK4_step_study(:,2:end) == 0) = "Nan";

Euler_step_study(Euler_step_study(:,2:end) == 0) = "Nan";

%% 3b: step size study

% Euler_time_step_study = figure;
% %set(Euler_time_step_study,"WindowState","maximized");
% set(findall(Euler_time_step_study,'-property','FontSize'),'FontSize',24);
% set(findall(Euler_time_step_study,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(Euler_time_step_study,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(Euler_time_step_study,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(Euler_time_step_study,'Position');
% set(Euler_time_step_study,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% 
% %plot(time__step(1,:),Euler_step_study(1,:),"r--");
% 
% semilogy(time__step(2,:),Euler_step_study(2,:),"b-.");
% hold on;
% semilogy(time__step(3,:),Euler_step_study(3,:),"#D95319:");
% legend('\Delta t = 0.5','\Delta t = 0.05','\Delta t = 0.005',location = 'best');
% 
% xlabel('Time (s)');
% ylabel('Displacement error (m)');
% grid on;
% hold off
% 
% saveas(Euler_time_step_study, 'Euler_time_step_study.svg');
% 

%% without large growth
Euler_time_step_study_2 = figure;
%set(Euler_time_step_study,"WindowState","maximized");
set(findall(Euler_time_step_study_2,'-property','FontSize'),'FontSize',24);
set(findall(Euler_time_step_study_2,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Euler_time_step_study_2,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Euler_time_step_study_2,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Euler_time_step_study_2,'Position');
set(Euler_time_step_study_2,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

semilogy(time__step(1,:),Euler_step_study(1,:),"r--");
hold on;
semilogy(time__step(2,:),Euler_step_study(2,:),"b-.");
semilogy(time__step(3,:),Euler_step_study(3,:),"g:");
legend('\Deltat =0.5','\Deltat =0.05','\Deltat =0.005',location = 'northwest');

xlabel('Time (s)');
ylabel('Displacement error (m)');
grid on;
hold off

saveas(Euler_time_step_study_2, 'Euler_time_step_study_2.svg');

%% 3c: step size study


RK4_time_step_study = figure;
%set(RK4_time_step_study,"WindowState","maximized");
set(findall(RK4_time_step_study,'-property','FontSize'),'FontSize',24);
set(findall(RK4_time_step_study,'-property','Interpreter'),'Interpreter','latex') 
set(findall(RK4_time_step_study,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(RK4_time_step_study,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(RK4_time_step_study,'Position');
set(RK4_time_step_study,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);



semilogy(time__step(1,:),RK4_step_study(1,:),"r--");
hold on;
semilogy(time__step(2,:),RK4_step_study(2,:),"b-.");
semilogy(time__step(3,:),RK4_step_study(3,:),"g:");
legend('\Deltat =0.5','\Deltat =0.05','\Deltat =0.005',location = 'south');

xlabel('Time (s)');
ylabel('Displacement error (m)');
grid on;
hold off

saveas(RK4_time_step_study, 'RK4_time_step_study.svg');


%% initial conditions effects


dt = 0.005;
Time__range = 0:0.005:10;
init_conditions = [5,4,2,0];
Y__Eular = zeros(length(init_conditions),length(Time__range));
Y__RK4 = zeros(length(init_conditions),length(Time__range));
analytical_res = zeros(length(init_conditions),length(Time__range));

itter = 1;
for i = init_conditions

    Y__0 = [0.1;i];
    temp_euler = Explicit_Euler(Y__0, A, dt, Time__range);
    Y__Eular(itter,:) = temp_euler(1,:);
    temp_Rk4 = RK4(Y__0, A, dt, Time__range);
    Y__RK4(itter,:) = temp_Rk4(1,:);
    

    analytical_res(itter, :) = 0.1 * cos(omega .* Time__range) + (i / omega) * sin(omega .* Time__range); 
    itter = itter + 1;

end

%%

initil_conditions_Euler = figure;

set(findall(initil_conditions_Euler,'-property','FontSize'),'FontSize',24);
set(findall(initil_conditions_Euler,'-property','Interpreter'),'Interpreter','latex') 
set(findall(initil_conditions_Euler,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(initil_conditions_Euler,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(initil_conditions_Euler,'Position');
set(initil_conditions_Euler,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);

Line_Style = ["-","-.","--",":",":"];

for i = 1:1:length(init_conditions)

    plot(Time__range,Y__Eular(i,:),'DisplayName', ['v_0=' num2str(init_conditions(i)) ','], LineStyle=Line_Style(i));

    hold on;

end  


T_____RANGE = 0:2:10;
y_sub = -1*ones(length(T_____RANGE),1);

plot(T_____RANGE, y_sub,'DisplayName', 'Submerged',LineWidth=2, Marker='x',Color="k")
hold off

legend;
xlabel('Time, (s)');
ylabel('Displacement, y (m)');
grid on;
legend(Orientation="horizontal",Location="northoutside");
hold off

saveas(initil_conditions_Euler, 'initil_conditions_Euler.svg');


%%

initil_conditions_RK4 = figure;

set(findall(initil_conditions_RK4,'-property','FontSize'),'FontSize',24);
set(findall(initil_conditions_RK4,'-property','Interpreter'),'Interpreter','latex') 
set(findall(initil_conditions_RK4,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(initil_conditions_RK4,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(initil_conditions_RK4,'Position');
set(initil_conditions_RK4,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);



for i = 1:1:length(init_conditions)

    plot(Time__range,Y__RK4(i,:),'DisplayName', ['v_0=' num2str(init_conditions(i)) ','], LineStyle=Line_Style(i));

    hold on;

end  

T_____RANGE = 0:2:10;
y_sub = -1*ones(length(T_____RANGE),1);

plot(T_____RANGE, y_sub,'DisplayName', 'Submerged',LineWidth=2, Marker='x',Color="k")
hold off


legend;
xlabel('Time, (s)');
ylabel('Displacement, y (m)');
grid on;
legend(Orientation="horizontal",Location="northoutside")
hold off
saveas(initil_conditions_RK4, 'initil_conditions_RK4.svg');


%%
%analytical initial conditions

initil_conditions_analytical = figure;

set(findall(initil_conditions_analytical,'-property','FontSize'),'FontSize',24);
set(findall(initil_conditions_analytical,'-property','Interpreter'),'Interpreter','latex') 
set(findall(initil_conditions_analytical,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(initil_conditions_analytical,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(initil_conditions_analytical,'Position');
set(initil_conditions_analytical,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)]);



for i = 1:1:length(init_conditions)

    plot(Time__range,analytical_res(i,:),'DisplayName', ['v_0=' num2str(init_conditions(i)) ','], LineStyle=Line_Style(i));

    hold on;

end  

T_____RANGE = 0:2:10;
y_sub = -1*ones(length(T_____RANGE),1);

plot(T_____RANGE, y_sub,'DisplayName', 'Submerged',LineWidth=2, Marker='x',Color="k")
hold off

legend;
xlabel('Time, (s)');
ylabel('Displacement, y (m)');
grid on;
legend(Orientation="horizontal",Location="northoutside");
hold off

saveas(initil_conditions_analytical, 'initil_conditions_analytical.svg');


%% Q4: Damping Effects

c = 0.2;

A = [0, 1; -omega_2, 0];

Y_RKA_ANALYTICAL = RK4(Y_0, A, delta_t, time_range);
A_damped = [0, 1; -omega_2, -c];

Y_RK4_damped = RK4(Y_0, A_damped, delta_t, time_range);

c_crit = 2*omega;

c_over = 2*c_crit;

A_crit = [0, 1; -omega_2, -c_crit];
Y_RK4_damped_crit = RK4(Y_0, A_crit, delta_t, time_range);

A_over = [0, 1; -omega_2, -c_over];
Y_RK4_damped_over = RK4(Y_0, A_over, delta_t, time_range);


Runge_Kutta_4_damped = figure;
%set(Runge_Kutta_4_damped,"WindowState","maximized");
set(findall(Runge_Kutta_4_damped,'-property','FontSize'),'FontSize',24);
set(findall(Runge_Kutta_4_damped,'-property','Interpreter'),'Interpreter','latex') 
set(findall(Runge_Kutta_4_damped,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
set(Runge_Kutta_4_damped,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
pos = get(Runge_Kutta_4_damped,'Position');
set(Runge_Kutta_4_damped,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])

plot(time_range,Y_RKA_ANALYTICAL(1,:),"k-");
hold on;
plot(time_range,Y_RK4_damped(1, :),"b-.");
plot(time_range,Y_RK4_damped_crit(1,:),"r-- ");
plot(time_range,Y_RK4_damped_over(1,:),"g: ");

%title("Figure 2: Runge-Kutta 4th Order, Damped Case: Displacement vs Time")
xlabel("Time, (s)")
ylabel("Displacement, (m)")
legend("c=  0.0","c= 0.2", "c= 2\omega", "c= 4\omega", Location = "northoutside",Orientation="horizontal");
grid minor
hold off



saveas(Runge_Kutta_4_damped,'Runge_Kutta_4_damped.svg');



%% Q4: Damping Effects
% 
% c = 0.2;
% 
% A = [0, 1; -omega_2, 0];
% 
% 
% 
% A_damped = [0, 1; -omega_2, -c];
% 
% Y_RK4_damped = RK4(Y_0, A_damped, delta_t, time_range);
% 
% c_crit = 2*omega;
% 
% c_over = 2*c_crit;
% 
% A_crit = [0, 1; -omega_2, -c_crit];
% Y_RK4_damped_crit = RK4(Y_0, A_crit, delta_t, time_range);
% 
% A_over = [0, 1; -omega_2, -c_over];
% Y_RK4_damped_over = RK4(Y_0, A_over, delta_t, time_range);
% 
% 
% Runge_Kutta_4_damped = figure;
% %set(Runge_Kutta_4_damped,"WindowState","maximized");
% set(findall(Runge_Kutta_4_damped,'-property','FontSize'),'FontSize',24);
% set(findall(Runge_Kutta_4_damped,'-property','Interpreter'),'Interpreter','latex') 
% set(findall(Runge_Kutta_4_damped,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% set(Runge_Kutta_4_damped,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% pos = get(Runge_Kutta_4_damped,'Position');
% set(Runge_Kutta_4_damped,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% 
% plot(time_range,Y_RKA_ANALYTICAL(1,:),"k-");
% hold on;
% plot(time_range,Y_RK4_damped(1, :),"b-.");
% plot(time_range,Y_RK4_damped_crit(1,:),"r-- ");
% plot(time_range,Y_RK4_damped_over(1,:),"g: ");
% 
% %title("Figure 2: Runge-Kutta 4th Order, Damped Case: Displacement vs Time")
% xlabel("Time, (s)")
% ylabel("Displacement, (m)")
% legend("c=  0.0","c= 0.2", "c= 2\omega", "c= 4\omega", Location = "northoutside",Orientation="horizontal");
% grid minor
% hold off
% 
% 
% 
% saveas(Runge_Kutta_4_damped,'Runge_Kutta_4_damped.svg');
% 
% 
% 



%%

% % hfig = figure;  % save the figure handle in a variable
% % t = 0:0.02:10; x = t.*sin(2*pi*t)+ 2*rand(1,length(t)); % data
% % plot(t,x,'k-','LineWidth',1.5,'DisplayName','$\Omega(t)$');
% % xlabel('time $t$ (s)')
% % ylabel('$\Omega$ (V)')
% % fname = 'myfigure';
% % 
% % picturewidth = 20; % set this parameter and keep it forever
% % hw_ratio = 0.65; % feel free to play with this ratio
% % set(findall(hfig,'-property','FontSize'),'FontSize',17) % adjust fontsize to your document
% % 
% % set(findall(hfig,'-property','Box'),'Box','off') % optional
% % set(findall(hfig,'-property','Interpreter'),'Interpreter','latex') 
% % set(findall(hfig,'-property','TickLabelInterpreter'),'TickLabelInterpreter','latex')
% % set(hfig,'Units','centimeters','Position',[3 3 picturewidth hw_ratio*picturewidth])
% % pos = get(hfig,'Position');
% % set(hfig,'PaperPositionMode','Auto','PaperUnits','centimeters','PaperSize',[pos(3), pos(4)])
% % %print(hfig,fname,'-dpdf','-painters','-fillpage')
% % print(hfig,fname,'-dpng','-vector')
% % 
% % saveas(hfig,"hfig.svg");

%% Runge-Kutta 4th order Function

function Y_RK4 = RK4(Y_0, A, delta_t, time_range)

    % Prealocating Array
    Y_RK4 = zeros(2,length(time_range));
    Y_RK4(:,1) = Y_0;
    
    % RK4 intergration steps
    for i = 2:length(time_range)
        k1 = A * Y_RK4(:, i-1);
        k2 = A * (Y_RK4(:, i-1) + k1 * delta_t/2);
        k3 = A * (Y_RK4(:, i-1) + k2 * delta_t/2);
        k4 = A * (Y_RK4(:, i-1) + k3 * delta_t);
        Y_RK4(:, i) = Y_RK4(:, i-1) + delta_t * (k1 + 2*k2 + 2*k3 + k4)/6;
    end
end
    
%% Explicit Euler Function

function Y_Euler = Explicit_Euler(Y_0, A, delta_t, time_range)
    
    % Prealocating Array
    Y_Euler = zeros(2, length(time_range));
    Y_Euler(:, 1) = Y_0;
    
    % Euler Intergration steps
    for i = 1:(length(time_range) - 1)
        Y_Dot = A * Y_Euler(:, i);
        Y_Euler(:, i+1) = Y_Euler(:, i) + delta_t * Y_Dot;
    end
end





