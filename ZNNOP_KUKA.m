% This is the code for ZNNOP for tracking cardiod by KUKA
%clc
clear
format long

%% initialize parameters
dt = 0.001 ; %time instant
t_end = 10; %end time for path
t_span = 0:dt:t_end; %time of simulation
t_sim = t_end; % end time for simulation

m = 3; %no. of equations
n = 7; % no of parametes

% performance constraints preliminaries
rho_inf = 0.005; %final value of bound

delta = 100; % exponential decay constant
gamma  = 100; %ZNN gain
k_pos = 2; % postional error gain
f = 2; % parameeter for initial values of PPC

%% Initial Conditions
%th(:,1) = [pi/4,-pi/2,-pi/4,pi/6,pi/3,pi/2,pi/6];
th(:,1) = [-0.2578,0.2946,0.5946,-0.9896,0.4993,0.3100,0.8620];
%example 2
%th(:,1) = [pi/14,pi/14,pi/14,pi/14,pi/14,pi/14];
th_dot = 0*ones(n,1);
lambda_0 = [1;1;1]; % initiual values of the lagragian coeffecients
y(:,1) = [th_dot;lambda_0]; % cascaded vector

%% parameters of KUKA robot
a8 = 0.126; %length of the end effector

d1 = 0.36;
d3 = 0.42;
d5 = 0.4;
 
%% position with respect to the body frame
% actual position of the robot
r_a(:,1) = forward_map_kuka(a8,d1,d3,d5,th(:,1)); % ac
% deviation of desired path from actual path
% initially
dev = [-0.04;0.035;0.025];
r_d(:,1) = r_a(:,1)+dev;  %initial position of end effector

%% constraints
%Constraints
% upper limit of joint velocities
th_dot_pos = 3.6*[0.0001;0.1;0.1;0.1;0.0001;0.0001;0.1]; 
% lower limit of joint velocities
th_dot_neg = -th_dot_pos;

string = '3rd + 5th + 6th Joint Fault';

%% Initializing Matricies
th_dot_pos_plot = th_dot_pos*ones(1,size(t_span,2)); % Upper bound for plot
th_dot_neg_plot = th_dot_neg*ones(1,size(t_span,2)); % lower boud for plot

error = zeros(m,size(t_span,2));

%% ZNN path tracking
for i=1:length(t_span)
   t = t_span(i);
   
    % Kinematics of the Robot
    J = Jacobian_matrix_kuka(a8,d1,d3,d5,th(:,i)); % Jacobian Matrix
    J_dot = Jacobian_dot_kuka(a8,d1,d3,d5,th(:,i),th_dot(:,i)); % Jacobian dot
    

    % desired path parameters
    c = 0.03;
    tilt = pi/6;
    shape = 'tricuspid';
    [r_d(:,i+1), r_d_dot ,r_d_ddot] = path_3D(c,tilt,r_d(:,1),t,t_end,shape);
  
    %error for drift free case
    error(:,i) = J*y(1:n,i)-r_d_dot + k_pos*(r_a(:,i)-r_d(:,i));           %error_function e(t)
    diff_error_t(:,i) = J_dot*y(1:n,i)-r_d_ddot + k_pos*(J*y(1:n,i)-r_d_dot);     %d(e(t))/dt
    
    %performance constraint
    rho_0 = f*ones(m,1);%error(:,1);
    rho_r = abs(rho_0)*exp(-delta*t)+rho_inf; 
    diff_rho_r_t = -delta*abs(rho_0)*exp(-delta*t);  %d(rho_r)/dt
    rho_l = -rho_r;
    diff_rho_l_t = -diff_rho_r_t;
    
    %neta transformations
    neta_1 = (y(1:n,i)-th_dot_neg)./(th_dot_pos-th_dot_neg);
    neta_2 = (error(:,i)-rho_l)./(rho_r-rho_l);
    
    %% chi transformations and g
    chi_1 = log(neta_1./(1-neta_1));
    chi_2 = log(neta_2./(1-neta_2));
    g = (chi_2.*chi_2)/2;
    
    % time and spartaial differentiation
    for j = 1:m
        diff_neta_2j_neta = neta_2(j)*(1-neta_2(j));
        diff_sec_2_neta_2 = (1+(2*neta_2(j)-1)*chi_2(j))/(diff_neta_2j_neta^2);
        diff_neta_2_t = ((rho_r(j)-rho_l(j))*(diff_error_t(j,i)-diff_rho_l_t(j))-(diff_rho_r_t(j)-diff_rho_l_t(j))*(error(j,i)-rho_l(j)))/((rho_r(j)-rho_l(j))^2);
        for k = 1:n
            diff_g_x(j,k)= J(j,k)*chi_2(j)/((rho_r(j)-rho_l(j))*(diff_neta_2j_neta));
            diff_sec_1_t = ((J_dot(j,k)*(rho_r(j)-rho_l(j)))-(J(j,k)*(diff_rho_r_t(j)-diff_rho_l_t(j))))/((rho_r(j)-rho_l(j))^2);
            diff_g_x_t(j,k) = (J(j,k)/(rho_r(j)-rho_l(j)))*diff_sec_2_neta_2*diff_neta_2_t+(chi_2(j)/diff_neta_2j_neta)*diff_sec_1_t;
        end
        diff_g_t(j,1) = chi_2(j)*diff_neta_2_t/(diff_neta_2j_neta);
    end
    for p = 1:n
        diff_neta_1p_neta = neta_1(p)*(1-neta_1(p));
        diff_delta_1_x(p,1) = chi_1(p)/((th_dot_pos(p)-th_dot_neg(p))*diff_neta_1p_neta); %d(delta_1)/dx ;;delta_1 = (1/2)*(chi_1)^T*chi_1
         for q =1:n
             if p==q
             diff_delta_1_x_x(p,q) = (1+(2*neta_1(p)-1)*chi_1(p))/(((th_dot_pos(p)-th_dot_neg(p))*diff_neta_1p_neta)^2);
            else
             diff_delta_1_x_x(p,q) = 0;
             end
            sum_glam_x = 0;
            for s = 1:m
             diff_chi_2s_neta = neta_2(s)*(1-neta_2(s));
             ele_glam_x = y(n+s,i)*J(s,p)*J(s,q)*(1+(2*neta_2(s)-1)*chi_2(s))/(((rho_r(s)-rho_l(s))*diff_chi_2s_neta)^2);
             sum_glam_x = sum_glam_x+ele_glam_x;
            end
            diff_glam_x(p,q) = sum_glam_x;
         end
    end
    h = [diff_delta_1_x+(diff_g_x'*y(n+1:end,i));
         g];
    diff_h_t = [diff_g_x_t'*y(n+1:end,i);
               diff_g_t];
    diff_h_y = [(diff_delta_1_x_x+diff_glam_x),  diff_g_x'
                diff_g_x, zeros(m,m)];
    
    y_dot(:,i) = -pinv(diff_h_y)*(gamma*(h)+diff_h_t);
      
    if i<=3 %euler
         y(:,i+1) = y(:,i) + y_dot(:,i)*dt;
        th_dot(:,i+1)= y(1:n,i+1);
        th(:,i+1) = th(:,i) + th_dot(:,i)*dt; 
    else %FIFD
        y(:,i+1) = (8/5)*y_dot(:,i)*dt + (3/5)*y(:,i) + (1/5)*y(:,i-1)+ (1/5)*y(:,i-2);
        th_dot(:,i+1)= y(1:n,i+1);
        th(:,i+1) = (8/5)*th_dot(:,i)*dt + (3/5)*th(:,i) + (1/5)*th(:,i-1)+ (1/5)*th(:,i-2); %FIFD
    end

    % errors
    error_norm(:,i) = norm(error(:,i)); % norm of error function    
    error_stat(:,i) = (J*y(1:n,i) - r_d_dot); % velocity level error
    error_trag(:,i) = r_a(:,i) - r_d(:,i); % end effector position error
    r_a(:,i+1) = forward_map_kuka(a8,d1,d3,d5,th(:,i+1));
    control_eff(:,i) = norm(th_dot(:,i)); %% control effort
    error_trag_norm(:,i) = norm(error_trag(:,i));
end

%% plots
figure(1)
plot(t_span,th_dot(1,1:end-1),t_span,th_dot(2,1:end-1),'c--',t_span,th_dot(3,1:end-1),'y:',t_span,th_dot(4,1:end-1),'g-.',t_span,th_dot(5,1:end-1),'m',t_span,th_dot(6,1:end-1),'k--',t_span,th_dot(7,1:end-1), 'LineWidth',1.5)
%title('Variation of the states of system with time')
hold on
plot(t_span,th_dot_pos_plot(2,:),'r--',t_span,th_dot_neg_plot(2,:),'r--')
hold on
h2 =legend('$\dot{\theta}_1 $','$\dot{\theta}_2$','$\dot{\theta}_3$','$\dot{\theta}_4$','$\dot{\theta}_5$','$\dot{\theta}_6$','$\dot{\theta}_7$','$\dot{\theta}^+ $','$\dot{\theta}^-$',[400 260 0 0]);
set(h2,'Interpreter', 'latex');
xlabel('t(s)')
hl = ylabel('$\dot{\theta}_{i}$');
set(hl,'Interpreter', 'latex');
ylim([th_dot_neg(2)+0.05*th_dot_neg(2) th_dot_pos(2)+0.05*th_dot_pos(2) ])

figure(2)
plot(t_span,th(1,1:end-1),t_span,th(2,1:end-1),'c--',t_span,th(3,1:end-1),'y:',t_span,th(4,1:end-1),'g-.',t_span,th(5,1:end-1),'m',t_span,th(6,1:end-1),'k--',t_span,th(7,1:end-1),'--', 'LineWidth',1.5)
%title('Variation of the states of system with time')
h2 =legend('${\theta}_1 $','${\theta}_2$','${\theta}_3$','${\theta}_4$','${\theta}_5$','${\theta}_6$','${\theta}_7$');
set(h2,'Interpreter', 'latex');
xlabel('t(s)')
hold on

% 
figure(3)
plot(t_span(1:(t_end/dt)+1),error_trag(1,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),error_trag(2,1:(t_end/dt)+1),'--',t_span(1:(t_end/dt)+1),error_trag(3,1:(t_end/dt)+1),'-.');
%title('Variation of Position tracking error with time')
xlabel('t(s)')
ylabel('Position tracking error(\epsilon) (m)')
legend('\epsilon_x','\epsilon_y','\epsilon_z')
hold on

figure(4)
plot(t_span(1:(t_end/dt)+1),error_stat(1,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),error_stat(2,1:(t_end/dt)+1),t_span(1:(t_end/dt)+1),error_stat(3,1:(t_end/dt)+1))
%title('Variation of velocity tracking error with time')
xlabel('t(s)')
ylabel('Velocity tracking error (m/s)')
legend('\epsilon_1','\epsilon_2','\epsilon_3')
hold on

figure(5)
plot(t_span(1:(t_sim/dt)+1),error_norm(1:(t_sim/dt)+1),'linewidth',1.2);
title('Variation of  ||e(t)||_2 with time ')
xlabel('t(s)')
ylabel('||e(t)||_2')
hold on


figure(6)
% simple tragectory
plot3(r_d(1,:),r_d(2,:),r_d(3,:),'r','linewidth',1.5)
hold on
plot3(r_a(1,:),r_a(2,:),r_a(3,:),'b--','linewidth',1.5)
hold on
xlabel('x(m)')
ylabel('y(m)')
zlabel('z(m)')
legend('Traced Tragectory','Desired Tragectory')

%% control effort
figure(7)
plot(t_span,control_eff)
cont_eff_total = norm(control_eff);
hold on


figure(8)
plot(t_span,y_dot(1,:),t_span,y_dot(2,:),'c--',t_span,y_dot(3,:),'y:',t_span,y_dot(4,:),'g-.',t_span,y_dot(5,:),'m',t_span,y_dot(6,:),'k--',t_span,y_dot(7,:),'--')
%title('Variation of the states of system with time')
h2 =legend('$\ddot{\theta}_1 $','$\ddot{\theta}_2$','$\ddot{\theta}_3$','$\ddot{\theta}_4$','$\ddot{\theta}_5$','$\ddot{\theta}_6$','$\ddot{\theta}_7$');
set(h2,'Interpreter', 'latex');
xlabel('t(s)')
hold on

figure(9)
plot(t_span(1:(t_sim/dt)+1),error_trag_norm(1:(t_sim/dt)+1),'linewidth',1.2);
title('Variation of  ||e(t)||_2 with time ')
xlabel('t(s)')
ylabel('||e(t)||_2')
hold on

figure(10)
semilogy(t_span(1:(t_sim/dt)+1),error_trag_norm(1:(t_sim/dt)+1),'linewidth',1.2);
title('Variation of  ||e(t)||_2 with time ')
xlabel('t(s)')
ylabel('||e(t)||_2')
hold on


% performance parameters
normal_control_eff = norm(control_eff)*sqrt(dt/(n*t_end))
mean_trac_err = norm(error_trag_norm)*sqrt(dt/(m*t_end))

% initial deviation percentages
x_per = dev(1)/(max(r_a(1,:))-min(r_a(1,:)))*100
y_per = dev(2)/(max(r_a(2,:))-min(r_a(2,:)))*100
z_per = dev(3)/(max(r_a(3,:))-min(r_a(3,:)))*100

