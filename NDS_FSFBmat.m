%% Nathalia De Souza -- AERO 560 -- HW 2

clear 
close all
clc

%% GEOCANSAT COEs and Mass Properties

% size
h = 10;         % height [m]
r = 1.5;        % radius [m]
mass = 1200;    % mass in [kg]

% moment of inertia matrix in body frame [kg * m^2]
Ixx = (1/12) * mass * (h^2 + (3 * r^2));
Iyy = Ixx;
Izz = (1/2) * mass * r^2;

J = [Ixx 0 0 ; 0 Iyy 0; 0 0 Izz];   % MOI matrix in body frame


% initial torques
T_d0 = zeros(3,1);
T_c0 = zeros(3,1);


% command (desired/target) quaternion
q_LVLH_b_c = [0 0 0 1]';

% COEs
mu = 398600;    % Earth standard gravitational parameter [km^3/s^2]
H = 129654.3;   % specific angular momentum [km/s^2]
e = 0;          % eccentricity [deg]
RAAN = 0;       % right ascension of ascending node [deg] (usually denoted by cap omega)
inc = 0;        % inclination [deg] (usually denoted by i)
AOP = 0;        % argument of perigee [deg] (usually denoted as little omega)
nu = 0;         % true anomaly [deg] (usually denoted by either nu or theta)


sim_time = 100; % simulation time [s]

%% Settling time, Damping, and Gains

DCR = 0.65; % damping coefficient ratio
ts = 30;    % 2% settling time [s]
wn = (-log(0.02 * sqrt(1 - DCR^2)))/DCR/(ts); % undamped natural frequency [rad/s]

w_n = diag([wn wn wn]);

kd = 2 .* J * DCR .* w_n;   % derivative gain
kp = 2 .* J .* (w_n).^2;    % proportional gain

%% Initial Conditions

% Note: I follow this reference frame naming convention: a transformation 
% from frame1 to frame2 is labeled x_frame2_frame1

eps_LVLH_b0 = [0.5; 0.5; 0.5];              % vector part of the initial quaternion from body to LVLH
eta_LVLH_b0 = 0.5;                          % scalar part of the initial quaternion from body to LVLH
w_LVLH_b0_b = (10^-5) .* [0.5; -7.27; 3.0]; % initial angular velocity from LVLH to body in body components

%% Frame Rotations

% initial rotation from ECI to LVLH 
[Q_peri2ECI] = EulerRotation313(RAAN, inc, AOP);
[rECI0, vECI0] = ECIstate(H, e, nu, Q_peri2ECI);

n = norm(vECI0)/norm(rECI0); % mean motion

% LVLH vectors in ECI components
zLVLH__ECI = -rECI0/norm(rECI0);
hLVLH__ECI = cross(rECI0, vECI0);
yLVLH__ECI = -(hLVLH__ECI/norm(hLVLH__ECI));
xLVLH__ECI = cross(yLVLH__ECI, zLVLH__ECI);


% rotation from ECI to LVLH
C_LVLH_ECI0 = [xLVLH__ECI'; yLVLH__ECI'; zLVLH__ECI'];
q_LVLH_ECI0 = C2quat(C_LVLH_ECI0);


% rotation from LVLH to body
q_b_LVLH0 = [eps_LVLH_b0; eta_LVLH_b0];
C_b_LVLH0 = quat2C(q_b_LVLH0);


% rotation from ECI to body
q_b_ECI0 = quatmult(q_LVLH_ECI0, q_b_LVLH0);
C_b_ECI0 = quat2C(q_b_ECI0);


% initial euler angles
E_LVLH_ECI0 = C2EA(C_LVLH_ECI0);        % intial euler angles from ECI to LVLH
E_b_LVLH0 = quat2EA(q_b_LVLH0);       % intial euler angles from LVLH to body
E_b_ECI0 = quat2EA(q_b_ECI0);           % intial euler angles from ECI to body


% initial angular velocities [rad/s]
w_ECI_b0_b = w_LVLH_b0_b;                     % initial rotation about ECI from body in body coordinates
w_LVLH0_b_b = hLVLH__ECI/((norm(rECI0))^2);   % rotation about LVLH from body in body coordinates
w_LVLH_ECI0_LVLH = C_LVLH_ECI0 * w_LVLH0_b_b; % rotation about LVLH from ECI in LVLH coordinates
w_b_LVLH0 = w_ECI_b0_b - w_LVLH_ECI0_LVLH;    % rotation about body from LVLH


%% Plots

out = sim('NDS_FSFB.slx');
time = out.tout;

figure(1)
sgtitle('Cubesat Reference Mission Dynamics', 'interpreter', 'latex')

% my naming convention is accidentally flipped for e to b

% body to ECI
w_ECI_b = squeeze(out.w_ECI_b.Data);
q_ECI_b = squeeze(out.q_ECI_b.Data);
E_ECI_b = squeeze(out.E_ECI_b.Data);

subplot(3,3,9)
plot(time, w_ECI_b, 'LineWidth', 1.25)
grid on
title('$\omega_{b->ECI}$ in $\mathcal{F}_b$ components', 'interpreter', 'latex')
legend('\omega_x', '\omega_y', '\omega_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Angular Velocity [rad/s]', 'interpreter', 'latex')


subplot(3,3,3)
plot(time, q_ECI_b, 'LineWidth', 1.25)
grid on
title('Quaternion from $\mathcal{F}_b$ to $\mathcal{F}_{ECI}$', 'interpreter', 'latex')
legend('\epsilon_x', '\epsilon_y', '\epsilon_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Quaternion Component', 'interpreter', 'latex')


subplot(3,3,6)
plot(time, E_ECI_b, 'LineWidth', 1.25)
grid on
legend('\phi_x', '\theta_y', '\psi_z')
title('Euler Angles from $\mathcal{F}_b$ to $\mathcal{F}_{ECI}$', 'interpreter', 'latex')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Euler Angles [rad]', 'interpreter', 'latex')



% LVLH to ECI in LVLH frame
w_LVLH_ECI_LVLH = squeeze(out.w_LVLH_ECI_LVLH);
q_ECI_LVLH_LVLH = squeeze(out.q_ECI_LVLH_LVLH);
E_ECI_LVLH_LVLH = squeeze(out.E_ECI_LVLH_LVLH);


subplot(3,3,8)
plot(w_LVLH_ECI_LVLH, 'LineWidth', 1.25)
grid on
title('$\omega_{LVLH->ECI}$ in $\mathcal{F}_{LVLH}$ components', 'interpreter', 'latex')
legend('\omega_x', '\omega_y', '\omega_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Angular Velocity [rad/s]', 'interpreter', 'latex')


subplot(3,3,2)
plot(q_ECI_LVLH_LVLH, 'LineWidth', 1.25)
grid on
title('Quaternion from $\mathcal{F}_{LVLH}$ to $\mathcal{F}_{ECI}$', 'interpreter', 'latex')
legend('\epsilon_x', '\epsilon_y', '\epsilon_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Quaternion Component', 'interpreter', 'latex')


subplot(3,3,5)
plot(E_ECI_LVLH_LVLH, 'LineWidth', 1.25)
grid on
legend('\phi_x', '\theta_y', '\psi_z')
title('Euler Angles from $\mathcal{F}_{LVLH}$ to $\mathcal{F}_{ECI}$', 'interpreter', 'latex')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Euler Angles [rad]', 'interpreter', 'latex')



% body to LVLH
w_LVLH_b = squeeze(out.w_LVLH_b);
q_LVLH_b = squeeze(out.q_LVLH_b);
E_LVLH_b = squeeze(out.E_LVLH_b);


subplot(3,3,7)
plot(w_LVLH_b, 'LineWidth', 1.25)
grid on
title('$\omega_{b->LVLH}$ in $\mathcal{F}_b$ components', 'interpreter', 'latex')
legend('\omega_x', '\omega_y', '\omega_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Angular Velocity [rad/s]', 'interpreter', 'latex')


subplot(3,3,1)
plot(q_LVLH_b, 'LineWidth', 1.25)
grid on
title('Quaternion from $\mathcal{F}_{b}$ to $\mathcal{F}_{LVLH}$', 'interpreter', 'latex')
legend('\epsilon_x', '\epsilon_y', '\epsilon_z')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Quaternion Component', 'interpreter', 'latex')


subplot(3,3,4) 
plot(E_LVLH_b, 'LineWidth', 1.25)
grid on
legend('\phi_x', '\theta_y', '\psi_z')
title('Euler Angles from $\mathcal{F}_b$ to $\mathcal{F}_{LVLH}$', 'interpreter', 'latex')
xlabel('Time [s]', 'interpreter', 'latex')
ylabel('Euler Angles [rad]', 'interpreter', 'latex')


%% Control Torque
Td = squeeze(out.Td.data);
Tc = squeeze(out.Tc.data);

figure(2)
plot(out.Tc, 'LineWidth', 1.25)
grid on
title("Tc in $\mathcal{F}_{ECI}$ components", 'interpreter', 'latex')
xlabel("Time [s]", 'interpreter', 'latex')
ylabel("Control torque (N * m)", 'interpreter', 'latex')
legend('T_{cx}', 'T_{cy}', 'T_{cz}')


%% Nadir Angle Calculation

rECI = squeeze(out.rECI.data);
vECI = squeeze(out.vECI.data);


nadirangle = zeros();
for i = 1:length(rECI)
    z_LVLH = -rECI(i,:)/norm(rECI(i,:));
    C_b_ECI = quat2C(q_ECI_b(:,i));
    C_ECI_b = C_b_ECI';
    zb = C_ECI_b(:,3);
    nadirangle(i,1) = acos(dot(zb, z_LVLH));

end


figure(3)
plot(linspace(1,100,length(rECI)), nadirangle, 'LineWidth', 1.25)
grid on
xlabel('Time [s]','interpreter', 'latex')
ylabel('Nadir Angle [rad/s]','interpreter', 'latex')
title('Body z-axis to Nadir vs Time','interpreter', 'latex')
ylim([-0.1 1.7])
