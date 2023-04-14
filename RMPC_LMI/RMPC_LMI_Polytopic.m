clc
clear
close all

kappa = .787;
A = @(alpha_k)[1 0.1;
      0 1 - 0.1 * alpha_k];

B = [0;
     0.1*kappa];

C = [1 0];

%% Figure 4a: Using nominal MPC with alpha_k = 1 sec^{-1};
A_nominal = A(1);
nx = size(A_nominal, 1);
nu = size(B, 2);
x0 = [0.05;0];
Q1 = C'*C;
R = .00002;

% Calculating F = YQ^-1
gamma = sdpvar(1, 1);
Q = sdpvar(nx, nx);
Y = sdpvar(nu, nx);


LMI1 = [1 x0';
        x0 Q];


LMI2 = [Q                       (A_nominal * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
        A_nominal * Q + B * Y   Q                           zeros(nx, nx)       zeros(nx, nu);
        sqrt(Q1)*Q              zeros(nx, nx)               gamma * eye(nx, nx) zeros(nx, nu);
        sqrt(R)*Y              zeros(nu, nx)               zeros(nu, nx)       gamma * eye(nu, nu)];
ZEROS = 1e-7;

Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, Q >= ZEROS];
Objective = gamma;
solvesdp(Constraints, Objective);
Y = double(Y);
Q = double(Q);

F_nor = Y / Q;
X_nor = [x0];
x_next = x0;
A_real = A(9);
for i = 1 : 40
    x_next = A_real * x_next + B * F_nor * x_next;
    X_nor = [X_nor x_next];
end

figure
plot(0:0.1:4,X_nor(1,:))
hold on;
plot(0:0.1:4,X_nor(2,:));
axis([0 4, -0.4 0.3]);
xlabel('time(sec)');
ylabel('$\theta (rad) and  \dot{\theta}(rad/sec)$','interpreter','latex');
title('Using nominal MPC with \alpha(k)=1 sec^{-1}');
h1 = legend('$\theta$','$\dot{\theta}$');
set(h1,'interpreter','latex');

%% Figure 4b Using robust LMI-based MPC 

A1 = A(0.1);
A2 = A(10);
nx = size(A1, 1);
nu = size(B, 2);
x0 = [0.05;0];
Q1 = C'*C;
R = .00002;

gamma = sdpvar(1, 1);
Q = sdpvar(nx, nx);
Y = sdpvar(nu, nx);


LMI1 = [1 x0';
        x0 Q];


LMI2 = [Q                       (A1 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
        A1 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
        sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
        sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];

LMI3 = [Q                       (A2 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
        A2 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
        sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
        sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];

ZEROS = 1e-7;

Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, Q >= ZEROS];
Objective = gamma;
solvesdp(Constraints, Objective);
Y = double(Y);
Q = double(Q);

F_robust = Y / Q;

X_robust = [x0];
x_next = x0;
A_real = A(9);
for i = 1 : 40
    x_next = A_real * x_next + B * F_robust * x_next;
    X_robust = [X_robust x_next];
end

figure
plot(0:0.1:4,X_robust(1,:))
hold on;
plot(0:0.1:4,X_robust(2,:));
axis([0 4,-0.35 0.1]);
xlabel('time(sec)');
ylabel('$\theta (rad) and  \dot{\theta}(rad/sec)$','interpreter','latex');
title('Using robust LMI-based MPC ');
h2 = legend('$\theta$','$\dot{\theta}$');
set(h2,'interpreter','latex');


%% Figure 5 Closed-loop responses for the time-varying system with input constraint; 

A1 = A(0.1);
A2 = A(10);
nx = size(A1, 1);
nu = size(B, 2);
x0 = [1;0];
Q1 = C'*C;
R = .00002;

% Robust receding horizon state-feedback

X_Receding_static = [x0];
x_next = x0;

gamma = sdpvar(1, 1);
X = sdpvar(nu, nu);
Q = sdpvar(nx, nx);
Y = sdpvar(nu, nx);
LMI1 = [1 x_next';
        x_next Q];


LMI2 = [Q                       (A1 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
        A1 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
        sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
        sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];

LMI3 = [Q                       (A2 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
        A2 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
        sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
        sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];

% Input constraints  |uk| <= 2
LMI4 = [X  Y;
        Y' Q];

ZEROS = 1e-7;
Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, LMI4 >= ZEROS, Q >= ZEROS];
for j = 1 : nu
    Constraints = [Constraints X(j, j) <= 4] ;
end

Objective = gamma;
solvesdp(Constraints, Objective)
Y = double(Y);
Q = double(Q);

F_Receding_Static = Y / Q;
U_Receding_static = [F_Receding_Static * x0];
for i = 1 : 100
    alpha = (10 - 0.1) * rand + 0.1;
    A_Receding = A(alpha);
    U_Receding_static = [U_Receding_static   F_Receding_Static * x_next];
    x_next = A_Receding * x_next + B * F_Receding_Static * x_next;
    X_Receding_static = [X_Receding_static x_next];
end






X_Receding = [x0];
x_next = x0;
U_Receding = [F_Receding_Static * x0];

for i = 1 : 100
    gamma = sdpvar(1, 1);
    X = sdpvar(nu, nu);
    Q = sdpvar(nx, nx);
    Y = sdpvar(nu, nx);
    LMI1 = [1 x_next';
            x_next Q];


    LMI2 = [Q                       (A1 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
            A1 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
            sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
            sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];
    
    LMI3 = [Q                       (A2 * Q + B * Y)'    (sqrt(Q1)*Q)'       (sqrt(R)*Y)'
            A2 * Q + B * Y          Q                    zeros(nx, nx)       zeros(nx, nu);
            sqrt(Q1)*Q              zeros(nx, nx)        gamma * eye(nx, nx) zeros(nx, nu);
            sqrt(R)*Y               zeros(nu, nx)        zeros(nu, nx)       gamma * eye(nu, nu)];
    
    % Input constraints  |uk| <= 2

    LMI4 = [X  Y;
            Y' Q];

    ZEROS = 1e-7;
    Constraints = [LMI1 >= ZEROS, LMI2 >= ZEROS, LMI3 >= ZEROS, LMI4 >= ZEROS, Q >= ZEROS];
    for j = 1 : nu
        Constraints = [Constraints X(j, j) <= 4] ;
    end

    Objective = gamma;
    solvesdp(Constraints, Objective);
    Y = double(Y);
    Q = double(Q);
    
    F_Receding = Y / Q;
        
    alpha = (10 - 0.1) * rand + 0.1;
    A_Receding = A(alpha);
    
    x_next = A_Receding * x_next + B * F_Receding * x_next;
    X_Receding = [X_Receding x_next];
    U_Receding = [U_Receding   F_Receding * x_next];
end

figure
plot(0:0.1:10,X_Receding_static(1,:))
hold on
plot(0:0.1:10,X_Receding(1,:))
axis([0 10,-0.2 1]);
xlabel('time(sec)');
ylabel('$\theta$ (rad) ','interpreter','latex');
title('Angular position $\theta$ (rad)','interpreter','latex');

figure
plot(0:0.1:10,U_Receding_static(1,:))
hold on 
plot(0:0.1:10,U_Receding(1,:))
axis([0 10,-2 0.5]);
xlabel('time(sec)');
ylabel('$u$ (volts) ','interpreter','latex');
title('Control signal $u$ (volts)','interpreter','latex');





















