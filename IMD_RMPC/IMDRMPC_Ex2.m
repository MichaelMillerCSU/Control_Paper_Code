clc
clear
close all
%% Not FIXED ! It is the same as Ex1.
global A0;
global A1;
global A2;
global B0;
global B1;
global B2;
global Q;
global R;
global L;
global NL;
global N;
global nx;
global nu;
global nc;
global F;
global G;
global Xk;
global Xk_m;
global Xk_p;
global Uk_m;
global P;
global K;
global Tnx;
global Tnu;
global Enx;
global Enu;
global AS;
global bS;
global PSI1;
global PSI2;

%% Numerical Example
% Offline
A0 = [0.8623 0 0.1472 0;
      0 0.8359 0 0.1052;
      0 0 0.6434 0;
      0 0 0 0.8130];

B0 = [0.0647 0;
      0 0.0647; 
      0 0.0680; 
      0.0619 0];
 
delta_A1 = [0 0 -0.1177 0;
            0 0 0 0;
            0 0 0.3163 0;
            0 0 0 0];
 
delta_A2 = [0 0 0.1177 0;
            0 0 0 0;
            0 0 -0.3163 0;
            0 0 0 0];

delta_B1 = [0 0;
            0 0;
            0 -0.0204;
            0 0];
delta_B2 = [0 0;
            0 0;
            0 0.0204;
            0 0];
A1 = A0 + delta_A1;
A2 = A0 + delta_A2;
B1 = B0 + delta_B1;
B2 = B0 + delta_B2;

L = 2;
NL = 4;
N = 2;
nx = 4;
nu = 2;
Q = diag([20, 0.01, 200, 0]);

R = 0.001*eye(2);

% gamma = 10000;
gamma = 14600;
% xp = [135; 35; 10; 8];
xp = [111; 24; 0.01; 5];

G1 = [
    -1/42.7506 0;
%      1/42.7506 0;
%      -1/17.2494 0;
     1/17.2494 0;
     0 -1/30.5563;
%      0 1/30.5563;
%      0 -1/29.4437;
     0 1/29.4437
     ];
 
F1 = [
    -1/6.4033 0 0 0;
%      1/6.4033 0 0 0;
%      -1/14.0967 0 0 0;
     1/14.0967 0 0 0;
     0 -1/3.0544 0 0;
%      0 1/3.0544 0 0;
%      0 -1/17.4456 0 0;
     0 1/17.4456 0 0;
     0 0 -1/0.1144 0;
%      0 0 1/0.1144 0;
%      0 0 -1/24.3856 0;
     0 0 1/24.3856 0;
     0 0 0 -1/2.5748;
%      0 0 0 1/2.5748;
%      0 0 0 -1/20.9252;
     0 0 0 1/20.9252
     ];
% F = [zeros(size(G1, 1), size(F1, 2)); F1];
% G = [G1;zeros(size(F1, 1), size(G1, 2))];
F = [F1; zeros(size(G1, 1), size(F1, 2))];
G = [zeros(size(F1, 1), size(G1, 2)); G1];
nc = size(F, 1);
% [Y0, Q0, theta_p, S] = LMI_Get_K(gamma, xp);
[Y0, Q0, theta_p, S] = Copy_of_LMI_Get_K(gamma, xp)
% [K, S, E] = dlqr(A0, B0, Q, R);
% K = -K;

K = Y0 * inv(Q0) 
% K = [-2.5141 -0.0808 -2.2485 0.1793;
%      0.2034 0.0093 -5.8007 0.3582];
% K = [-1.6614 -0.0879 -0.2392 -0.3194;
%      -0.0518 0.0537 -4.2735 0.0693];

Phi1 = A1 + B1 * K;
Phi2 = A2 + B2 * K;

Tnx = [zeros(nx, 1 * nx), eye(nx, nx);
       zeros(nx, 2 * nx)];

Tnu = [zeros(nu, 1 * nu), eye(nu, nu);
       zeros(nu, 2 * nu)];

Enx = [eye(nx) zeros(nx, (N - 1) * nx)];

Enu = [eye(nu) zeros(nu, (N - 1) * nu)];

PSI1 = [Tnx zeros(N * nx, N * nu + nx + N * nu + N * nx);
        zeros(N * nu, N * nx) Tnu zeros(N * nu, nx + N * nu + N * nx);
        zeros(nx, N * nx + N * nu), Phi1, B1, zeros(nx, N * nu - nu), Enx, zeros(nx, N * nx - N * nx);
        zeros(N * nu, N * nx + N * nu + nx), Tnu, zeros(N * nu, N * nx);
        zeros(N * nx, N * nx + N * nu + nx + N * nu) Tnx;];

PSI2 = [Tnx zeros(N * nx, N * nu + nx + N * nu + N * nx);
        zeros(N * nu, N * nx) Tnu zeros(N * nu, nx + N * nu + N * nx);
        zeros(nx, N * nx + N * nu), Phi2, B2, zeros(nx, N * nu - nu), Enx, zeros(nx, N * nx - N * nx);
        zeros(N * nu, N * nx + N * nu + nx), Tnu, zeros(N * nu, N * nx);
        zeros(N * nx, N * nx + N * nu + nx + N * nu) Tnx;];

Q_bar_assist = [Enx zeros(nx, N * nu), eye(nx), zeros(nx, N * nu) zeros(nx, N * nx);
                zeros(nu, N * nx) Enu K Enu zeros(nu, N * nx)];
Q_bar = Q_bar_assist' * [Q zeros(nx, nu); zeros(nu, nx) R] * Q_bar_assist;

% P = dlyap(PSI1', PSI1, Q_bar);

% Compute P satisfy Lemma 3
P = LMI_Get_P(PSI1, PSI2, Q_bar);

% Calculating the polyhedron RPI set Omega

Ac = [eye(nx); -eye(nx); K ; -K];
Bc = [14.0967; 17.4456; 24.3856; 20.9252; -(-6.4033); -(-3.0544); -(-0.1144); -(-2.5748); 17.2494; 29.4437; -(-42.7506); -(-30.5563)];
[AS, bS] = RPI_Calculate([A1 A2], [B1 B2], Ac, Bc, K) % Then we got the Omega / Set AS x <= BS

% X = [];
% for x1 = -5 : 1 : 20
%     for x2 = -5 : 1 : 20
%         for x3 = -5 : 1 : 20
%             for x4 = -5 : 1 : 20
%                 if AS * [x1; x2; x3; x4] <= bS
%                     X = [X [x1; x2; x3; x4]];
%                 end
%             end
%         end
%     end
% end

% X = [];
% for x2 = -5 : 0.03 : 30
%     for x3 = -5 : 0.03 : 30
%         if AS * [1; x2; x3; 2] <= bS
%             X = [X [1; x2; x3; 2]];
%         end
%     end
% end

% for x1 = -7 : 0.03 : 20
%     for x4 = -7 : 0.03 : 20
%         if AS * [x1; 7; 5; x4] <= bS
%             X = [X [x1; 7; 5; x4]];
%         end
%     end
% end


% scatter(X(2, :), X(3, :))
% axis([-7 20 -2 21])


%% Initialize the data set
Xk = [];
Uk_m = [];
zero = 1e-208
x0 = [zero; zero; zero; zero];
xk = zero;

for i = 1 : NL + 1
    uk = zero.*ones(nu, 1);
    p1 = rand;
    A = p1 .* A1 + (1 - p1) .* A2;
    B = p1 .* B1 + (1 - p1) .* B2;
    xk = A * xk + B * uk;
    Xk = [xk Xk];
    Uk_m = [uk Uk_m];
end
Xk_m = Xk(:, 2: NL + 1);
Xk_p = Xk(:, 1 : NL);
Uk_m = Uk_m (:, 1 : NL);
%% Get the feasible region (not fixed)
% XX = [];
% for x2 = -5 : 1 : 30
%     for x3 = -5 : 1 : 30
%         [uk, cost, EXITFLAG, check] = solve([[1; x2; x3; 2]]);
%         if size(check, 1) ~= 0 && EXITFLAG == 1
%             XX = [XX [1; x2; x3; 2]];
%         end
%     end
% end
% figure
% scatter(XX(2, :), XX(3, :))
% axis([-7 20 -2 21])




%% MPC Solve
maxStep = 50;
% p1 = @(i)1 ./ exp(3.5 - 2.5.*i);
% p1 = @(i) rand;
% p1 = @(i) 0.5 ;
xk = [2; 1; 1; 2];
% xk = [10;7;6;12];
% xk = [5;6;4;5];
Xk_p = add(xk, Xk_p);
X_log = [xk];
Cost_log = [];
% for i = 0 : maxStep - 1
for i = 0 : maxStep - 1
    p1 = rand;
    A = p1 .* A1 + (1 - p1) .* A2;
    B = p1 .* B1 + (1 - p1) .* B2;
    if i ~= 0
        xk = xk_next;
    end
    [uk, cost, EXITFLAG, check] = solve(xk);
    if size(check, 1) == 0  || EXITFLAG ~= 1
        error("No solution!");
    end
    xk_next = A*xk + B*uk
    Xk_p = add(xk_next, Xk_p);
    Xk_m = add(xk, Xk_m);
    Uk_m = add(uk, Uk_m);
    X_log = [X_log xk_next];
    Cost_log = [Cost_log cost];
end
figure 
plot(X_log(1, :));
hold on 
plot(X_log(2, :));
hold on 
plot(X_log(3, :));
hold on 
plot(X_log(4, :));
hold on 
grid on 
legend('x1', 'x2', 'x3', 'x4')

figure 
plot(Cost_log)


%% Add function
function Set = add(vector, set)
    Set = [vector set];
    Set(:, end) = [];
end






























