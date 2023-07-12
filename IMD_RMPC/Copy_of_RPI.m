clc
clear
close all

% Uncertainty system
A1 = [1 0.1; 0 1];
B1 = [0; 1];
A2 = [1 0.2; 0 1];
B2 = [0; 1.5];
K = [-0.5 -0.3];
A_cl1 = A1 + B1*K;
A_cl2 = A2 + B2*K;
A_CL = (A_cl1 + A_cl2) / 2.0; 
% Constraints [-10 -10] <= x <= [10 10] -1 <= u <= 1

Ac = [eye(2); -eye(2); K ; -K];
bc = [10; 10; 10; 10; 1; 0.4];

% Without uncertainty
L = 2;
As = [A1 A2];
Bs = [B1 B2];
tic
[AS, bS] = RPI_Calculate(As, Bs, Ac, bc, K);
toc
X = [];
for x1 = -6 : 0.02 : 6
    for x2 = -6 : 0.02 : 6
        if AS * [x1; x2] <= bS
            X = [X [x1; x2]];
        end
    end
end

scatter(X(1, :), X(2, :))



