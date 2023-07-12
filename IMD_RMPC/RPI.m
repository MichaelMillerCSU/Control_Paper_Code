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

AS = [-eye(2); eye(2); -K; K];
bS = [10; 10; 10; 10; 0.4; 1];

% Without uncertainty
% L = 2;

% Testing for the set
% X = [];
% for k = 1 : 50
%     x = 10 * rand(2, 1) - 5;
%     for i = 1 : 1000
%         X = [X x];
%         if rand > 0.5
%             x = A_cl1 * x;
%         else
%             x = A_cl2 * x;
%         end
%     end
% end
% scatter(X(1, :), X(2, :))
% hold on 
% i = 1;
% while 1
%     if i > size(AS, 1)
%         break;
%     end
%     a = AS(i, :);
%     b = bS(i, :);
%     for j = 1 : L
%         if j == 1
%             Phi_i = A_cl1;
%         else
%             Phi_i = A_cl2;
%         end
%         f = a*Phi_i;
%         x = linprog(-f, AS, bS);
%         c_i = f * x - b;
%         if c_i > 0
%             AS = [AS; f];
%             bS = [bS; b];
%         end
%     end
%     
%     i = i + 1;
%     
% end
[AS, bS] = RPI_Calculate([A1 A2], [B1 B2], AS, bS, K) % Then we got the Omega / Set AS x <= BS

X = [];
for x1 = -6 : 0.1 : 6
    for x2 = -6 : 0.1 : 6
        if AS * [x1; x2] <= bS
            X = [X [x1; x2]];
        end
    end
end

scatter(X(1, :), X(2, :))


