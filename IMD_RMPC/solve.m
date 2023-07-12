function [uk, cost, EXITFLAG, kesi, output] = solve(xk)
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
    global Ac;
    global bc;
    global PSI1;
    global PSI2;
    Ncol = N*(NL+nu) - N*(N-1) / 2 + 1;
    Nrow = 2*N*(nx+nu)+nx;
    Dk = [Xk_m(:, 1 : NL) zeros(nx, Ncol - NL);
          zeros(nx, NL) Xk_p(:, 2 : NL) zeros(nx, Ncol - NL - (NL - 1));
          Uk_m(:, 1 : NL) zeros(nu, Ncol - NL);
          zeros(nu, NL) Uk_m(:, 1 : NL - 1) zeros(nu, Ncol - NL - (NL - 1));
          -Xk_m(:, 1 : NL) zeros(nx, Ncol - NL - 1 - nu*N ) xk zeros(nx, nu*N);
          zeros(nu * N, Ncol - nu*N) eye(nu*N);
          Xk_p(:, 1 : NL) -Xk_p(:, 2 : NL) zeros(nx, Ncol - NL - (NL - 1));
          zeros(nx, NL) Xk_p(:, 1 : NL - 1) zeros(nx, Ncol - NL - (NL - 1))];
    Mk = Dk' * P * Dk;
%     [x,fval,exitflag,output,lambda] =quadprog(H,f,A,b,[],[],lb);
%   Variables : kesi = col[theta_sub(k), 1, c_sub(0|k)]; so that
%   theta_sub(k) = col[theta(0|k), theta(1|k), theta(2|k),...,theta(N-1|k)];
%   theta(i|k)\in {R^(NL - i)} NL + NL - 1 + NL - 2 + ... + NL - (N - 1)
%   = N * NL - N  * (N - 1) / 2
%   c_sub(i|k) = col[c(i|k), c(i+1|k),...,c(i+N-1|k)] c(i|k)\in {R^nu}
%   = N * nu
%   So kesi = col[theta_sub(k), 1, c_sub(0|k)]; that is N * NL - N  * (N - 1) / 2 + N * nu + 1
%     H = 2.*Mk;
    H = (Mk + Mk');
    symmetric = norm(H-H',inf);
    f = zeros(1, N*(NL+nu) - N*(N - 1) / 2 + 1);
%   zeta(0) = Dk * kesi(k);
%   zeta(0) = col[z_sub(0|k),v_sub(0|k),e(0|k),c_sub(0|k),W_sub(0|K)];
%   xxx_sub(i|k) = col[xxx(i|k),xxx(i + 1|k),...,xxx(i + N - 1|k)]

%% 初始版本
    Aeq_cons_0 = [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * Dk;
    beq_cons_0 = [xk];
    
    Aineq_cons0 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ])* Dk;
    bineq_cons0 = [ones(nc, 1)];
%% 这个约束写法等价
%     Aineq_consF = F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * Dk;
%     bineq_consF = [ones(nx * 2, 1);zeros(nc - nx * 2, 1)];
%     Aineq_consG = G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]* Dk;
%     bineq_consG = [zeros(nc - nu * 2, 1) ;ones(nu * 2, 1)];
%     Aineq_cons0 = [Aineq_consF;Aineq_consG];
%     bineq_cons0 = [bineq_consF;bineq_consG];
%% 修改了一些写法
%     Aineq_cons1 = (F * [zeros(nx, nx) eye(nx) zeros(nx, N * nu) A1+B1*K B1 zeros(nx, nu) eye(nx) zeros(nx, nx)]...
%                    + G * [zeros(nu, N * nx + nu) eye(nu) K*(A1+B1*K) K*B1 eye(nu) K*eye(nx) zeros(nu, nx)]) *Dk;
%     bineq_cons1 = [ones(nc, 1)];
    Aineq_cons1 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI1 * Dk;
    bineq_cons1 = [ones(nc, 1)];
%     Aineq_cons2 = (F * [zeros(nx, nx) eye(nx) zeros(nx, N * nu) A2+B2*K B2 zeros(nx, nu) eye(nx) zeros(nx, nx)]...
%                    + G * [zeros(nu, N * nx + nu) eye(nu) K*(A2+B2*K) K*B2 eye(nu) K*eye(nx) zeros(nu, nx)]) *Dk;
%     bineq_cons2 = [ones(nc, 1)];
    Aineq_cons2 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI2 * Dk;
    bineq_cons2 = [ones(nc, 1)];
    Aineq_cons3_1 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI1 * PSI2 * Dk;
    bineq_cons3_1 = [ones(nc, 1)];
    Aineq_cons3_2 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI2 * PSI1 * Dk;
    bineq_cons3_2 = [ones(nc, 1)];
    Aineq_cons3_3 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI1 * PSI1 * Dk;
    bineq_cons3_3 = [ones(nc, 1)];
    Aineq_cons3_4 = (F * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]...
                   + G * [zeros(nu, N * nx) eye(nu) zeros(nu, nu) K eye(nu) zeros(nu, nu + N * nx) ]) * PSI2 * PSI2 * Dk;
    bineq_cons3_4 = [ones(nc, 1)];
%     cons3 已经包含在了AS * x <= bS 里面
%     Aineq_cons3_11 = ((F + G*K) * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]) * PSI1 * PSI1 * Dk;
%     bineq_cons3_11 = [ones(nc, 1)];
%     Aineq_cons3_21 = ((F + G*K) * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]) * PSI2 * PSI1 * Dk;
%     bineq_cons3_21 = [ones(nc, 1)];
%     Aineq_cons3_12 = ((F + G*K) * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]) * PSI1 * PSI2 * Dk;
%     bineq_cons3_12 = [ones(nc, 1)];
%     Aineq_cons3_22 = ((F + G*K) * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)]) * PSI2 * PSI2 * Dk;
%     bineq_cons3_22 = [ones(nc, 1)];
    Aineq_terminal1_1 = AS * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * PSI1 * PSI1 * Dk;
%     Aineq_terminal1_1 = AS * [zeros(nx, N * nx) zeros(nx, N * nu) (A1+B1*K)*(A1+B1*K) (A1+B1*K)*B1 B1 (A1+B1*K)*eye(nx) eye(nx, nx)]*Dk;
    bineq_terminal1_1 = [bS];
    Aineq_terminal2_1 = AS * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * PSI2 * PSI1 * Dk;
%     Aineq_terminal2_1 = AS * [zeros(nx, N * nx) zeros(nx, N * nu) (A2+B2*K)*(A1+B1*K) (A2+B2*K)*B1 B2 (A2+B2*K)*eye(nx) eye(nx, nx)]*Dk;
    bineq_terminal2_1 = [bS];
    Aineq_terminal1_2 = AS * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * PSI1 * PSI2 * Dk;
%     Aineq_terminal1_2 = AS * [zeros(nx, N * nx) zeros(nx, N * nu) (A1+B1*K)*(A2+B2*K) (A1+B1*K)*B2 B1 (A1+B1*K)*eye(nx) eye(nx, nx)]*Dk;
    bineq_terminal1_2 = [bS];
    Aineq_terminal2_2 = AS * [eye(nx) zeros(nx, N * nu + nx) eye(nx) zeros(nx, N * nu + N * nx)] * PSI2 * PSI2 * Dk;
%     Aineq_terminal2_2 = AS * [zeros(nx, N * nx) zeros(nx, N * nu) (A2+B2*K)*(A2+B2*K) (A2+B2*K)*B2 B2 (A2+B2*K)*eye(nx) eye(nx, nx)]*Dk;
    bineq_terminal2_2 = [bS];
    A_eq = [Aeq_cons_0];
    b_eq = [beq_cons_0];
    A_ineq = [Aineq_cons0;
              Aineq_cons1;
              Aineq_cons2;
              Aineq_cons3_1;
              Aineq_cons3_2;
              Aineq_cons3_3;
              Aineq_cons3_4;
              Aineq_terminal1_1;
              Aineq_terminal2_1;
              Aineq_terminal1_2;
              Aineq_terminal2_2];
    b_ineq = [bineq_cons0;
              bineq_cons1;
              bineq_cons2;
              bineq_cons3_1;
              bineq_cons3_2;
              bineq_cons3_3;
              bineq_cons3_4;
              bineq_terminal1_1;
              bineq_terminal2_1;
              bineq_terminal1_2;
              bineq_terminal2_2];

    [kesi,feval,EXITFLAG,output] = quadprog(H,f,A_ineq,b_ineq,A_eq,b_eq,[],[],[]);

    cost = feval;
    theta_0 = kesi(1 : NL)
    e_0 = [zeros(nx, N * nu + N * nx) eye(nx) zeros(nx, N * nu + N * nx)] * Dk * kesi
    c_0 = [zeros(nu, N * nu + N * nx + nx) eye(nu) zeros(nu, nu + N * nx)] * Dk * kesi
    uk = (Uk_m * theta_0) + K * e_0 + c_0


end
