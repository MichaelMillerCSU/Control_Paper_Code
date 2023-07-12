function [AS, bS] = RPI_Calculate(As, Bs, Ac, Bc, K)
%     nx = size(As, 1);
%     L = size(As, 2) / nx;
%     nu = size(Bs, 2) / L;
%     Ac = [eye(nx); -eye(nx); K ; -K];
%     bc = [repmat(ubx, nx, 1); -repmat(lbx, nx, 1); repmat(ubu, nu, 1); -repmat(lbu, nu, 1)];
%     AS = Ac;
%     bS = bc;
%     i = 1;
%     while 1
%         if i > size(AS, 1)
%             break;
%         end
%         a = AS(i, :);
%         b = bS(i, :);
%         for j = 0 : L - 1
%             Phi = As(:, j * nx + 1 : j * nx + nx) + Bs(:, j * nu + 1 : j * nu + nu) * K;
%             f = a*Phi;
%             x = linprog(-f, AS, bS);
%             c_i = f * x - b;
%             if c_i > 0
%                 AS = [AS; f];
%                 bS = [bS; b];
%             end
%         end
%         i = i + 1;
%     end
    nx = size(As, 1);
    L = size(As, 2) / nx;
    nu = size(Bs, 2) / L;
%     Ac = [eye(nx); -eye(nx); K ; -K];
%     bc = [repmat(ubx, nx, 1); -repmat(lbx, nx, 1); repmat(ubu, nu, 1); -repmat(lbu, nu, 1)];
%   Ac x <= Bc
    AS = Ac;
    bS = Bc;
    i = 1;
    while 1
        if i > size(AS, 1)
            break;
        end
        a = AS(i, :);
        b = bS(i, :);
        for j = 0 : L - 1
            Phi = As(:, j * nx + 1 : j * nx + nx) + Bs(:, j * nu + 1 : j * nu + nu) * K;
            f = a * Phi;
            [x, fval] = linprog(-f, AS, bS);
            c_i = -fval - b;
            if c_i > 0
                AS = [AS; f];
                bS = [bS; b];
            end
        end
        i = i + 1;
    end
end