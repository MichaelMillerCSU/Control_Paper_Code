function [Y0, Q0, theta_p] = LMI_Get_K(g, xp)
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
    global ncu;
    global ncx;
    global F;
    global G;
    global G1;
    global F1;
    theta_p = sdpvar(1, 1);
    Q0 = sdpvar(nx, nx);
    Y0 = sdpvar(nu, nx);
%     Y1 = sdpvar(nu, nx);
%     Y2 = sdpvar(nu, nx);
    S = sdpvar(nc, nc);
%     S1 = sdpvar(ncu, ncu);
    S2 = sdpvar(ncx, ncx);
    S3 = sdpvar(ncu, ncu);
    N1 = [1 theta_p.*xp';
          theta_p.*xp Q0];
    
    M1 = [Q0 (A1*Q0+B1*Y0)' (sqrt(Q)*Q0)' (sqrt(R)*Y0)';
          A1*Q0+B1*Y0 Q0 zeros(nx, nx) zeros(nx, nu);
          sqrt(Q)*Q0 zeros(nx, nx) g*eye(nx) zeros(nx, nu);
          sqrt(R)*Y0 zeros(nu, nx) zeros(nu, nx) g*eye(nu)];
      
    M2 = [Q0 (A2*Q0+B2*Y0)' (sqrt(Q)*Q0)' (sqrt(R)*Y0)';
          A2*Q0+B2*Y0 Q0 zeros(nx, nx) zeros(nx, nu);
          sqrt(Q)*Q0 zeros(nx, nx) g*eye(nx) zeros(nx, nu);
          sqrt(R)*Y0 zeros(nu, nx) zeros(nu, nx) g*eye(nu)];

%     N2 = [S G*Y0 F*Q0;
%           (G*Y0)' Q0 zeros(nx, nx);
%           (F*Q0)' zeros(nx, nx) Q0];
%     N2 = [S F*Q0 + G*Y0;
%           (F*Q0 + G*Y0)' Q0];
    N2 = [S2 F1*Q0;
          (F1*Q0)' Q0];

    N3 = [S3 G1*Y0;
          (G1*Y0)' Q0];

%     N3 = [S2 F1*(A1 * Q0 + B1 * Y0);
%          (F1*(A1 * Q0 + B1 * Y0))' Q0];
% 
%     N4 = [S3 F1*(A2 * Q0 + B2 * Y0);
%          (F1*(A2 * Q0 + B2 * Y0))' Q0];
    Constraints = [M1 >= 0; M2 >= 0; N1 >= 0; N2 >= 0; N3 >= 0];
%     Constraints = [M1 >= 0; M2 >= 0; N1 >= 0; N2 >= 0];
%     for i = 1 : ncu
%         Constraints = [Constraints; S1(i, i) <= 1];
%     end
%     for i = 1 : nc
%         Constraints = [Constraints; S(i, i) <= 1];
%     end
    for i = 1 : ncx
        Constraints = [Constraints; S2(i, i) <= 1];
    end
    for i = 1 : ncu
        Constraints = [Constraints; S3(i, i) <= 1];
    end
    obj = -theta_p;
    options = sdpsettings('solver','mosek')
    solvesdp(Constraints, obj,options)
      
    theta_p = double(theta_p)
    Q0 = double(Q0)
    Y0 = double(Y0)
%     S = double(S)
    
end