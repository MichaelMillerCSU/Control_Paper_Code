function [Y0, Q0, theta_p, S] = Copy_of_LMI_Get_K(g, xp)
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
    setlmis([]);
    theta_p = lmivar(1, [1, 1]);
    Y0 = lmivar(2, [nu, nx]);
    Q0 = lmivar(2, [nx, nx]);
    
    S = lmivar(1, [nc, 1]);
    
    lmiterm([-1, 1, 1, 0], 1);
    lmiterm([-1, 1, 2, theta_p], 1, xp');
    lmiterm([-1, 2, 2, Q0], 1, 1);
    
    lmiterm([-2, 1, 1, Q0], 1, 1);
    lmiterm([-2, 2, 2, Q0], 1, 1);
    lmiterm([-2, 3, 3, 0], g.*eye(nx));
    lmiterm([-2, 4, 4, 0], g.*eye(nu));
    lmiterm([-2, 2, 1, Q0], A1, 1);
    lmiterm([-2, 2, 1, Y0], B1, 1);
    lmiterm([-2, 3, 1, Q0], sqrt(Q), 1);
    lmiterm([-2, 4, 1, Y0], sqrt(R), 1);
    
    lmiterm([-3, 1, 1, Q0], 1, 1);
    lmiterm([-3, 2, 2, Q0], 1, 1);
    lmiterm([-3, 3, 3, 0], g.*eye(nx));
    lmiterm([-3, 4, 4, 0], g.*eye(nu));
    lmiterm([-3, 2, 1, Q0], A2, 1);
    lmiterm([-3, 2, 1, Y0], B2, 1);
    lmiterm([-3, 3, 1, Q0], sqrt(Q), 1);
    lmiterm([-3, 4, 1, Y0], sqrt(R), 1);
    
    lmiterm([-4, 1, 1, S], 1, 1);
    lmiterm([-4, 1, 2, Q0], F, 1);
    lmiterm([-4, 1, 3, Y0], G, 1);
    lmiterm([-4, 2, 2, Q0], 1, 1);
    lmiterm([-4, 3, 3, Q0], 1, 1);
    
    for i = 1 : nc 
        e = zeros(1, nc);
        e(i) = 1;
        lmiterm([4 + i, 1, 1, S], e, e');
        lmiterm([-4 - i, 1, 1, 0], 1);
    end
    
    
    
    lmisys=getlmis;
    n = decnbr(lmisys);
    c=zeros(n,1);
    for j = 1 : n
        [theta_pj]=defcx(lmisys,j,theta_p)
        %  [uj]=defcx(lmisys,j,u);
        c(j) = -theta_pj;
    end
    [copt,xopt]=mincx(lmisys,c);
    
    theta_p = dec2mat(lmisys,xopt,theta_p)
    Y0 = dec2mat(lmisys,xopt,Y0)
    Q0 = dec2mat(lmisys,xopt,Q0)
    S = dec2mat(lmisys,xopt,S)
    K = Y0 * inv(Q0)
end