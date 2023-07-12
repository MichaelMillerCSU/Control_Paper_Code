function P = LMI_Get_P(PSI1, PSI2, Q_BAR)
    setlmis([]);
    X = lmivar(1, [28, 1]);
    lmiterm([1 1 1 X], PSI1', PSI1);
    lmiterm([1 1 1 X], -1, 1);
    lmiterm([1 1 1 0], Q_BAR);
    
    lmiterm([2 1 1 X], PSI2', PSI2);
    lmiterm([2 1 1 X], -1, 1);
    lmiterm([2 1 1 0], Q_BAR);
    
    lmiterm([-3 1 1 X], 1, 1);
    
    lmis = getlmis;
    [tmin, xfeas] = feasp(lmis);
    if tmin < 0
        P = dec2mat(lmis,xfeas,X);
    else
        error("LMI NO SOLUTION!")
    end
end