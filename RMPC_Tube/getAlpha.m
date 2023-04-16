function alpha = getAlpha()
    global A;
    global B;
    global s;
    global W_Set;
    global epsilon;
    global K;
    alpha = -inf;
    n = size(W_Set.A, 1);
    As = (A + B * K)^s;

    for i = 1 : n
        fi = W_Set.A(i, :)';
        gi = W_Set.b(i);
        alpha = max([alpha, W_Set.support(As * fi) / gi]);
    end
