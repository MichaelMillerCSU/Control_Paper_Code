function Ms = getMs()
    global A;
    global B;
    global s;
    global W_Set;
    global epsilon;
    global K;

    hWpos = 0;
    hWneg = 0;

    for i = 1 : s 
        Ai = (A+B*K)^i;
        hWpos = hWpos + W_Set.support(Ai);
        hWneg = hWneg + W_Set.support(-Ai);
    end

    Ms = max([hWpos; hWneg]);










