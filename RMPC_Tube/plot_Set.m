function P_Z = plot_Set(A, color, linewidth)
    Vertex = A.V;
    Vertex = round(Vertex, 5);
    idx = convhull(Vertex(:, 1), Vertex(:, 2));
%     A = Polyhedron([A.V(idx, 1), A.V(idx, 2)]);
%     plot(A.V(:, 1), A.V(:, 2));
    P_Z = plot(A.V(idx, 1), A.V(idx, 2), 'color', color, 'LineWidth', linewidth);
    hold on
end