% Calculate the shortest path of any two nodes with 'roadArray'
% 'N' is the number of nodes
function [path] = FloydShortestPath(roadArray, N)

    cost = ones(N, N) .* 1e10;
    path = zeros(N, N);
    for i = 1:length(roadArray)
        cost(roadArray(i).nodeStart, roadArray(i).nodeEnd) = roadArray(i).length;
        path(roadArray(i).nodeStart, roadArray(i).nodeEnd) = -1;  % indicates a direct path
    end
    
    for k = 1:N
        for i = 1:N
            for j = 1:N
                if (cost(i, j) > cost(i, k) + cost(k, j))
                    cost(i, j) = cost(i, k) + cost(k, j);
                    path(i, j) = k;
                end
            end
        end
    end
end
