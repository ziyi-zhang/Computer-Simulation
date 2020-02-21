% mass spring
M = [1.5, 2, 1, 2]; % Mass 1*N
L = [-4, 4, 0, 0;
     -1, -1, 4, -3;
     0, 0, 0, 0]; % location 3*N
V = [0, 0, 0, 0;
     0, 0, 0, 0;
     0, 0, 0, 0]; % velocity 3*N
N = length(M);
S = ones(N, N); % spring constant
R0 = calcR(L) ./ 4;

% Constants
G = 9.8;
% simulation parameters
tmax = 200 * 24 * 60 * 60;
clockmax = 1e9;
dt = tmax / clockmax;
% visualization parameters
c = [0, 0, 1;
     0.4, 0.8, 0.5;
     0.8, 0, 0.5;
     0.9, 0.4, 0];
dequeSize = 200;
XHist = zeros(N, dequeSize);
YHist = zeros(N, dequeSize);
ZHist = zeros(N, dequeSize);
dequePtr = 1;
% graphics output
set(gcf, 'double', 'on')
hold on
h = scatter3(L(1, :), L(2, :), L(3, :), 25.*ones(1, N), c, 'filled');
for i = 1:N
    htrail(i) = plot3(XHist(i, :), YHist(i, :), ZHist(i, :), 'Color', c(i, :));
end
hold off

% simulation starts here
for clock = 1:clockmax

    t = clock * dt;
    R = calcR(L);
    T = calcT(R, R0, S);
    for i = 1:N
        a = zeros(3, 1);
        for j = 1:N

            if (i == j), continue;end
            a = a + T(i, j) .* (L(:, j) - L(:, i)) ./ R(i, j) ./ M(i);
        end
        % update velocity
        V(:, i) = V(:, i) + a .* dt;
    end
    for i = 1:N  % updating the location of i-th object

        L(:, i) = L(:, i) + V(:, i) .* dt;
    end
    % re-plot
    set(h, 'xdata', L(1, :), 'ydata', L(2, :), 'zdata', L(3, :));
   
    XHist(:, dequePtr) = L(1, :)';
    YHist(:, dequePtr) = L(2, :)';
    ZHist(:, dequePtr) = L(3, :)';
    dequePtr = dequePtr + 1;
    if (dequePtr > dequeSize), dequePtr = 1;end
    for i = 1:N
        set(htrail(i), 'xdata', [XHist(i, dequePtr:dequeSize), XHist(i, 1:dequePtr-1)],...
                       'ydata', [YHist(i, dequePtr:dequeSize), YHist(i, 1:dequePtr-1)],...
                       'zdata', [ZHist(i, dequePtr:dequeSize), ZHist(i, 1:dequePtr-1)]);
    end
    drawnow
    axis([-5, 5, -5, 5])
    pause(0.01);
end


function [R] = calcR(L)

    N = size(L, 2);
    R = zeros(N, N);
    for i = 1:N
        for j = 1:N
        
            R(i, j) = norm(L(:, i) - L(:, j));
        end
    end
end


function [T] = calcT(R, R0, S)

    T = S .* (R - R0);
end
