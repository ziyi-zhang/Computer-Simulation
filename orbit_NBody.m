% N-body initialization
G = 9.8;
Ms = 1.989e30;
Me = 5.972e24;
Mm = 7.348e22;
Dse = 147.5e9;
Dem = 384.4e6;
M = [Ms, Me, Mm]; % Mass
L = [0, Dse, Dse+Dem;
     0, 0, 0;
     0, 0, 0]; % Location 3*N
V = [0, 0, 0;
     0, sqrt(G*Ms/Dse), sqrt(G*Me/Dem)+sqrt(G*Ms/Dse);
     0, 0, 0]; % Velocity 3*N
N = length(M);
% simulation parameters
tmax = 3 * 365.25 * 24 * 60 * 60;
clockmax = 1e8;
dt = tmax / clockmax;
% visualization parameters
c = [0, 0, 1;
     0.4, 1, 0.5;
     1, 0.3, 0];
XHist = zeros(N, clockmax);
YHist = zeros(N, clockmax);
ZHist = zeros(N, clockmax);
% graphics output
set(gcf, 'double', 'on')
hold on
h = scatter3(L(1, :), L(2, :), L(3, :), [30, 15, 7].*2, c, 'filled');
htrail = plot3(XHist, YHist, ZHist);
hold off

% simulation starts here
for clock = 1:clockmax

    t = clock * dt;
    R = calcR(L);
    for i = 1:N % now updating the velocity of i-th object
        
        a = zeros(3, 1);
        for j = 1:N
        
            if (i==j), continue;end
            a = a + G .* M(j) .* (L(:, j) - L(:, i)) ./ R(i, j)^3;
        end
        % update velocity
        V(:, i) = V(:, i) + a .* dt;
    end
    for i = 1:N % now updating the location of i-th object
    
        L(:, i) = L(:, i) + V(:, i) .* dt;
    end
    % save value and re-plot
    set(h, 'xdata', L(1, :), 'ydata', L(2, :), 'zdata', L(3, :));
    XHist(:, i) = L(1, :);
    YHist(:, i) = L(2, :);
    ZHist(:, i) = L(3, :);
    set(htrail, 'xdata', XHist, 'ydata', YHist, 'zdata', ZHist);
    drawnow
    axis([-1.2*Dse, 1.2*Dse, -1.2*Dse, 1.2*Dse])
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
