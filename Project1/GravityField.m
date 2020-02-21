% gravity field test
% Constants
G = 6.674e-11;
% N-body initialization
M = [1e10, 1e6]; % Mass
L = [0, 1e3;
     0, 0]; % Location 2*N
V = [0, 0;
     0, sqrt(G*M(1)/L(1, 2))]; % Velocity 2*N
N = length(M);
% simulation parameters
tmax = 365.25 * 24 * 60 * 60;
clockmax = 1e5;
dt = tmax / clockmax;
% visualization parameters
c = [0, 0, 1;
     0.4, 0.8, 0.5;
     0.8, 0, 0.5;
     0.9, 0.4, 0];
dequeSize = 50;
XHist = ones(N, dequeSize) .* L(1, :)';
YHist = ones(N, dequeSize) .* L(2, :)';
dequePtr = 1; % deque to show trace
%{
axisXMin = min(L(1, :)) - 0.2*(max(L(1, :))-min(L(1, :)));
axisXMax = max(L(1, :)) + 0.2*(max(L(1, :))-min(L(1, :)));
axisYMin = min(L(2, :)) - 0.2*(max(L(2, :))-min(L(2, :)));
axisYMax = max(L(2, :)) + 0.2*(max(L(2, :))-min(L(2, :)));
%}
axisXMin = -1.2 * max(L, [], 'all');
axisXMax =  1.2 * max(L, [], 'all');
axisYMin = -1.2 * max(L, [], 'all');
axisYMax =  1.2 * max(L, [], 'all');
fieldResolution = 40;
axisXSpace = (axisXMax - axisXMin) / (fieldResolution-1);
axisYSpace = (axisYMax - axisYMin) / (fieldResolution-1);
[surfX, surfY] = meshgrid(axisXMin:axisXSpace:axisXMax, axisYMin:axisYSpace:axisYMax);
% graphics output
set(gcf, 'double', 'on')
subplot(1, 2, 1)
hold on
h = scatter(L(1, :), L(2, :), 25.*ones(1, N), c(1:N, :), 'filled');
for i = 1:N
    htrail(i) = plot(XHist(i, :), YHist(i, :), 'Color', c(i, :));
end
hold off
subplot(1, 2, 2)
surf(surfX, surfY, zeros(length(surfX), length(surfY)));


% simulation starts here
for clock = 1:clockmax

    t = clock * dt;
    R = calcR(L);
    for i = 1:N % now updating the velocity of i-th object
        
        a = zeros(2, 1);
        for j = 1:N
        
            if (i==j), continue;end
            a = a + G .* M(j) .* (L(:, j) - L(:, i)) ./ R(i, j)^3;
        end
        % update velocity
        V(:, i) = V(:, i) + a .* dt;
    end
    for i = 1:N  % updating the location of i-th object

        L(:, i) = L(:, i) + V(:, i) .* dt;
    end
    % re-plot subimage-1
    set(h, 'xdata', L(1, :), 'ydata', L(2, :));
   
    XHist(:, dequePtr) = L(1, :)';
    YHist(:, dequePtr) = L(2, :)';
    dequePtr = dequePtr + 1;
    if (dequePtr > dequeSize), dequePtr = 1;end
    for i = 1:N
        set(htrail(i), 'xdata', [XHist(i, dequePtr:dequeSize), XHist(i, 1:dequePtr-1)],...
                       'ydata', [YHist(i, dequePtr:dequeSize), YHist(i, 1:dequePtr-1)]);
    end
    % replot subimage-2
    gradient = zeros(fieldResolution, fieldResolution, 2);
    XArr = axisXMin:axisXSpace:axisXMax;
    YArr = axisYMin:axisYSpace:axisYMax;
    for i = 1:length(XArr)
        for j = 1:length(YArr)
            for k = 1:N % for point (i, j) calc k-th planet's force
                
                p = [XArr(i); YArr(i)];
                r = norm(L(:, k) - p);
                gradient(i, j, 1) = gradient(i, j, 1) + G * M(k) / r^3 .* (L(1, k) - p(1));
                gradient(i, j, 2) = gradient(i, j, 2) + G * M(k) / r^3 .* (L(2, k) - p(2));
            end
        end
    end
    % draw now
    drawnow
    subplot(1, 2, 1)
    xlim([axisXMin, axisXMax]);
    ylim([axisYMin, axisYMax]);
    % pause(0.01);
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
