G = 9.8;
Ms = 1.989 * 10^30;
X = 152e9;
Y = 0;
U = 0;
V = 0.5*sqrt(G*Ms/X); % circular
tmax = 3 * 365.25 * 24 * 60 * 60;
clockmax = 1e9;
dt = tmax / clockmax;

XHist = zeros(clockmax, 1);
YHist = zeros(clockmax, 1);
tHist = zeros(clockmax, 1);

set(gcf, 'double', 'on')
plot(0, 0, 'r*')
hold on
h = plot(X, Y, 'bo');
htrail = plot(X, Y);
hold off

for clock = 1:clockmax

    t = clock * dt;
    R = sqrt(X^2 + Y^2);
    U = U - dt * G * Ms * X / R^3;
    V = V - dt * G * Ms * Y / R^3;
    X = X + dt * U;
    Y = Y + dt * V;

    set(h, 'xdata', X, 'ydata', Y);
    set(htrail, 'xdata', XHist(1:clock-1), 'ydata', YHist(1:clock-1));
    drawnow
    axis([-1.2*X, 1.2*X, -1.2*X, 1.2*X])
    axis manual
    pause(0.01)
    XHist(clock) = X;
    YHist(clock) = Y;
    tHist(clock) = t;
end

