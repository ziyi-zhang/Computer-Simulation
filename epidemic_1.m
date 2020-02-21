% epidemic-1
a = 0.01;
b = 0.0002;
S = 1e10;
I = 10;
R = 0;
N = S + I + R;
% simulation parameters
clockmax = 1e5;
dt = 1;
% visualization parameters
SHist = zeros(1, clockmax);
IHist = zeros(1, clockmax);
RHist = zeros(1, clockmax);

% simulation starts here
for clock = 1:clockmax

    t = clock * dt;
    SI = dt * a * (I/N) * S;
    IR = dt * b * I;
    S = S - SI;
    I = I + SI - IR;
    R = R + IR;
    % save data
    SHist(clock) = S;
    IHist(clock) = I;
    RHist(clock) = R;
end

% plot
figure
tArr = dt .* (1:clockmax);
subplot(2, 1, 1);
hold on
plot(tArr, SHist, tArr, IHist, tArr, RHist);
hold off
subplot(2, 1, 2);
plot(tArr, IHist);
