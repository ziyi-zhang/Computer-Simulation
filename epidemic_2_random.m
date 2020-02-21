% epidemic-2-random
a = 0.01;
b = 0.0004;
S = 1e5;
I = 10;
R = 0;
N = S + I + R;
% simulation parameters
clockmax = 2e4;
dt = 1;
% visualization parameters
SHist = zeros(1, clockmax);
IHist = zeros(1, clockmax);
RHist = zeros(1, clockmax);

% simulation starts here
for clock = 1:clockmax

    t = clock * dt;
    % SI = dt * a * (I/N) * S;
    SI = 0;
    for i = 1:S
        if (rand < a*(I/N)*dt)
            SI = SI + 1;
        end
    end
    % IR = dt * b * I;
    IR = 0;
    for i = 1:I
        if (rand < b*dt)
            IR = IR + 1;
        end
    end
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
plot(tArr, IHist, tArr, SHist, tArr, RHist);
hold off
subplot(2, 1, 2);
plot(tArr, IHist);
