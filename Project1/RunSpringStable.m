% RunSpringStable

% N-body initialization
M = [10 10 5 20 10  0.1 0.1 0.1 75/2]; % Mass
s3 = sqrt(3);
L = [0, 0,  7.5,     0,  -7.5,   0, 2.5*s3,-2.5*s3,  1;
     0, 10, 7.5*s3,  -5, 0,      5, -2.5,  -2.5,    -s3]; % Location 2*N
L(:, 6:8) = L(:, 6:8) + repmat(L(:, 2), 1, 3);
V = [0, 10, 30*s3,  -5, 0,      50, -25,   -25,   -4*s3;
     0, 0, -30,     0, 30,      0, -25*s3,25*s3,  -4]; % Velocity 2*N
V(:, 6:8) = V(:, 6:8) + repmat(V(:, 2), 1, 3);
A = [0, 0, 0, 0, 0, 0, 0, 0, 0]; % Whether anchor point 0/1
N = length(M);
% Spring connection and property
% suggested R0=10
Spring = [
    % IdxA idxB R0   S
      1    2    5   20;
      1    3    5   120;
      1    4    2.5 40;
      1    5    5   480;
      2    6    4   50;
      2    7    4   50;
      2    8    4   50;
      1    9    1   1200;
    ];
% simulation parameters
tmax = 2.5;
dt = 0.0002;
clockmax = tmax / dt;
% visualization parameters
visRealtime = false;
storeData = true;
% visRealtime parameters
pauseTime = 0.0;
dequeSize = 600;

% call SpringSimulation
SpringSimulation
