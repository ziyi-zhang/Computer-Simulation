% RunSpringTribody

% N-body initialization
M = [1, 1, 1]; % Mass
L = [10, -5, -5;
     0, 8.66, -8.66]; % Location 2*N
V = [0, 0, 0;
     0, 0, 0]; % Velocity 2*N
A = [0, 0, 0]; % Whether anchor point 0/1
N = length(M);
% Spring connection and property
% suggested R0=10
Spring = [
    % IdxA idxB R0   S
      1    2    10   1;
      2    3    10   1;
      3    1    10   1;
    ];
% simulation parameters
tmax = 10;
clockmax = 2000;
dt = tmax / clockmax;
% visualization parameters
visRealtime = true;
storeData = true;
% visRealtime parameters
pauseTime = 0.0;
dequeSize = 200;

% call SpringSimulation
SpringSimulation
