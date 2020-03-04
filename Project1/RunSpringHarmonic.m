% RunSpringHarmonic

% N-body initialization
M = [1 1]; % Mass
L = [0, 15;
     0, 0]; % Location 2*N
V = [0, 0;
     0, 0]; % Velocity 2*N
A = [1, 0]; % Whether anchor point 0/1
N = length(M);
% Spring connection and property
% suggested R0=10
Spring = [
    % IdxA idxB R0   S
      1    2    10   1;
    ];
% simulation parameters
tmax = 20;
clockmax = 2000;
dt = tmax / clockmax;
% visualization parameters
visRealtime = true;
storeData = true;
% visRealtime parameters
pauseTime = 0.0;
dequeSize = 1;

% call SpringSimulation
SpringSimulation
