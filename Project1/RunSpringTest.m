% RunSpringTest

% N-body initialization
M = [1.5, 2, 3]; % Mass
L = [-7, 7, 0;
     5, 4, -6]; % Location 2*N
V = [0, 0, 0;
     0, 0, 0]; % Velocity 2*N
A = [1, 0, 1]; % Whether anchor point 0/1
N = length(M);
% Spring connection and property
% suggested R0=10
Spring = [
    % IdxA idxB R0   S
      1    2    10   1;
%      1    3    10   1;
      2    3    10   1
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

% call GravityField
GravityField
