% runTrafficSimulationCongestion

o.roadArray = struct(...
    ...% road idx  1   2   3   4   5   6   7   8   9   10   11   12   13
    'nodeStart',  {1,  2,  3,  4,  5,  6,  7,  8,  2,  2,   1,   4,   4},...
    'nodeEnd',    {2,  1,  1,  1,  1,  2,  2,  2,  8,  7,   4,   5,   3},...
    'length',     {200,200,283,200,283,283,200,283,283,200, 200, 200, 200},...
    'imageRoad',  {2,  1,  0,  11, 0,  0,  10, 9,  8,  7,   4,   0,   0}...
);

o.nodeArray = struct(...
    ...% node idx  1    2    3    4    5    6    7    8
    'spawnRate',  {0.1, 0.1, 0.4, 0.4, 0.4, 0.5, 0.5, 0.5},...
    'destChance', {1,   1,   2,   2,   2,   1,   2,   2}...
);

% visArray can be omitted when simulating, but it is necessary for
% visualization
o.visArray = struct(...
    ...% node idx  1    2    3    4    5    6    7    8
    'x',          {0.0, 0.0, -0.6,0.0, 0.6, -0.6,0,   0.6},...
    'y',          {-0.3,0.3, -0.9,-0.9,-0.9,0.9, 0.9, 0.9}...
);

% car speed
o.vmax = 22.2;  % meters/second
o.dmin = 5;  % meters
o.dmax = 100;  % meters
% Simulation speed
o.simulationTime = 200;  % in seconds
o.dt = 0.1;  % in seconds
o.fastForward = 2;

% Run Simulation
TrafficSimulation(o);
