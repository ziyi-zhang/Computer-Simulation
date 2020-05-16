% runTrafficSimulationPath

o.roadArray = struct(...
    ...% road idx  1   2   3  
    'nodeStart',  {1,  1,  3},...
    'nodeEnd',    {2,  3,  2},...
    'length',     {400,283,283},...
    'imageRoad',  {0,  0,  0}...
);

o.nodeArray = struct(...
    ...% node idx  1    2    3  
    'spawnRate',  {2.0, 0,   0},...
    'destChance', {0,   1,   0}...
);

% visArray can be omitted when simulating, but it is necessary for
% visualization
o.visArray = struct(...
    ...% node idx  1    2    3 
    'x',          {0.0, 0.0,-0.5},...
    'y',          {0.5,-0.5, 0}...
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
