% runTrafficSimulationTest

o.roadArray = struct(...
    ...% road idx  1  2 
    'nodeStart',  {1, 2},...
    'nodeEnd',    {2, 1},...
    'length',     {1e3, 1e3},...
    'imageRoad',  {2, 1}...
);

o.nodeArray = struct(...
    ...% node idx  1  2
    'spawnRate',  {1, 1},...
    'destChance', {1, 1}...
);

% visArray can be omitted when simulating, but it is necessary for
% visualization
o.visArray = struct(...
    ...% node idx  1     2
    'x',          {-0.5, 0.5},...
    'y',          {0,    0}...
);

% car speed
o.vmax = 22.2;  % meters/second
o.dmin = 5;  % meters
o.dmax = 100;  % meters
% Simulation speed
o.simulationTime = 10;  % in seconds
o.dt = 0.1;  % in seconds

% Run Simulation
TrafficSimulation(o);
