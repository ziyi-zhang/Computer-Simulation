% runTrafficSimulationCityEL

o.roadArray = struct(...
    ...% road idx  1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26
    'nodeStart',  {1,  1,  1,  1,  2,  2,  3,  3,  3,  4,   4,   5,   5,   6,   6,   6,   6,   7,   7,   8,   8,   9,   10,  11,  11,  12,...
                   2,  3,  4,  5,  6,  7,  4,  8,  9,  5,   10,  11,  18,  11,  12,  13,  14,  8,   14,  9,   15,  16,  17,  12,  18,  13,...
                   20, 15, 16, 17, 18, 19, 13, 14,...
                   15, 16, 17, 18, 19, 13, 14, 20},...
    'nodeEnd',    {2,  3,  4,  5,  6,  7,  4,  8,  9,  5,   10,  11,  18,  11,  12,  13,  14,  8,   14,  9,   15,  16,  17,  12,  18,  13,...
                   1,  1,  1,  1,  2,  2,  3,  3,  3,  4,   4,   5,   5,   6,   6,   6,   6,   7,   7,   8,   8,   9,   10,  11,  11,  12,...
                   15, 16, 17, 18, 19, 13, 14, 20,...
                   20, 15, 16, 17, 18, 19, 13, 14},...
    'length',     {120,133,167,167,137,243,137,211,167,200, 170, 203, 224, 316, 213, 236, 233, 240, 203, 307, 236, 180, 120, 180, 194, 120,...
                   120,133,167,167,137,243,137,211,167,200, 170, 203, 224, 316, 213, 236, 233, 240, 203, 307, 236, 180, 120, 180, 194, 120,...
                   120,317,97, 184,164,154,177,195,...
                   120,317,97, 184,164,154,177,195},...
    'imageRoad',  {27, 28, 29, 30, 31, 32, 33, 34, 35, 36,  37,  38,  39,  40,  41,  42,  43,  44,  45,  46,  47,  48,  49,  50,  51,  52,...
                   1   2   3   4   5   6   7   8   9   10   11   12   13   14   15   16   17   18   19   20   21   22   23   24   25   26,...
                   61, 62, 63, 64, 65, 66, 67, 68,...
                   53, 54, 55, 56, 57, 58, 59, 60}...
);

o.nodeArray = struct(...
    ...% node idx  1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17   18    :  19    20
    'spawnRate',  {1.6, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,1.2, 0.1, 0.1, 0.1,1.2, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1,      0,    0},...
    'destChance', {8,   1,   1,   1,   1,   1,   1,   8,   1,   1,   1,   8,   2,   2,   3,   3,   2,   2,      0,    0}...
);

% visArray can be omitted when simulating, but it is necessary for
% visualization
o.visArray = struct(...
    ...% node idx  1    2    3     4    5     6      7      8     9     10    11    12    13    14     15    16     17     18    :  19    20
    'x',          {0.0, -0.2,0,    0.27,0.27, -0.27, -0.67, -0.4, 0.2,  0.6,  0.33, 0,    -0.2, -0.73, -0.73, 0.53, 0.73,  0.67,    0.4,  -0.97},...
    'y',          {0.0, 0.13,-0.27,-0.2,0.2,  0.4,   0,     -0.4, -0.53,-0.13,0.6,  0.73, 0.87, 0.4,   -0.73, -0.67,-0.33, 0.4,     0.97, -0.33}...
);

% car speed
o.vmax = 22.2;  % meters/second
o.dmin = 5;  % meters
o.dmax = 100;  % meters
% Simulation speed
o.simulationTime = 200;  % in seconds
o.dt = 0.1;  % in seconds
o.fastForward = 3;

% Run Simulation
TrafficSimulation(o);
