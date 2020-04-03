% Main simulation file
function [] = TrafficSimulation(o)

%% Phase input variables
checkInput();
vmax = o.vmax;
dmin = o.dmin;
dmax = o.dmax;
simulationTime = o.simulationTime;
dt = o.dt;
clockMax = ceil(simulationTime / dt);
roadArray = o.roadArray;
nodeArray = o.nodeArray;
%% Initialization
% carArray
carArray.velocity = 0;
carArray.roadNum = 0;
carArray.position = 0;
carArray.frontCar = 0;
carArray.backCar = 0;
carArray.alive = 0;
carArray = repmat(carArray, 1, 10);  % initialize a car pool with ten inactivated cars
%% Simulation
for clock = 1:clockMax
    t = clock * dt;
    
end


%% nested function: do some basic check of whether input variables are valid
function [] = checkInput()
    
    % check imageRoad
    for i = 1:length(o.roadArray)
        if (o.roadArray(i).imageRoad > 0)
            j = o.roadArray(i).imageRoad;
            if (o.roadArray(i).nodeStart ~= o.roadArray(j).nodeEnd ||...
                o.roadArray(i).nodeEnd ~= o.roadArray(j).nodeStart ||...
                o.roadArray(i).length ~= o.roadArray(j).length)
                error("roadArray input invalid: 'imageRoad'.");
            end
        end
    end
    % check nodeArray
    for i = 1:length(o.nodeArray)
        if (o.nodeArray(i).spawnRate < 0 || o.nodeArray(i).spawnRate > 1)
            error("nodeArray input invalid: 'spawnRate'.");
        end
        if (o.nodeArray(i).destChance < 0)
            error("nodeArray input invalid: 'destChance'.");
        end
    end
end
end
