% Main simulation file
% Ziyi. April 2020.
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
carArray.roadArray = [];
carArray = repmat(carArray, 1, 10);  % initialize a car pool with ten inactivated cars
% roadArray
for i = 1:length(roadArray)
    roadArray(i).frontCar = 0;
    roadArray(i).backCar = 0;
    roadArray(i).queue = [];
end
path = FloydShortestPath(roadArray, length(nodeArray));
% data recording
trafficRecord = InitTrafficRecord();
carRecord = cell(1, clockMax);
% random number generator
rng(1);

%% Simulation
for clock = 1:clockMax
    % generate cars / set destination to idle cars
    GenerateNewRoutes();
    % update road.queue to get new cars on road
    UpdateRoadQueue();
    % update velocity of cars currently on road
    UpdateVelocity();
    % move cars & delete arrived cars
    MoveCar();

    % data recording
    carRecord{clock} = carArray;
end
  
trafficRecord.carRecord = carRecord;
filename = string(datestr(now,'yyyy-mm-dd-HH-MM')) + ".mat";
save(filename, 'trafficRecord');
disp("Result saved to file "+filename);
clear TrafficSimulation

%==========================================================================

%% nested function: do some basic check of whether input variables are valid
function [] = checkInput()
    
    % check imageRoad
    for ii = 1:length(o.roadArray)
        if (o.roadArray(ii).imageRoad > 0)
            j = o.roadArray(ii).imageRoad;
            if (o.roadArray(ii).nodeStart ~= o.roadArray(j).nodeEnd ||...
                o.roadArray(ii).nodeEnd ~= o.roadArray(j).nodeStart ||...
                o.roadArray(ii).length ~= o.roadArray(j).length)
                error("roadArray input invalid: 'imageRoad'.");
            end
        end
    end
    % check nodeArray
    for ii = 1:length(o.nodeArray)
        if (o.nodeArray(ii).spawnRate < 0 || o.nodeArray(ii).spawnRate > 1)
            error("nodeArray input invalid: 'spawnRate'.");
        end
        if (o.nodeArray(ii).destChance < 0)
            error("nodeArray input invalid: 'destChance'.");
        end
    end
    % check roadArray
    for ii = 1:length(o.roadArray)
        if (o.roadArray(ii).length < o.dmin * 2)
            error("roadArray input invalid: 'length' too short compared with 'dmin'.");
        end
    end
end


%% nested function: initialize traffic record object
function [trafficRecord_] = InitTrafficRecord()
    trafficRecord_.vmax = vmax;
    trafficRecord_.dmin = dmin;
    trafficRecord_.dmax = dmax;
    trafficRecord_.simulationTime = simulationTime;
    trafficRecord_.dt = dt;
    trafficRecord_.clockMax = clockMax;
    trafficRecord_.roadArray = roadArray;
    trafficRecord_.nodeArray = nodeArray;
    if (isfield(o, 'visArray') && length(o.visArray)==length(nodeArray))
        trafficRecord_.visArray = o.visArray;
    end
end


%% nested function: at the beginning of every tick, generate new routes. If
% there are idle cars, use them first, otherwise create new cars.
function [] = GenerateNewRoutes()
    
    % Generate routes for every node
    for ii = 1:length(nodeArray)
        targetNumCars = nodeArray(ii).spawnRate * dt; % note: this number can be smaller than 1
        % first generate integer part of 'targetNumCars'
        while (targetNumCars >= 1)
            targetNumCars = targetNumCars - 1;
            GenerateOneNewRoute(ii);
        end
        % then generate decimal part of 'targetNumCars'
        randNum = rand(1);
        if (randNum < targetNumCars)
            GenerateOneNewRoute(ii);
        end
    end
    
    
    %% double nested function: Generate one route from 'origin'
    function [] = GenerateOneNewRoute(origin)
        
        persistent destChanceSum
    
        % Calculate 'destChanceSum' if this is the first call of this
        % function
        if isempty(destChanceSum)
            destChanceSum = 0;
            for iii = 1:length(nodeArray)
                destChanceSum = destChanceSum + nodeArray(iii).destChance;
            end
        end
    
        localDestChanceSum = destChanceSum - nodeArray(origin).destChance;
        sum = 0;
        randNumDest = rand(1);
        for iii = 1:length(nodeArray)
            if (iii == origin), continue;end
            sum = sum + nodeArray(iii).destChance;
            if (sum / localDestChanceSum > randNumDest)
                % the iii-th node is chosen to be destination
                carIndex = GetIdleCar();
                carPath = GenerateRoadArray(origin, iii);
                carArray(carIndex).roadNum = carPath(1);
                carArray(carIndex).alive = 2;
                carArray(carIndex).roadArray = carPath;
                % append to the waiting queue of this road
                roadIdx = carPath(1);
                roadArray(roadIdx).queue(end+1) = carIndex;
                % do not generate another route
                break;
            end
        end
    end
    
    
    %% double nested function: Get the index of an idle car. Generate a new
    %  car if there is none.
    function [res] = GetIdleCar()
        for iii = 1:length(carArray)
            if (carArray(iii).alive == 0)
                res = iii;
                return;
            end
        end
        % no idle car left, generate some new cars
        num = length(carArray);
        carArray(num+1:num*2) = carArray(1:num);  % double the size every time
        for iii = num+1:num*2
            carArray(iii).velocity = 0;
            carArray(iii).roadNum = 0;
            carArray(iii).position = 0;
            carArray(iii).frontCar = 0;
            carArray(iii).backCar = 0;
            carArray(iii).alive = 0;
            carArray(iii).roadArray = [];
        end

        res = num + 1;
    end
    
    
    %% Generate an array of path (represented by road number) from 'origin' 
    %  node and 'dest' node
    function [res] = GenerateRoadArray(origin, dest)
        if (path(origin, dest) == -1)
            res = GetRoadIndex(origin, dest);
            return;
        else
            res = GenerateRoadArray(origin, path(origin, dest));
            res = [res, GenerateRoadArray(path(origin, dest), dest)];
        end
    end
end


%% nested function: get road index from 'nodeStart' and 'nodeEnd'
function [res] = GetRoadIndex(nodeStart, nodeEnd)
    
    for ii = 1:length(roadArray)
        if (roadArray(ii).nodeStart == nodeStart && roadArray(ii).nodeEnd == nodeEnd)
            res = ii;
            return;
        end
    end
    res = -1;
end


%% nested function: if the road is not full, let go one car from queue
function [] = UpdateRoadQueue()

    for ii = 1:length(roadArray)
        if isempty(roadArray(ii).queue)
            continue;
        end
        queueCarIdx = roadArray(ii).queue(1);
        backCarIdx = roadArray(ii).backCar;
        if (backCarIdx > 0 &&...
            carArray(backCarIdx).position * roadArray(ii).length < dmin)
            % this road is full now
            continue;
        else
            % good to go!
            carArray(queueCarIdx).frontCar = backCarIdx;
            carArray(queueCarIdx).alive = 1;
            if (backCarIdx > 0)
                carArray(backCarIdx).backCar = queueCarIdx;
            end
            roadArray(ii).backCar = queueCarIdx;
            if (roadArray(ii).frontCar == 0)
                roadArray(ii).frontCar = queueCarIdx;
            end
            roadArray(ii).queue(1) = [];  % pop from queue
        end
    end
end


%% nested function: update the velocity of cars on roads
function [] = UpdateVelocity()
    
    % update velocity
    for ii = 1:length(roadArray)
        carIdx = roadArray(ii).backCar;
        while (carIdx > 0)
            frontCarIdx = carArray(carIdx).frontCar;
            if (frontCarIdx == 0)
                % no car ahead (on current road)
                if (length(carArray(carIdx).roadArray) == 1)  % if near destination
                    carArray(carIdx).velocity = Velocity((dmax + dmin)/2);
                else
                    nextRoadIdx = carArray(carIdx).roadArray(2);
                    nextCarIdx = roadArray(nextRoadIdx).backCar;
                    if (nextCarIdx == 0)  
                        % if no car ahead (even taking next road into consideration)
                        carArray(carIdx).velocity = Velocity((dmax + dmin)/2);
                    else
                        % calculate the dist to the car on next road
                        currentRoadIdx = carArray(carIdx).roadArray(1);
                        dist = (1 - carArray(carIdx).position) * roadArray(currentRoadIdx).length;
                        dist = dist + carArray(nextCarIdx).position * roadArray(nextRoadIdx).length;
                        carArray(carIdx).velocity = Velocity(dist);
                    end
                end
            else
                % car ahead
                dist = roadArray(ii).length * (carArray(frontCarIdx).position - carArray(carIdx).position);
                carArray(carIdx).velocity = Velocity(dist);
            end
            carIdx = frontCarIdx;
        end
    end
    
    
    %% double nested function: Calculate velocity from distance to front car
    function [res] = Velocity(d)

        if(d < dmin)
            res=0;
        elseif(d < dmax)
            res = vmax * log(d/dmin) / log(dmax/dmin);
        else
            res = vmax;
        end
    end
end


%% nested function: move cars & delete arrived cars
function [] = MoveCar()
    
    for ii = 1:length(carArray)
        if (carArray(ii).alive ~= 1), continue;end
        deltaDist = carArray(ii).velocity * dt;
        roadLength = roadArray(carArray(ii).roadNum).length;
        carArray(ii).position = carArray(ii).position + deltaDist / roadLength;
        
        if (carArray(ii).position > 1)  % end of current road
            
            carArray(ii).roadArray(1) = [];
            if isempty(carArray(ii).roadArray)  % destination arrived
                currentRoadIdx = carArray(ii).roadNum;
                backCarIdx = carArray(ii).backCar;
                % update car
                DeleteCar(ii);
                % update back car
                if (backCarIdx > 0)
                    carArray(backCarIdx).frontCar = 0;
                end
                % update road
                roadArray(currentRoadIdx).frontCar = backCarIdx;
                if (roadArray(currentRoadIdx).backCar == ii)
                    roadArray(currentRoadIdx).backCar = 0;
                end
            else  % move to next road
                currentRoadIdx = carArray(ii).roadNum;
                nextRoadIdx = carArray(ii).roadArray(1);
                backCarIdx = carArray(ii).backCar;
                frontCarIdx = roadArray(nextRoadIdx).backCar;
                % update this car
                carArray(ii).roadNum = nextRoadIdx;
                carArray(ii).position = (carArray(ii).position - 1) * roadArray(currentRoadIdx).length / roadArray(nextRoadIdx).length;
                carArray(ii).frontCar = frontCarIdx;
                carArray(ii).backCar = 0;
                % update back car
                if (backCarIdx > 0)
                    carArray(backCarIdx).frontCar = 0;
                end
                % update front car
                if (frontCarIdx > 0)
                    carArray(frontCarIdx).backCar = ii;
                end
                % update current road
                roadArray(currentRoadIdx).frontCar = backCarIdx;
                if (roadArray(currentRoadIdx).backCar == ii)
                    roadArray(currentRoadIdx).backCar = 0;
                end
                % update next road
                roadArray(nextRoadIdx).backCar = ii;
                if (roadArray(nextRoadIdx).frontCar == 0)
                    roadArray(nextRoadIdx).frontCar = ii;
                end
            end
        end
    end
    
    
    % double nested function: inactivate the car once it arrives at
    % destination
    function [] = DeleteCar(idx)
        
        carArray(idx).velocity = 0;
        carArray(idx).roadNum = 0;
        carArray(idx).position = 0;
        carArray(idx).frontCar = 0;
        carArray(idx).backCar = 0;
        carArray(idx).alive = 0;
        carArray(idx).roadArray = [];
    end
end
end

% END of TrafficSimulation.m
