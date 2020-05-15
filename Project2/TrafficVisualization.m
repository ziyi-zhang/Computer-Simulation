% Main visualization file
% Ziyi. April 2020.
function [] = TrafficVisualization(trafficRecord, fastForward)

%% Phase input variables
if (nargin<2), fastForward=1;end
if (~isfield(trafficRecord, 'visArray'))
    error("TrafficVisualization needs a visArray to plot the figure. It might be missing or invalid.");
end
vmax = trafficRecord.vmax;
dmin = trafficRecord.dmin;
dmax = trafficRecord.dmax;
simulationTime = trafficRecord.simulationTime;
dt = trafficRecord.dt;
clockMax = trafficRecord.clockMax;
roadArray = trafficRecord.roadArray;
roadArray(1).xyStart = [];  % Add two new fields to 'roadArray'
roadArray(1).unitVec = [];
nodeArray = trafficRecord.nodeArray;
carRecord = trafficRecord.carRecord;
visArray = trafficRecord.visArray;
assert(length(visArray) == length(nodeArray));

%% Fast forward
carRecord = carRecord(1:fastForward:end);

%% Draw background (road & node)
% everything drawn in [-1, 1]*[-1, 1]
figure
set(gcf, 'Position',  [200, 200, 1000, 1000])
axis manual
axis([-1 1 -1 1]);
pbaspect([1 1 1]);
hold on
title("Traffic Simulation Visualization");

nodeRadius = 0.03;
roadWidth = nodeRadius * 0.7;
% roads
for i = 1:length(roadArray)
    if (roadArray(i).imageRoad == 0)
       % one-way road
       DrawOnewayRoad(i);
    else
       if (roadArray(i).imageRoad < i), continue;end
       % two-way road
       DrawTwowayRoad(i);
    end
end
% nodes
for i = 1:length(visArray)
    coord = [visArray(i).x-nodeRadius, visArray(i).y-nodeRadius, 2*nodeRadius, 2*nodeRadius];
    rectangle('Position', coord, 'Curvature', [1, 1], 'FaceColor', 'white', 'EdgeColor', 'black');
end


%% Prepare handles of cars
carSize = 10;
for i = 1:length(roadArray)
    hcar(i) = scatter(0, 0, carSize, 'black', 'filled');
end


%% Start Animation
pauseTime = 0.05;  % pauseTime/FPS control
if (fastForward>1), pauseTime=0;end  % do not pause if fast forward
for clock = 1:length(carRecord)
    
    pause(pauseTime);
    carArray = carRecord{clock};
    roadCar = carArray([carArray.alive] == 1);
    queueCar = carArray([carArray.alive] == 2);
    % plot cars on road (alive == 1)
    for i = 1:length(roadArray)
        cars = roadCar([roadCar.roadNum] == i);
        if isempty(cars)
            set(hcar(i), 'xdata', [], 'ydata', []);
            continue;
        end
        loc = repmat(roadArray(i).xyStart, length(cars), 1);
        loc = loc + roadArray(i).unitVec .* [cars.position]';
        velocity = [cars.velocity];
        % plot
        set(hcar(i), 'xdata', loc(:, 1), 'ydata', loc(:, 2), 'cdata', Velocity2Color(velocity));
    end
    % plot cars in queue (alive == 2)
    for i = 1:length(roadArray)
        
    end
    
    % update figure
    drawnow
end


%% nested function: draw one-way road
function [] = DrawOnewayRoad(idx)

    x0 = visArray(roadArray(idx).nodeStart).x;
    y0 = visArray(roadArray(idx).nodeStart).y;
    x1 = visArray(roadArray(idx).nodeEnd).x;
    y1 = visArray(roadArray(idx).nodeEnd).y;
    dist = norm([x0-x1, y0-y1]);
    orthVec = [y0-y1, x1-x0]/dist;  % unit vector, anti-clockwise 90 degrees + from start to end
    orthVec = orthVec .* roadWidth ./ 2;
    copperColormap = copper(100);
    colormap(copperColormap(1:70, :));
    patch([x0+orthVec(1), x1+orthVec(1)], [y0+orthVec(2), y1+orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x0-orthVec(1), x1-orthVec(1)], [y0-orthVec(2), y1-orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    % set new fields
    roadArray(idx).xyStart = [x0, y0];
    roadArray(idx).unitVec = [y0-y1, x1-x0]/dist;
end


%% nested function: draw two-way road
function [] = DrawTwowayRoad(idx)

    x0 = visArray(roadArray(idx).nodeStart).x;
    y0 = visArray(roadArray(idx).nodeStart).y;
    x1 = visArray(roadArray(idx).nodeEnd).x;
    y1 = visArray(roadArray(idx).nodeEnd).y;
    dist = norm([x0-x1, y0-y1]);
    orthVec = [y0-y1, x1-x0]/dist;  % unit vector, anti-clockwise 90 degrees + from start to end
    dispVec = orthVec .* roadWidth .* 0.2;
    orthVec = orthVec .* roadWidth;
    copperColormap = copper(100);
    colormap(copperColormap(1:70, :));
    patch([x0-dispVec(1), x1-dispVec(1)], [y0-dispVec(2), y1-dispVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x0-orthVec(1)-dispVec(1), x1-orthVec(1)-dispVec(1)], [y0-orthVec(2)-dispVec(2), y1-orthVec(2)-dispVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x1+dispVec(1), x0+dispVec(1)], [y1+dispVec(2), y0+dispVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x1+orthVec(1)+dispVec(1), x0+orthVec(1)+dispVec(1)], [y1+orthVec(2)+dispVec(2), y0+orthVec(2)+dispVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    % set new fields
    roadArray(idx).xyStart = [x0, y0] - dispVec - 0.5 .* orthVec;
    roadArray(idx).unitVec = [x1-x0, y1-y0]/dist;
    imageIdx = roadArray(idx).imageRoad;
    roadArray(imageIdx).xyStart = [x1, y1] + dispVec + 0.5 .* orthVec;
    roadArray(imageIdx).unitVec = -[x1-x0, y1-y0]/dist;
end


%% nested function: return a color based on speed
function [res] = Velocity2Color(speed)

    persistent carColor gap
    
    if isempty(carColor)
        carColor = parula(101);
        gap = vmax / 100;
    end

    res = carColor(ceil(speed / gap), :);
end
end

% END of TrafficVisualization.m
