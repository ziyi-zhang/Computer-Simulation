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
nodeArray = trafficRecord.nodeArray;
carRecord = trafficRecord.carRecord;
visArray = trafficRecord.visArray;
assert(length(visArray) == length(nodeArray));

%% Fast forward
carRecord = carRecord(1:fastForward:end);

%% Draw background (road & node)
% everything drawn in [-1, 1]*[-1, 1]
figure
axis manual
axis([-1 1 -1 1]);
pbaspect([1 1 1]);
hold on
title("Traffic Simulation Visualization");

nodeRadius = 0.03;
roadWidth = nodeRadius * 0.5;
% roads
for i = 1:length(roadArray)
    if (roadArray(i).imageRoad == 0)
       % one-way road
       DrawOnewayRoad(i);
    else
       % two-way road
       DrawTwowayRoad(i);
    end
end
% nodes
for i = 1:length(visArray)
    coord = [visArray(i).x-nodeRadius, visArray(i).y-nodeRadius, 2*nodeRadius, 2*nodeRadius];
    rectangle('Position', coord, 'Curvature', [1, 1], 'FaceColor', 'white', 'EdgeColor', 'black');
end


%% Start Animation
for clock = 1:length(carRecord)
    
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
    colormap(copper);
    patch([x0+orthVec(1), x1+orthVec(1)], [y0+orthVec(2), y1+orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x0-orthVec(1), x1-orthVec(1)], [y0-orthVec(2), y1-orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
end


%% nested function: draw two-way road
function [] = DrawTwowayRoad(idx)

    x0 = visArray(roadArray(idx).nodeStart).x;
    y0 = visArray(roadArray(idx).nodeStart).y;
    x1 = visArray(roadArray(idx).nodeEnd).x;
    y1 = visArray(roadArray(idx).nodeEnd).y;
    dist = norm([x0-x1, y0-y1]);
    orthVec = [y0-y1, x1-x0]/dist;  % unit vector, anti-clockwise 90 degrees + from start to end
    orthVec = orthVec .* roadWidth ./ 2;
    colormap(copper);
    patch([x0-orthVec(1), x1-orthVec(1)], [y0-orthVec(2), y1-orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x0-orthVec(1)*3, x1-orthVec(1)*3], [y0-orthVec(2)*3, y1-orthVec(2)*3], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x1+orthVec(1), x0+orthVec(1)], [y1+orthVec(2), y0+orthVec(2)], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
    patch([x1+orthVec(1)*3, x0+orthVec(1)*3], [y1+orthVec(2)*3, y0+orthVec(2)*3], [0.1, 0.9], 'FaceColor','none','EdgeColor','interp');
end
end

% END of TrafficVisualization.m
