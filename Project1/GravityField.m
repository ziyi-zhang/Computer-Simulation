function [] = GravityField(filename, focusIdx, fastForward)
% Visualize the gravity field based on the stored data
% 'filename' should be a 'mat' file created by 'SaveData.m'

    if nargin<2
        if strcmp(filename, 'harmonic'), focusIdx=2;fastForward=4;end
        if strcmp(filename, 'tribody'), focusIdx=1;fastForward=4;end
        if strcmp(filename, 'square'), focusIdx=1;fastForward=4;end
        if strcmp(filename, 'stable'), focusIdx=1;fastForward=30;end
    end
    %% read in data
    fprintf('Loading file <%s> with focus point index %d.\n', filename, focusIdx);
    data = load(filename);
    data = data.data;
    M = data.M;
    L = data.L;
    V = data.V;
    A = data.A;
    N = length(M);
    Spring = data.Spring;
    colorSprIdx = data.colorSprIdx;
    tmax = data.tmax;
    clockmax = data.clockmax;
    dt = tmax / clockmax;
    if ~isfield(data, 'XHist') || ~isfield(data, 'YHist')
        warning('Please run simulation before running this script.');
        return;
    end
    XHist = data.XHist;
    YHist = data.YHist;
    % fastForward
    idx = 1:fastForward:clockmax;
    clockmax = length(idx);
    XHist = XHist(:, idx);
    YHist = YHist(:, idx);
    colorSprIdx = colorSprIdx(:, idx);
    
    %% visRealtime
    % marker color
    cMarker = colorcube;
    s = RandStream('mt19937ar','Seed',0);
    RandStream.setGlobalStream(s);
    idx = randperm(size(cMarker, 1));
    cMarker = cMarker(idx, :);
    cMarker = [lines(7); cMarker];
    % spring color
    cSpringConst = 55;
    cSpring = whitejet(cSpringConst);
    % size
    mSize = calcSize(M);
    % trail deque
    dequeSize = 600;
    if strcmp(filename, 'harmonic')
        dequeSize = 1;
    elseif strcmp(filename, 'tribody')
        dequeSize = 200;
    elseif strcmp(filename, 'stable')
        dequeSize = 200;
    end
    Xtrail = ones(N, dequeSize) .* L(1, :)';
    Ytrail = ones(N, dequeSize) .* L(2, :)';
    % figrue
    figure('Position', [1739 50 1600 900])
    set(gcf, 'double', 'on');
    leftPanel = subplot(1, 2, 1);
    hold on
    % plot anchor and normal points
    for i = 1:N
        if A(i)
            h(i) = scatter(leftPanel, L(1, i), L(2, i), mSize(i)+10, 's', 'MarkerFaceColor', cMarker(i, :),...
                        'MarkerEdgeColor', 'black', 'LineWidth', 2);
        else
            h(i) = scatter(leftPanel, L(1, i), L(2, i), mSize(i), 'o', 'MarkerFaceColor', cMarker(i, :), 'MarkerEdgeColor', 'none');
            htrail(i) = plot(leftPanel, Xtrail(i, :), Ytrail(i, :), 'Color', cMarker(i, :));
        end
    end
    % plot springs
    for i = 1:size(Spring, 1)
        hspr(i) = plot(leftPanel, L(1, [Spring(i, 1), Spring(i, 2)]), L(2, [Spring(i, 1), Spring(i, 2)]),...
                      '-.', 'Marker', 'none', 'LineWidth', 1.5);
    end
    axis equal
    hold off
    % color bar
    cmap = colormap(leftPanel, cSpring);
    caxis(leftPanel, [-9 9]);
    caxis(leftPanel, 'manual');
    cbar = colorbar();
    cbar.Label.String = 'String force from compression to extension';
    cbar.Label.FontSize = 15;
    % axis lim
    halfLength = max([max(XHist, [], 'all')-min(XHist, [], 'all'), max(YHist, [], 'all')-min(YHist, [], 'all')]);
    if halfLength<5, halfLength=5;end
    axisXMin = -halfLength*1.1;
    axisXMax = halfLength*1.1;
    axisYMin = -halfLength*1.1;
    axisYMax = halfLength*1.1;
    if strcmp(filename, 'harmonic')
        axisXMin=-2;axisXMax=20;axisYMin=-7;axisYMax=7;
    elseif strcmp(filename, 'tribody')
        axisXMin=-10;axisXMax=17;axisYMin=-8;axisYMax=8;
    elseif strcmp(filename, 'square')
        axisXMin=-20;axisXMax=20;axisYMin=-20;axisYMax=20;
    elseif strcmp(filename, 'stable')
        axisXMin=-22;axisXMax=22;axisYMin=-22;axisYMax=22;
    end
    xlim(leftPanel, [axisXMin, axisXMax]);
    ylim(leftPanel, [axisYMin, axisYMax]);

    %% Gravity Field
    % springFocus to accelerate
    SpringFocus = [];
    SpringFocusConnected = false(N, 1);
    SpringFocusConnected(focusIdx) = true;
    for i = 1:size(Spring, 1)
        if Spring(i, 1)==focusIdx
            SpringFocus = [SpringFocus; Spring(i, 2:4)];  %#ok
            SpringFocusConnected(Spring(i, 2)) = true;
        end
        if Spring(i, 2)==focusIdx
            SpringFocus = [SpringFocus; Spring(i, 1), Spring(i, 3:4)];  %#ok
            SpringFocusConnected(Spring(i, 1)) = true;
        end
    end
    % plot hfieldmass
    fieldResolution = 50;
    axisXSpace = (axisXMax - axisXMin) / (fieldResolution-1);
    axisYSpace = (axisYMax - axisYMin) / (fieldResolution-1);
    [surfX, surfY] = meshgrid(axisXMin:axisXSpace:axisXMax, axisYMax:-axisYSpace:axisYMin);
    rightPanel = subplot(1, 2, 2);
    hold on
    hfield = surf(rightPanel, surfX, surfY, zeros(length(surfX), length(surfY)), 'FaceColor', 'flat', 'EdgeAlpha', 0.1);
    for i = 1:N
        if ~SpringFocusConnected(i), continue;end  % only plot connected masses
        hfieldmass(i) = scatter3(rightPanel, L(1, i), L(2, i), 0, 40, cMarker(i, :), 'filled', 'MarkerEdgeColor', 'b');
    end
    % legend
    fakeTarget = scatter3(rightPanel, 0, 0, -200, 2, cMarker(focusIdx, :), 'filled');
    legend_ = legend(fakeTarget, 'Surface reconstructed using this mass');
    legend_.FontSize = 15;
    legend_.AutoUpdate = false;
    % done
    grid on
    hold off
    view(3)
    force = zeros(fieldResolution, fieldResolution, 2);
    XArr = axisXMin:axisXSpace:axisXMax;
    YArr = axisYMax:-axisYSpace:axisYMin;
    xlim(rightPanel, [axisXMin, axisXMax]);
    ylim(rightPanel, [axisYMin, axisYMax]);
    zlimMin = -100;
    zlimMax = 100;
    if strcmp(filename, 'harmonic')
        zlimMin = -30;zlimMax = 40;
    elseif strcmp(filename, 'tribody')
        zlimMin = -100;zlimMax = 150;
    elseif strcmp(filename, 'square')
        zlimMin = -150;zlimMax = 50;
    elseif strcmp(filename, 'stable')
        zlimMin = -3e5;zlimMax = 4.5e5;
    end
    zlim(rightPanel, [zlimMin, zlimMax]);
    massZCorrection = 0.01 * (zlimMax - zlimMin);

    % pauseTime/FPS control
    pauseTime = 0.0;

    %% Start Animation
    for t = 1:clockmax
        pause(pauseTime);
        %% Left Panel
        % update points
        for i = 1:N
            if A(i), continue;end
            set(h(i), 'xdata', XHist(i, t), 'ydata', YHist(i, t));
        end

        for i = 1:N
            if A(i), continue;end
            idxStart = t-dequeSize;
            if idxStart<1, idxStart=1;end
            set(htrail(i), 'xdata', XHist(i, idxStart:t),...
                           'ydata', YHist(i, idxStart:t));
        end
        % update spring, update spring color
        for i = 1:size(Spring, 1)
            set(hspr(i), 'xdata', XHist([Spring(i, 1), Spring(i, 2)], t),...
                         'ydata', YHist([Spring(i, 1), Spring(i, 2)], t),...
                         'Color', cSpring(colorSprIdx(i, t), :));
        end
        %% Right Panel
        force = zeros(fieldResolution, fieldResolution, 2);
        for i = 1:length(YArr)
            for j = 1:length(XArr)
                for k = 1:size(SpringFocus, 1)

                    p = [XArr(j); YArr(i)];  % position of probe
                    k_ = SpringFocus(k, 1);  % idx of connected point
                    Lk_ = [XHist(k_, t); YHist(k_, t)];  % location of k_
                    r = norm([XHist(k_, t); YHist(k_, t)] - p);
                    forceElement = SpringFocus(k, 3) .* (r - SpringFocus(k, 2)) .* (Lk_ - p) ./ r;
                    force(i, j, 1) = force(i, j, 1) + forceElement(1);
                    force(i, j, 2) = force(i, j, 2) + forceElement(2);
                end
            end
        end
        force = toSmooth(force);
        Z = -g2s(force(:,:,1), force(:,:,2), XArr', YArr');

        Z(Z>zlimMax) = zlimMax;
        Z(Z<zlimMin) = zlimMin;
        
        set(hfield, 'zdata', Z); % update field
        massZ = interp2(surfX, surfY, Z, XHist(:, t), YHist(:, t));
        for i = 1:N
            if ~SpringFocusConnected(i), continue;end
            set(hfieldmass(i), 'xdata', XHist(i, t), 'ydata', YHist(i, t), 'zdata', massZ(i)+massZCorrection);
        end
        
        %% update figure
        drawnow
    end
end


%%
function [mSize] = calcSize(M)

    low = 32;
    high = 90;
    minM = min(M);
    mSize = low .* (M ./ minM);
    mSize(mSize > high) = high;
end


function [force] = toSmooth(force)
    
    forceAbs = abs(force);
    forceMax = quantile(forceAbs(:), 0.98);
    force(force>forceMax) = forceMax;
    force(force<-forceMax) = -forceMax;
end


function [massZ] = calcZ(Z, L, XArr, YArr)
    
    massZ = zeros(1, size(L, 2));
    for i = 1:size(L, 2)
        xx = 1;
        yy = 1;
        for j = 1:length(XArr)
            if XArr(j)>L(1, i)
               xx = j;
               break;
            end
        end
        for j = 1:length(YArr)
            if YArr(j)<L(2, i)
               yy = j;
               break;
            end
        end
        massZ(i) = Z(yy, xx) + 0.2;
    end
end
