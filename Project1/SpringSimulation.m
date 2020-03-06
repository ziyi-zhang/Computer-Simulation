%% spring system simulation
% The main script to run simulation
% Use 'GravityField' to visualize the gravity field after running this
% simulation

%{
% visualization parameters
c = [0, 0, 1;
     0.4, 0.8, 0.5;
     0.8, 0, 0.5;
     0.9, 0.4, 0];
dequeSize = 150;
XHist = ones(N, dequeSize) .* L(1, :)';
YHist = ones(N, dequeSize) .* L(2, :)';
dequePtr = 1; % deque to show trace
%{

%}
axisXMin = -1.5 * max(L, [], 'all');
axisXMax =  1.5 * max(L, [], 'all');
axisYMin = -1.5 * max(L, [], 'all');
axisYMax =  1.5 * max(L, [], 'all');
fieldResolution = 50;
axisXSpace = (axisXMax - axisXMin) / (fieldResolution-1);
axisYSpace = (axisYMax - axisYMin) / (fieldResolution-1);
[surfX, surfY] = meshgrid(axisXMin:axisXSpace:axisXMax, axisYMax:-axisYSpace:axisYMin);
% graphics output
set(gcf, 'double', 'on')
subplot(1, 2, 1)
hold on
h = scatter(L(1, :), L(2, :), 25.*ones(1, N), c(1:N, :), 'filled');
for i = 1:N
    htrail(i) = plot(XHist(i, :), YHist(i, :), 'Color', c(i, :));
end
hold off
subplot(1, 2, 2)
hold on
hfield = surf(surfX, surfY, zeros(length(surfX), length(surfY)), 'FaceColor', 'flat', 'EdgeAlpha', 0.1);
for i = 1:N
    hfieldmass(i) = scatter3(L(1, i), L(2, i), 0, 40, c(i, :), 'filled', 'MarkerEdgeColor', 'b');
end
hold off
view(3)
%}
%
%% Check validity of input
A = logical(A);
assert(nnz(V(:, A)) == 0);  % anchor points must have velocity zero
assert(size(Spring, 2) == 4);
assert(nnz(M<=0) == 0);
for i = 1:size(Spring, 1)
    assert(~( A(Spring(i, 1)) && A(Spring(i, 2)) ));  % should not have string connecting two anchor points
end

%% visualization
if visRealtime
    
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
    Xtrail = ones(N, dequeSize) .* L(1, :)';
    Ytrail = ones(N, dequeSize) .* L(2, :)';
    dequePtr = 1; % deque to show trace
    % figrue
    figure('Position', [1739 50 1600 900])
    hold on
    set(gcf, 'double', 'on');
    % plot anchor and normal points
    for i = 1:N
        if A(i)
            h(i) = scatter(L(1, i), L(2, i), mSize(i)+10, 's', 'MarkerFaceColor', cMarker(i, :),...
                        'MarkerEdgeColor', 'black', 'LineWidth', 2);
        else
            h(i) = scatter(L(1, i), L(2, i), mSize(i), 'o', 'MarkerFaceColor', cMarker(i, :), 'MarkerEdgeColor', 'none');
            htrail(i) = plot(Xtrail(i, :), Ytrail(i, :), 'Color', cMarker(i, :));
        end
    end
    % plot springs
    for i = 1:size(Spring, 1)
        hspr(i) = plot(L(1, [Spring(i, 1), Spring(i, 2)]), L(2, [Spring(i, 1), Spring(i, 2)]),...
                      '-.', 'Marker', 'none', 'LineWidth', 1.5);
    end
    axis equal
    hold off
    % color bar
    %{
    colormap(cSpring);
    caxis([-9 9]);
    cbar = colorbar;
    cbar.Label.String = 'String force from compression to extension';
    %}
    % axis lim
    halfLength = max([max(L(1, :))-min(L(1, :)), max(L(2, :))-min(L(2, :))]);
    if halfLength<5, halfLength=5;end
    axisXMin = -halfLength*1.5;
    axisXMax = halfLength*1.5;
    axisYMin = -halfLength*1.5;
    axisYMax = halfLength*1.5;
end

if storeData
    
    XHist = zeros(N, clockmax);
    YHist = zeros(N, clockmax);
    colorSprIdx = zeros(size(Spring, 1), clockmax);
end

%% simulation starts here
count = 1;
for clock = 1:clockmax

    t = clock * dt;
    F = calcF(L, Spring);
    % now updating the velocity of i-th object
    for i = 1:N

        if A(i), continue;end  % skip anchor
        a = zeros(2, 1);
        for j = 1:N
            if i==j, continue;end  % divide by 0
            a = a + F(i, j) .* (L(:, j)-L(:, i)) ./ norm(L(:, i)-L(:, j));
        end
        a = a ./ M(i);

        % update velocity
        V(:, i) = V(:, i) + a .* dt;
    end
    % updating the location of i-th object
    L = L + V .* dt;

    % update realtime figrue
    if visRealtime
        for i = 1:N
            if A(i), continue;end
            set(h(i), 'xdata', L(1, i), 'ydata', L(2, i));
        end

        Xtrail(:, dequePtr) = L(1, :)';
        Ytrail(:, dequePtr) = L(2, :)';
        dequePtr = dequePtr + 1;
        if (dequePtr > dequeSize), dequePtr = 1;end
        for i = 1:N
            if A(i), continue;end
            set(htrail(i), 'xdata', [Xtrail(i, dequePtr:dequeSize), Xtrail(i, 1:dequePtr-1)],...
                           'ydata', [Ytrail(i, dequePtr:dequeSize), Ytrail(i, 1:dequePtr-1)]);
        end
    end
    % update spring, update spring color
    for i = 1:size(Spring, 1)
        dev = ceil(abs(3*F(Spring(i, 1), Spring(i, 2))));
        if dev>27, dev=27;end
        if F(Spring(i, 1), Spring(i, 2))>0
            color = cSpring(28+dev, :);
            if storeData, colorSprIdx(i, clock)=28+dev;end
        else
            color = cSpring(28-dev, :);
            if storeData, colorSprIdx(i, clock)=28-dev;end
        end
        if visRealtime
            set(hspr(i), 'xdata', L(1, [Spring(i, 1), Spring(i, 2)]),...
                         'ydata', L(2, [Spring(i, 1), Spring(i, 2)]),...
                         'Color', color);
        end
    end
    if visRealtime
        % drawnow
        drawnow
        xlim([axisXMin, axisXMax]);
        ylim([axisYMin, axisYMax]);
        pause(pauseTime);
    end
    % store data
    if storeData
        XHist(:, clock) = L(1, :)';
        YHist(:, clock) = L(2, :)';
    end
    % replot subimage-2
    %{
    gradient_ = zeros(fieldResolution, fieldResolution, 2);
    XArr = axisXMin:axisXSpace:axisXMax;
    YArr = axisYMax:-axisYSpace:axisYMin;
    for i = 1:length(XArr)
        for j = 1:length(YArr)
            for k = 1:N % for point (i, j) calc k-th planet's force

                p = [XArr(j); YArr(i)]; % position of probe
                r = norm(L(:, k) - p);
                gradient_(i, j, 1) = gradient_(i, j, 1) + 1 * (r - R0_) * (L(1, k) - p(1)) / r;
                gradient_(i, j, 2) = gradient_(i, j, 2) + 1 * (r - R0_) * (L(2, k) - p(2)) / r;
            end
        end
    end
    Z = -g2s(gradient_(:,:,1), gradient_(:,:,2), XArr', YArr');
    ZMax = quantile(Z(:), 0.95);
    Z(Z>ZMax) = ZMax;
    ZMin = quantile(Z(:), 0.05);
    Z(Z<ZMin) = ZMin;
    set(hfield, 'zdata', Z); % update field
    massZ = calcZ(Z, L, XArr, YArr);
    for i = 1:N
        set(hfieldmass(i), 'xdata', L(1, i), 'ydata', L(2, i), 'zdata', massZ(i));
    end
    %}
    if ~visRealtime
        if count > 500
            fprintf('%.2f\n', clock/clockmax);
            count = 0;
        end
        count = count + 1;
    end
end

%%
function [F] = calcF(L, Spring)

    N = size(L, 2);
    F = zeros(N, N);
    for i = 1:size(Spring, 1)
        f = Spring(i, 4) * (norm(L(:, Spring(i, 1))-L(:, Spring(i, 2))) - Spring(i, 3));
        F(Spring(i, 1), Spring(i, 2)) = f;
        F(Spring(i, 2), Spring(i, 1)) = f;
    end
end


function [mSize] = calcSize(M)

    low = 32;
    high = inf;
    minM = min(M);
    mSize = low .* (M ./ minM);
    mSize(mSize > high) = high;
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
