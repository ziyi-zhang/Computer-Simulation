% Main Simulation File
% Ziyi. May, 2020.
function [] = CellularAutomata3D(o)

%% Phase input variables
size_ = o.size;
init = o.init;
rule = o.rule;
assert(length(rule) == 64);
assert(nnz(rule == 0 | rule == 1) == 64);
MaxIteration = o.it;

%% Initialization
rule(1) = 0;  % if no cell around, stay deactivated
grid_ = zeros(size_, size_, size_);
state = zeros(size_, size_, size_);
for i = 1:size(init, 2)
    Activate(init(1, i), init(2, i), init(3, i));
    state(init(1, i), init(2, i), init(3, i)) = 1;
end

%% Plot
figure
hold on
h = scatter3(1, 1, 1, 'filled');
axis([0 size_+1 0 size_+1 0 size_+1]);
pbaspect([1 1 1]);
view(3);
grid on
set(gcf, 'Position',  [200, 200, 800, 800])

%% Start simulation
pauseTime = 0.3;
for it = 1:MaxIteration

    pause(pauseTime);
    
    % Update the figure
    UpdatePlot();
    % Update the cells (in parallel)
    UpdateCell();
end
disp("CellularAutomata3D Done.");

%==========================================================================

%% nested function: Activate one cell
function [] = Activate(i, j, k)
    
    % front
    if (i<size_)
        grid_(i+1, j, k) = grid_(i+1, j, k)+4;
    end
    % right
    if (j<size_)
        grid_(i, j+1, k) = grid_(i, j+1, k)+8;
    end
    % back
    if (i>1)
        grid_(i-1, j, k) = grid_(i-1, j, k)+1;
    end
    % left
    if (j>1)
        grid_(i, j-1, k) = grid_(i, j-1, k)+2;
    end
    % up
    if (k<size_)
        grid_(i, j, k+1) = grid_(i, j, k+1)+16;
    end
    % down
    if (k>1)
        grid_(i, j, k-1) = grid_(i, j, k-1)+32;
    end
end


%% nested function: Deactivate one cell
function [] = Deactivate(i, j, k)
    
    % front
    if (i<size_)
        grid_(i+1, j, k) = grid_(i+1, j, k)-4;
    end
    % right
    if (j<size_)
        grid_(i, j+1, k) = grid_(i, j+1, k)-8;
    end
    % back
    if (i>1)
        grid_(i-1, j, k) = grid_(i-1, j, k)-1;
    end
    % left
    if (j>1)
        grid_(i, j-1, k) = grid_(i, j-1, k)-2;
    end
    % up
    if (k<size_)
        grid_(i, j, k+1) = grid_(i, j, k+1)-16;
    end
    % down
    if (k>1)
        grid_(i, j, k-1) = grid_(i, j, k-1)-32;
    end
end


%% nested function: Transform 1D index to {x, y, z} index
function [res] = idx2xyz(list)

    if (size(list, 1)~=1), list = list';end
    t = mod(list, size_*size_);  % used to determine x and y
    t = t - 1;
    t(t < 0) = size_ * size_ - 1;
    % x = mod(t, size_) + 1
    % y = ceil((t+1) / size_)
    res = [mod(t, size_) + 1; ceil((t+1) / size_); ceil(list/(size_*size_))];
end


%% nested function: implemented with parallel computation
function [] = UpdatePlot()

    persistent cmap
    if isempty(cmap)
        cmap = parula(102);
    end

    acList = find(state);  % who is activated now
    if isempty(acList), it=MaxIteration;end
    xyzList = idx2xyz(acList);
    dist = sum((xyzList - size_/2).^2);  % distance, used to determine color
    dist = sqrt(dist) ./ (size_/2 * 1.74);
    
    set(h, 'xdata', xyzList(1, :), 'ydata', xyzList(2, :), 'zdata', xyzList(3, :), 'cdata', cmap(floor(dist / 0.01)+1, :));
end


%% nested function: update cell activity
function [] = UpdateCell()

    newState = rule(grid_ + 1);  % apply rule based on neighbors
    % activate & deactivate
    activateList = find((state == 0) & (newState == 1));
    deactivateList = find((state == 1) & (newState == 0));
    state = newState;
    xyzList = idx2xyz(activateList);
    for ii = 1:size(xyzList, 2)
        Activate(xyzList(1, ii), xyzList(2, ii), xyzList(3, ii));
    end
    xyzList = idx2xyz(deactivateList);
    for ii = 1:size(xyzList, 2)
        Deactivate(xyzList(1, ii), xyzList(2, ii), xyzList(3, ii));
    end
end
end

% End of CellularAutomata3D.m
