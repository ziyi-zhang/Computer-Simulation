% runCellularAutomata3D_5

clear o
o.size = 40;  % N*N*N grid
o.init = repmat(2:39, 3, 1);  % initially activated cells
rng(7)
o.rule = randi(7, 1, 64) <= 2;  % 64 bits 0/1
o.it = 100;  % number of iterations
CellularAutomata3D(o);
