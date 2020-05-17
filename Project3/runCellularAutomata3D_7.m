% runCellularAutomata3D_7

clear o
o.size = 40;  % N*N*N grid
o.init = randi(40, 3, 100);  % initially activated cells
rng(13)
o.rule = randi(7, 1, 64) <= 2;  % 64 bits 0/1
o.it = 100;  % number of iterations
CellularAutomata3D(o);
