% runCellularAutomata3D_4

clear o
o.size = 40;  % N*N*N grid
o.init = [19, 21, 20, 20, 20, 20;...
          20, 20, 19, 21, 20, 20;...
          20, 20, 20, 20, 19, 21];  % initially activated cells
rng(3)
o.rule = randi(5, 1, 64) <= 2;  % 64 bits 0/1
o.it = 100;  % number of iterations
CellularAutomata3D(o);
