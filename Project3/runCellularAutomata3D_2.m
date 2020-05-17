% runCellularAutomata3D_2

clear o
o.size = 40;  % N*N*N grid
o.init = [19, 21, 20, 20, 20, 20;...
          20, 20, 19, 21, 20, 20;...
          20, 20, 20, 20, 19, 21];  % initially activated cells
rng(4)
o.rule = randi(2, 1, 64)-1;  % 64 bits 0/1
o.it = 100;  % number of iterations
CellularAutomata3D(o);
