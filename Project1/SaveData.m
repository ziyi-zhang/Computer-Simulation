% Save data to 'mat' file after running spring simulation
% The data will then be used by 'GravityField' to create movie of gravity
% field

data.M = M;
data.L = L;
data.V = V;
data.A = A;
data.Spring = Spring;
data.colorSprIdx = colorSprIdx;

data.tmax = tmax;
data.clockmax = clockmax;

if storeData
    data.XHist = XHist;
    data.YHist = YHist;
end

save([filename, '.mat'], 'data');
