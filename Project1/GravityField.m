function [] = GravityField(filename)
% Visualize the gravity field based on the stored data
% 'filename' should be a 'mat' file created by 'SaveData.m'

    % read in data
    data = load(filename);
    data = data.data;
    M = data.M;
    L = data.L;
    V = data.V;
    A = data.A;
    Spring = data.Spring;
    tmax = data.tmax;
    clockmax = data.clockmax;
    dt = tmax / clockmax;
    if ~isfield(data, 'XHist') || ~isfield(data, 'YHist')
        warning('Please run simulation before running this script.');
        return;
    end
    XHist = data.XHist;
    YHist = data.YHist;
    
    % 
end
