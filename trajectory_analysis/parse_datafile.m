function [x,y] = parse_datafile(filename)
    x = readmatrix(filename,'Range','Z:Z');
    x = x(4:end);
    y = readmatrix(filename,'Range','AA:AA');
    y = y(4:end);
end