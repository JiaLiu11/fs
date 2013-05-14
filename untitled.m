clear all;
close all;
cubic_ori = load('blockTable_cubic.dat');
cubic_set0 = load('blockTable_cubic_set0.dat');
linear = load('blockTable_linear.dat');
sum(sum(cubic_ori));
sum(sum(cubic_set0));
sum(sum(linear));
rel_err = abs(linear - cubic_ori)./(abs(linear) + abs(cubic_ori) + 1e-16);
ave_rel_err = sum(sum(rel_err))/261/261

rel_err_set0_linear = abs(linear - cubic_set0)./(abs(linear) ...
    + abs(cubic_set0) + 1e-16);
ave_rel_err_set0_linear = sum(sum(rel_err_set0_linear))/261/261

rel_err_cubics = abs(cubic_set0 - cubic_ori)./(abs(cubic_set0) + abs(cubic_ori) + 1e-16);
ave_rel_err = sum(sum(rel_err_cubics ))/261/261

x= -13:0.1:13;
y=-13:0.1:13;
[X,Y]=meshgrid(x,y);
surf(X,Y,rel_err')

figure
surf(X,Y,rel_err_set0_linear')

