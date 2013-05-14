%Plot Epx to the nth order in the same plot

%Load epx data file
epx_data=load('Epxn_mckln.dat');

order_begin = 2;
order_end = size(epx_data,2)-1;
