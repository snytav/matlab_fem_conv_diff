function [num_years,depths_uniform,distances_uniform] = read_Ameland()
distances_uniform = dlmread('dist_unif.txt');
num_years = dlmread('num_years.txt');

depths_uniform = dlmread('dep.txt');
depths_uniform = reshape(depths_uniform,size(distances_uniform,1),size(num_years,1));
[X,Y] = ndgrid(distances_uniform,num_years);
surf(X,Y,depths_uniform);
xlabel('distance');
ylabel('Year');
end