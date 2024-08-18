function [C,D] = readCD(C_fname,D_fname)

C=dlmread(C_fname);
D=dlmread(D_fname);
end