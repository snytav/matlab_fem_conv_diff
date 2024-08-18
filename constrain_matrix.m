function [A_ff,A_fp,A_pf,A_pp]=constrain_matrix(A,dof_constrained)

% Constrain a matrix
N=length(A);
p=dof_constrained;
f_aus=1:N;
p_aus=zeros(1,N);
p_aus(p)=p;
f=f_aus-p_aus;
f=find(f);

A_ff=A(f,f);
dlmwrite('c:\ocean\conv_diff\A_ff_cons_m.txt',A_ff,'delimiter','\n','precision',15);
A_fp=A(f,p);
dlmwrite('c:\ocean\conv_diff\A_fp_cons_m.txt',A_fp,'delimiter','\n','precision',15);

A_pf=A(p,f);
dlmwrite('c:\ocean\conv_diff\A_pf_cons_m.txt',A_pf,'delimiter','\n','precision',15);
A_pp=A(p,p);
dlmwrite('c:\ocean\conv_diff\A_pp_cons_m.txt',A_pp,'delimiter','\n','precision',15);

end