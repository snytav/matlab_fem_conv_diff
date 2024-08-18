function [y] = u_0_fun(x)
    u_max=5/7;                              % Peak of Gauss hille
    x_0=2/15;
    l=7*sqrt(2)/300;                        % Width of Gauss hill
    y = u_max*exp(-((x-x_0)/l).^2);
end