%Script solving the stochastic particle moving in a potential field with
%friction
function fevals = Potential_friction_motion(eval_pts)
addpath('/home/hemaditya/Documents/chebfun'); %Adding chebfun if not added


f = 2; %Larger this value, faster the steady state is reached
N = chebop(@(x, u) diff(u,2) + f*diff(u) + 35*0.5*u.^3 - 15*0.5*u, [0, 10]); %ODE
ic = 0.05 + 0.2*(eval_pts);%Initial condition

fevals = zeros(length(eval_pts),1);
for iter=1:length(ic)
N.lbc = [ic(iter); 0];
u = N\0;
fevals(iter,1) = u(end);
end


end