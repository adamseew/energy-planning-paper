function [u_theta] = gvf_control_2D(p, dot_p, ke, kd, Path, direction)

e = Path.e;
H = Path.Hessian;
n = Path.grad;
E = [0 1; -1 0];
tau = direction*E*n;

dot_pd = tau - ke*e*n; % (7)
ddot_pd = (E - ke*e)*H*dot_p - ke*n'*dot_p*n; % (10)
ddot_pdhat = -E*(dot_pd*dot_pd')*E*ddot_pd / norm(dot_pd)^3; % (9)

dot_Xid = ddot_pdhat'*E*dot_pd/norm(dot_pd); % (13)

u_theta = dot_Xid + kd*dot_p'*E*dot_pd/(norm(dot_p)*norm(dot_pd)); % (16)

end

