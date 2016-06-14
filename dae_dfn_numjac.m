function ydot = dae_dfn_numjac(t,y,cur,p)

x = y(1:size(p.f_x,1));
z = y(size(p.f_x,1)+1 : end);

[f,g] = dae_dfn_federico_scott(x,z,cur,p);
ydot = [f; g];