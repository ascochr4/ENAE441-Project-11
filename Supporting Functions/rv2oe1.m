function oe = rv2oe1(r,v,mu)
% Takes input position and velocity and returns Keplerian Orbital Elements
a = -mu/2*(norm(v)^2/2-mu/norm(r))^-1;
h = cross(r,v);
e = cross(v,h)/mu-r/norm(r);
ecc = norm(e);
i = acos(h(3)/norm(h))*180/pi;
n = cross([0 0 1],h);
if n(2) > 0
    Om = acos(n(1)/norm(n))*180/pi;
else 
    Om = (2*pi -acos(n(1)/norm(n)))*180/pi;
end
if e(3)>0 
    om = acos(dot(n,e)/(norm(n)*ecc))*180/pi;
else
    om = (2*pi-acos(dot(n,e)/(norm(n)*ecc)))*180/pi;
end
if dot(r,v) > 0
    nu = acos(dot(e,r)/(norm(r)*ecc))*180/pi;
else
    nu = (2*pi -acos(dot(e,r)/(norm(r)*ecc)))*180/pi;
end
oe = [a, ecc, i, Om, om, nu];
end