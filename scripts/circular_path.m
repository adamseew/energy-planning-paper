function Path = circular_path(p, xo, yo, r)

xel = p(1) - xo;
yel = p(2) - yo;

phi = xel*xel + yel*yel - r*r;

Path.e = phi;
Path.grad = [2*xel; 2*yel];
Path.Hessian = 2*eye(2);

end