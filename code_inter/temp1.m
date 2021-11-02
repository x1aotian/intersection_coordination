addpath('D:\Control_System\casadi')
import casadi.*

P = [-2.0,-1.0,0.0,1.0,2.0];
X = [-2,-1,0,1,1];
Y = [1,1,2,3,4];

px = casadi.interpolant('PX','bspline',{P},X);
py = casadi.interpolant('PY','bspline',{P},Y);

x = MX.sym('x');
y = MX.sym('y');
f = x^2+y^2;
g = px(x)+py(y);
prob = struct('f',f,'g',g,'x',[x;y]);
F = nlpsol('F','ipopt',prob);
r = F('x0',[1;1],'ubg',5,'lbg',-inf);