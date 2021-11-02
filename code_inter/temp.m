B = 5;
N = 20;
S = 20;
R = 15;
L = 2*S + 0.5*pi*R;
l = L/(2*N);
P = [];
X = [];
Y = [];
for i = -N:N
    p = i*l;
    if p < -0.25*pi*R
        x = -2*B+(p+0.25*pi*R);
        y = -B;
    elseif p > 0.25*pi*R
        x = B;
        y = 2*B+(p-0.25*pi*R);
    else
        theta = (p+0.25*pi*R)/R;
        x = -2*B+3*B*sin(theta);
        y = 2*B-3*B*cos(theta);
    end
    
    P=[P,p];
    X=[X,x];
    Y=[Y,y];
    
end