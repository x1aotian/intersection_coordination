function [pltx, plty] = ellip_gen(x,y,rad,ls,ss)
theta = 0:pi/12:2*pi;
if rad ~= 0
    RT = [ cos(rad) -sin(rad);
           sin(rad)  cos(rad)];
    pltx = ls*cos(theta);
    plty = ss*sin(theta);
    pltxy = RT * [pltx; plty];
    pltx = pltxy(1,:);
    plty = pltxy(2,:);
else
    pltx = ls*cos(theta);
    plty = ss*sin(theta);
end

if x~=0
    pltx = pltx + x;
end
if y~=0
    plty = plty + y;
end

end