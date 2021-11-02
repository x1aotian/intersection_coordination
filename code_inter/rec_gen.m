function [x_list,y_list] = rec_gen(x,y,rad,ls,ss)
    v_ss = [-ss*sin(rad),ss*cos(rad)];
    v_ls = [-ls*cos(rad),-ls*sin(rad)];
    points = [[x,y] - 0.5*v_ss - 0.5*v_ls;...
              [x,y] + 0.5*v_ss - 0.5*v_ls;...
              [x,y] + 0.5*v_ss + 0.5*v_ls;...
              [x,y] - 0.5*v_ss + 0.5*v_ls];
    x_list = points(:,1);
    y_list = points(:,2);

end