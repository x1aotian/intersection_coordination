function output = io_gene(ee_i,Ls,Lb)
    R_set = [1,8;3,2;5,4;7,6];
    L_set = [1,4;3,6;5,8;7,2];
    if ismember(ee_i,R_set,'rows')
        d = 0.25*pi*0.25*Ls;
        p_in = -d;
        p_out = d;
    elseif ismember(ee_i,L_set,'rows')
        d = 0.75*pi*0.25*Ls;
        p_in = -d;
        p_out = d;
    else
        p_in = -0.5*Ls;
        p_out = 0.5*Ls;
    end
    output = [p_in-Lb,p_out+Lb];
end