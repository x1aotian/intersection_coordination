function output1 = MXcell(NN0,NN1,x,nn)
    addpath('D:\Control_System\casadi')
    import casadi.*
    output1 = cell(NN1-NN0+1,1);
    for i = NN0:NN1
        output1{i} = MX.sym([x int2str(i)],nn);
    end
end  