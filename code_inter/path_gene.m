function output = path_gene(ee_i,Ls,smp)
    simu_L = 9;
    addpath('D:\Control_System\casadi')
    import casadi.*
    point_cors = [-0.25*Ls,0.5*Ls; 0.25*Ls,0.5*Ls; 0.5*Ls,0.25*Ls; 0.5*Ls,-0.25*Ls;...
                  0.25*Ls,-0.5*Ls; -0.25*Ls,-0.5*Ls; -0.5*Ls,-0.25*Ls; -0.5*Ls,0.25*Ls];
    R_cens = [-0.5*Ls,0.5*Ls;0.5*Ls,0.5*Ls;0.5*Ls,-0.5*Ls;-0.5*Ls,-0.5*Ls];
    R_rads = [0;-0.5*pi;-pi;-1.5*pi];
    R_set = [1,8;3,2;5,4;7,6];
    L_cens = [0.5*Ls,0.5*Ls;0.5*Ls,-0.5*Ls;-0.5*Ls,-0.5*Ls;-0.5*Ls,0.5*Ls];
    L_rads = [pi;0.5*pi;0;-0.5*pi];
    L_set = [1,4;3,6;5,8;7,2];
    unitd_set = [0,-1;0,1;-1,0;1,0;0,1;0,-1;1,0;-1,0];
    pipo = io_gene(ee_i,Ls,0);
    p_in = pipo(1);
    p_out = pipo(2);
    
    thes = [1.5*pi;pi;0.5*pi;0];
    
    if ismember(ee_i,R_set,'row')
        P = [-simu_L*Ls+p_in:simu_L*Ls/smp:p_in-simu_L*Ls/smp, p_in:(p_out-p_in)/(smp-1):p_out, p_out+simu_L*Ls/smp:simu_L*Ls/smp:p_out+simu_L*Ls];
        C = zeros(3,3*smp);
        for i = 1:smp
            C(1:2,i) = point_cors(ee_i(1),:) + (P(i)-p_in)*unitd_set(ee_i(1),:);
            C(3,i) = thes((ee_i(1)+1)/2);
        end
        for i = smp+1:2*smp
            rad = (P(i) - p_in)/(0.25*Ls);
            C(1,i) = R_cens((ee_i(1)+1)/2,1) + 0.25*Ls*cos(-rad+R_rads((ee_i(1)+1)/2));
            C(2,i) = R_cens((ee_i(1)+1)/2,2) + 0.25*Ls*sin(-rad+R_rads((ee_i(1)+1)/2));
            C(3,i) = thes((ee_i(1)+1)/2)-rad;
        end
        for i = 2*smp+1:3*smp
            C(1:2,i) = point_cors(ee_i(2),:) + (P(i)-p_out)*unitd_set(ee_i(2),:);
            C(3,i) = thes((ee_i(1)+1)/2)-0.5*pi;
        end
        output = [P;C];
        
    elseif ismember(ee_i,L_set,'row')
        P = [-simu_L*Ls+p_in:simu_L*Ls/smp:p_in-simu_L*Ls/smp, p_in:(p_out-p_in)/(smp-1):p_out, p_out+simu_L*Ls/smp:simu_L*Ls/smp:p_out+simu_L*Ls];
        C = zeros(3,3*smp);
        for i = 1:smp
            C(1:2,i) = point_cors(ee_i(1),:) + (P(i)-p_in)*unitd_set(ee_i(1),:);
            C(3,i) = thes((ee_i(1)+1)/2);
        end
        for i = smp+1:2*smp
            rad = (P(i) - p_in)/(0.75*Ls);
            C(1,i) = L_cens((ee_i(1)+1)/2,1) + 0.75*Ls*cos(rad+L_rads((ee_i(1)+1)/2));
            C(2,i) = L_cens((ee_i(1)+1)/2,2) + 0.75*Ls*sin(rad+L_rads((ee_i(1)+1)/2));
            C(3,i) = thes((ee_i(1)+1)/2)+rad;
        end
        for i = 2*smp+1:3*smp
            C(1:2,i) = point_cors(ee_i(2),:) + (P(i)-p_out)*unitd_set(ee_i(2),:);
            C(3,i) = thes((ee_i(1)+1)/2)+0.5*pi;
        end       
        output = [P;C];    
        
    else
        direct = unitd_set(ee_i(1),:);
        P = [-simu_L*Ls+p_in:simu_L*Ls/smp:p_in-simu_L*Ls/smp, p_in:(p_out-p_in)/(smp-1):p_out, p_out+simu_L*Ls/smp:simu_L*Ls/smp:p_out+simu_L*Ls];
        C = zeros(3,3*smp);
        for i = 1:3*smp
            C(1:2,i) = (point_cors(ee_i(2),:)+point_cors(ee_i(1),:))/2 + direct*P(i);
            C(3,i) = thes((ee_i(1)+1)/2);
        end
        output = [P;C];
    end

end