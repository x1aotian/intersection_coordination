addpath('D:\Control_System\casadi')
import casadi.*

%% parameters

% structure of roads
Nr = [1;2;3;4;5;6;7;8];
Er_pi = [1, 1,2; 2, 2,3; 3, 2,4; 4, 4,5; 5, 4,6; 6, 5,7; 7, 5,8];
Er = Er_pi(:,2:3);
num_r = size(Er,1);
Nm = [2; 4; 5]; % merge points
Nm_p = [2; -10; -14]; % position of merge points
N_rs = [1, 2,3; 3, 4,5; 4, 6,7]; % crossed roads at the merge points
num_m = size(Nm,1);
Theta_m = [1/6*pi; 1/3*pi; 1/3*pi]; % Angle at merge points

%initial state of vehicles
% V0 = [1, 2, 4, -2, 8; % [index; N_ahead, N_behind, p_i, v_i]
%       2, 2, 4, -7.5, 8;
%       3, 2, 3, -5, 8;
%       4, 4, 5, -11, 8;
%       5, 4, 6, -14, 8;
%       6, 5, 7, -18, 8;
%       7, 5, 8, -21, 8]; 
% a_min = -20; a_max = 20; v_min = -0; v_max = 15; v_ref = 8; ds = 3;
V0 = [1, 2, 4, -2, 5; 
      2, 2, 4, -7, 8;
      3, 2, 3, -5, 8;
      4, 4, 5, -11, 4;
      5, 4, 6, -13, 7;
      6, 5, 7, -15, 3;
      7, 5, 8, -17, 9]; 
a_min = -100; a_max = 100; v_min = -0; v_max = 30; v_ref = 8; ds = 3;



num_V = size(V0,1);
V_ini = V0;
x_ini = reshape(V_ini(:,4:5)',2*num_V,1);

% vehicle paras


% MPC paras
T = 0.5; Nk = 10; ts = T/Nk;  
NK = 20;
A = [1,ts;0,1]; B = [0.5*ts^2;ts];

% reserve space for trajectories
xux_all = zeros(num_V*(3*NK+2),1); 
xux_all(1:2*num_V) = x_ini;
% pxy_trj = zeros(2*num_V,NK+1);

p_trj = zeros(num_V,NK+1);
p_trj(:,1) = V0(:,4);
v_trj = zeros(num_V,NK+1);
v_trj(:,1) = V0(:,5);
a_trj = zeros(num_V,NK);

%% MPC scheme
for KK = 1:NK
    
    %% construct tree topology
    R = cell(num_r,1); % R{i}: all vehilces in road i
    for i = 1:num_V
        r_i = V_ini(i,2:3);
        temp_bool = ismember(Er,r_i,'rows');
        ind_i = find(temp_bool==1);
        R{ind_i} = [R{ind_i}; V_ini(i,:)];
    end

    R_ab = zeros(num_r,2); % R_ab(i) = [front est, behind est]
    for i = 1:num_r  
        if size(R{i},1) ~= 0
            R{i} = sortrows(R{i},4);   
            R_ab(i,1) = R{i}(end,1);
            R_ab(i,2) = R{i}(1,1);
        end
    end

    N = V0(:,1);
    E_p = [];

    for i = 1:num_r  % edges in the same road
        num_vinr =  size(R{i},1);
        if num_vinr > 1
            for j = 1:num_vinr-1
                E_p = [E_p; R{i}(j+1,1),R{i}(j,1),0,0];
            end
        end
    end

    for i = 1:num_m  % edges in the merge point
        ris =  N_rs(i,:);
        aa = R_ab(ris(1),2); bb = R_ab(ris(2),1); cc = R_ab(ris(3),1);
        temp_m = [aa; bb; cc];
        temp_m(find(temp_m==0)) = [];
        size_m = size(temp_m,1);
        temp_M = [];
        if size_m >= 2
            for j = 1:size_m
                temp_M = [temp_M; temp_m(j), V_ini(temp_m(j),4)];
            end
            temp_M = sortrows(temp_M,2);
            for k = 1:size_m-1
                if (temp_M(k+1,1)==bb && temp_M(k,1)==cc) || (temp_M(k,1)==bb && temp_M(k+1,1)==cc)
                    E_p = [E_p; temp_M(k+1,1), temp_M(k,1),Theta_m(i),Nm_p(i)];
                else
                    E_p = [E_p; temp_M(k+1,1), temp_M(k,1),0,0];
                end
            end
        end
    end
    E = E_p(:,1:2);
    
    %% construct variables
    % reserve space for w and w0
    w = cell(num_V,2*Nk+1);
    lbw = zeros(num_V*(3*Nk+2),1);
    ubw = zeros(num_V*(3*Nk+2),1);
    w0 = zeros(num_V*(3*Nk+2),1);
    g = {}; lbg = []; ubg = []; J = 0;
    
    % initial condition
    w(:,1) = MXcell(1,num_V,['X' int2str(0) '_'],2);
    lbw(1:2*num_V) = x_ini;
    ubw(1:2*num_V) = x_ini;
    
    % construct lbw and ubw
    for k = 1:Nk
        lbw(3*num_V*k-num_V+1:3*num_V*k) = a_min*ones(num_V,1);
        ubw(3*num_V*k-num_V+1:3*num_V*k) = a_max*ones(num_V,1);
        lbw(3*num_V*k+1:2:3*num_V*k+2*num_V-1) = -inf;
        ubw(3*num_V*k+1:2:3*num_V*k+2*num_V-1) = +inf;
        lbw(3*num_V*k+2:2:3*num_V*k+2*num_V) = v_min;
        ubw(3*num_V*k+2:2:3*num_V*k+2*num_V) = v_max;
    end
    
    %% predict horizon
    for kk = 1:Nk
        w(:,2*kk) = MXcell(1,num_V,['U' int2str(kk-1) '_'],1);
        w(:,2*kk+1) = MXcell(1,num_V,['X' int2str(kk) '_'], 2);
        
        for h = 1:num_V
            g = [g {w{h,2*kk+1} - (A * w{h,2*kk-1} + B * w{h,2*kk})}];
            lbg = [lbg; 0; 0];
            ubg = [ubg; 0; 0];
            J = J + 0.2*w{h,2*kk}^2 + (w{h,2*kk+1}(2)-v_ref)^2;
        end
       
        num_E = size(E,1);
        
        if kk >= 5
            for e = 1:num_E
                ai = E(e,1); bi = E(e,2); theta_ab = E_p(e,3); p_cross = E_p(e,4);
                ap = w{ai,2*kk+1}(1); bp = w{bi,2*kk+1}(1);
                if theta_ab == 0
                    g = [g {ap-bp}];
                else
                    side_a = p_cross - ap; side_b = p_cross - bp;
                    %d_temp = (cos(theta_ab)*side_b-side_a)^2 + (sin(theta_ab)*side_b)^2;
                    d_temp = (side_b-side_a)*cos(0.5*theta_ab);
                    g = [g {d_temp}];
                end
                lbg = [lbg;ds];
                ubg = [ubg;+inf];
            end    
        end
    end 
    
    % create an NLP solver
    prob = struct('f',J,'x',vertcat(w{:}),'g',vertcat(g{:}));
    setting.ipopt.print_level = 5;
    solver = nlpsol('solver', 'ipopt', prob, setting);

    % solve the NLP
    sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
    w_opt = full(sol.x);
    
    x_ini = w_opt(3*num_V+1:5*num_V);
    
    % use x_ini to update V_ini
    for i = 1:num_V
        V_ini(i,4:5) = x_ini(2*i-1:2*i);
        if ismember(V_ini(i,2),Nm)
            m_ind = find(Nm==V_ini(i,2));
            m_p = Nm_p(m_ind);
            if V_ini(i,4) > m_p
                new_ind = find(Er(:,2)==V_ini(i,2));
                V_ini(i,2:3) = Er(new_ind,:);
            end
        end
    end
    
    xux_all(num_V*(5*KK-2)+1: num_V*5*KK) = x_ini;
    xux_all(num_V*(3*KK-1)+1: num_V*3*KK) = w_opt(2*num_V+1:3*num_V);
    a_trj(:,KK) = w_opt(2*num_V+1:3*num_V);
    p_trj(:,KK+1) = x_ini(1:2:end);
    v_trj(:,KK+1) = x_ini(2:2:end);
end

%% plot
t_v = [0:ts:ts*NK];
t_a = [0:ts:ts*(NK-1)];
style_set = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y';'-.k';];

figure(1);
for i = 1:num_V
    plot(t_v,p_trj(i,:),style_set(i));
    hold on;
end
title('trajectories of scalar position');
ylabel('s-m');
xlabel('t-s');
grid on;


figure(2);
for i = 1:num_V
    plot(t_v,v_trj(i,:),style_set(i));
    hold on;
end
title('velocity');
ylabel('v-m/s');
xlabel('t-s');

figure(3);
for i = 1:num_V
    plot(t_a,a_trj(i,:),style_set(i));
    hold on;
end
title('acceleration');
ylabel('a-m/s^2');
xlabel('t-s');


%% 