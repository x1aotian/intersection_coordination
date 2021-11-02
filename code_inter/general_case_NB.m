function general_case_NB(Ls,Do,Di,Nv,V,Bdy,T,N,NK)

    %% Description
    % Take the center of the intersection as the origin;
    % Ls = side length;
    % Do = safe distance outside;
    % Di = sade distance inside;
    % Nv = number of vehicles;
    % V = [enter index, exit index, start position, start velocity,
    %      priority index] = all information of the vehicles;
    % e.g. V = [1,4,-10,8, 5,0;
    %           3,6,-7,10, 3,0];
    % Bdy = [v_min,v_max;a_min,a_max];
    % T = one horizon predication time;
    % N = number of iteration in one horizon predication time;
    % Na = number of apllied steps
    % NK = number of big iterations;
    
    %% Model
    addpath('D:\Control_System\casadi')
    import casadi.*
    
    % dynamic: x_k+1 = A*x_k + B*u_k
    Ts = T/N;
    A = [1,Ts;0,1];
    B = [0.5*Ts^2;Ts];
    
    Na = 1;
        
    % priority set and inside/outside index set
    ee_set = V(:,1:2);
    p_set = V(:,5); % default 1-Nv
    
    % implement the io_gene function
    io_cors = zeros(Nv,2);
    for i = 1:Nv
        io_cors(i,:) = io_gene(ee_set(i,:),Ls);
    end
    
    % boundary
    v_min = Bdy(1,1);
    v_max = Bdy(1,2);
    v_ref = (v_max+v_min)/2;
    a_min = Bdy(2,1);
    a_max = Bdy(2,2);
    
    % the intial states
    x_ini = V(:,3:4);
    x_last = reshape(x_ini',2*Nv,1);
    
    % reserve some space for trajectories
    xux_all = cell(Nv,2*NK+1);
    x_trj = zeros(2*Nv,NK+1);
    appl_trj = zeros(2*Nv,1);
    pxy_trj = zeros(2*Nv,NK+1);
    
    for i = 1:Nv
        pv_temp = x_ini(i,:)';
        xux_all{i,1} = pv_temp;
        xy_temp = full(path_gene(ee_set(i,:),Ls,pv_temp(1),20));
        pxy_trj(2*i-1:2*i,1) = xy_temp';
    end
    x_trj(:,1) = x_last;
    
    %% MPC loop
    for K = 1:NK
        
        % reserve some space to w and w0
        w = cell(Nv,2*N+1);
        w0 = zeros(3*Nv*N+2*Nv,1);
        lbw = zeros(3*Nv*N+2*Nv,1);
        ubw = zeros(3*Nv*N+2*Nv,1);
        g = {};
        lbg = [];
        ubg = [];
        J = 0;
        
        % initial conditions       
        w(:,1) = MXcell(1,Nv,['X' int2str(K-1) '_'],2);
        w0(1:2*Nv) = x_last;
        lbw(1:2*Nv) = x_last;
        ubw(1:2*Nv) = x_last;
        for k = 1:N
            
            lbw(3*Nv*k-2:3*Nv*k) = a_min*ones(3,1);
            ubw(3*Nv*k-2:3*Nv*k) = a_max*ones(3,1);
            lbw(3*Nv*k+1:2:3*Nv*k+5) = -inf;
            ubw(3*Nv*k+1:2:3*Nv*k+5) = +inf;
            lbw(3*Nv*k+2:2:3*Nv*k+6) = v_min;
            ubw(3*Nv*k+2:2:3*Nv*k+6) = v_max;
            
        end
            
        % collision avoidance priciples             
        
        I_set = [];  % set of inside vehicles' index
        R_set = cell(8,1);  % set of outside vehicles' index, cell size = 8*1
        I_col_set = [];
        
        for i = 1:Nv
            in_index = ee_set(i,1);
            out_index = ee_set(i,2);
            if x_last(2*i-1) < io_cors(i,1)
                R_set{in_index,1} = [R_set{in_index,1},i];
            elseif x_last(2*i-1) > io_cors(i,2)
                R_set{out_index,1} = [R_set{out_index,1},i];
            else
                I_set = [I_set,i];
            end
        end
        
        if size(I_set) > 1
            I_col_set = nchoosek(I_set,2);
        end
        R_col_set = cell(8,1);
        margin_set = [];
        for i = 1:8
            num_v = length(R_set{i,1});
            for j = 1 : num_v-1
                R_col_set{i,1} = [R_col_set{i,1};R_set{i,1}(1,j:j+1)];
            end
            if size(R_set{i,1},1) > 0
                if ismember(i,[1,3,5,7])
                    margin_set = [margin_set,R_set{i,1}(1)];
                else
                    margin_set = [margin_set,R_set{i,1}(num_v)];
                end
            end
        end
        [mm,nn] = meshgrid(I_set',margin_set);
        IR_col_set = [reshape(mm,[],1), reshape(nn,[],1)];
        
        cnt = 0;
        for k = 0+K-1:N+K-2
            cnt = cnt + 1;
            w(:,2*cnt) = MXcell(1,Nv,['U' int2str(k) '_'],1);       
            
            % next X
            w(:,2*cnt+1) = MXcell(1,Nv,['X' int2str(k+1) '_'],2);
            
            % use dynamic to make sure state continous
            for h = 1:Nv
                g = [g {w{h,2*cnt+1} - (A * w{h,2*cnt-1} + B * w{h,2*cnt})}];
                lbg = [lbg; 0; 0];
                ubg = [ubg; 0; 0];
                % reference
                J = J + w{h,2*cnt}^2 + (w{h,2*cnt+1}(2)-v_ref)^2;
            end
            
            % inside collision avoidance
            I_len = size(I_col_set,1);
            for i = 1:I_len
                a_index = I_col_set(i,1);
                b_index = I_col_set(i,2);
                a_cor = path_gene(ee_set(a_index,:),Ls,w{a_index,2*cnt+1}(1),20);
                b_cor = path_gene(ee_set(b_index,:),Ls,w{b_index,2*cnt+1}(1),20);
                g = [g {norm(a_cor-b_cor)}];
                lbg = [lbg; Di];
                ubg = [ubg; +inf];
            end
            
            % outside collision avoidance
            R_col_mat = cell2mat(R_col_set);
            R_len = size(R_col_mat,1);
            for i = 1:R_len
                a_index = R_col_mat(i,1);
                b_index = R_col_mat(i,2);
%                 a_cor = path_gene(ee_set(a_index,:),Ls,w{a_index,2*cnt+1}(1),20);
%                 b_cor = path_gene(ee_set(b_index,:),Ls,w{b_index,2*cnt+1}(1),20);
%                 g = [g {norm(a_cor-b_cor)}];
                g = [g { w{a_index,2*cnt+1}(1)-w{b_index,2*cnt-1}(1) }];
                lbg = [lbg; Do];
                ubg = [ubg; +inf];                   
            end
            
            % in and out collision avoidance
            IR_len = size(IR_col_set,1);
            for i = 1:IR_len
                a_index = IR_col_set(i,1);
                b_index = IR_col_set(i,2);
                a_cor = path_gene(ee_set(a_index,:),Ls,w{a_index,2*cnt+1}(1),20);
                b_cor = path_gene(ee_set(b_index,:),Ls,w{b_index,2*cnt+1}(1),20);
                g = [g {norm(a_cor-b_cor)}];
                lbg = [lbg; Di];
                ubg = [ubg; +inf];
            end       
        end
        
        % create an NLP solver
        prob = struct('f',J,'x',vertcat(w{:}),'g',vertcat(g{:}));
        setting.ipopt.print_level = 5;
        solver = nlpsol('solver','ipopt',prob,setting);
        
        % solve the NLP
        sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
        w_opt = full(sol.x);
        

        for i = 1:Nv
            xux_all{i,2} = w_opt(3*Nv-3+i); %u0
            xux_all{i,21} = w_opt(3*Nv+2*i-1: 3*Nv+2*i); %x1
            appl_trj(2*i-1:2*i,1) = w_opt(3*Nv+2*i-1: 3*Nv+2*i);
        end
        
        x_last = appl_trj;
  
        for i = 1:Nv
            x_trj(:,K+1) = x_last;
            xy_temp = full(path_gene(ee_set(i,:),Ls,x_last(2*i-1),20));
            pxy_trj(2*i-1:2*i,K+1) = xy_temp';
            
        end
        
    end
end
    
