function general_case(Ls,Do,Di,Nv,V,Bdy,T,N,Na,NK)

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
        
    % priority set and inside/outside index set
    ee_set = V(:,1:2);
    p_set = V(:,5); % default 1-Nv
    
    % implement the io_gene function
    io_cors = zeros(Nv,2);
    for i = 1:Nv
        io_cors(i,:) = io_gene(ee_set(i,:),Ls);
    end
    
    % boundary
    v_min = Bdy(1);
    v_max = Bdy(2);
    v_ref = (v_max+v_min)/2;
    a_min = Bdy(3);
    a_max = Bdy(4);
    
    % the intial states
    x_ini = V(:,3:4);
    x_last = x_ini;
    
    % reserve some space for trajectories
    trj = zeros(NK,3*Nv);
    trj_pred = zeros(N,3*Nv);
    trj_appl = zeros(Na,3*Nv);
    
    %% MPC loop
    for K = 1:NK
        
        % reserve some space to w and w0
        w = cell(Nv,2*N+1);
        w0 = zeros(Nv,3*N+2);
        lbw = zeros(Nv,3*N+2);
        ubw = zeros(Nv,3*N+2);
        g = cell(Nv,N);
        lbg = zeros(Nv,2*N);
        ubg = zeros(Nv,2*N);
        J = 0;
        
        % initial conditions       
        w(:,1) = MXcell(1,Nv,['X' int2str(K-1) '_'],2);
        w0(:,1:2) = x_last;
        lbw(:,1:2) = x_last;
        ubw(:,1:2) = x_last;
        lbw(:,3:end) = -inf;
        ubw(:,3:end) = +inf;
        
        % collision avoidance priciples             
        
        I_set = [];  % set of inside vehicles' index
        R_set = cell(8,1);  % set of outside vehicles' index, cell size = 8*1
        I_col_set = [];
        
        for i = 1:Nv
            in_index = ee_set(i,1);
            out_index = ee_set(i,2);
            if x_last(i,1) < io_cors(i,1)
                R_set{in_index,1} = [R_set{in_index,1},i];
            elseif x_last(i,1) > io_cors(i,2)
                R_set{out_index,1} = [R_set{out_index,2},2];
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
        [IR_col_set(:,1),IR_col_set(:,2) ] = deal( reshape(mm,[],1), reshape(nn,[],1) );
        
        cnt = 0;
        for k = 0+K-1:N+K-2
            cnt = cnt + 1;
            w(:,2*cnt) = MXcell(1,Nv,['U' int2str(k) '_'],1);       
            
            % next X
            w(:,2*cnt+1) = MXcell(1,Nv,['X' int2str(k+1) '_'],2)
            
            % use dynamic to make sure state continous
            for h = 1:Nv
                g{h,cnt} = w{h,2*cnt+1} - (A * w{h,2*cnt-1} + B * w{h,2*cnt});
                % boundary
                J = J - log(w{h,2*cnt} - a_min) - log(-w{h,2*cnt} + a_max); % a
                J = J - log(w{h,2*cnt+1}(2) - v_min) - log(-w{h,2*cnt+1}(2) + v_max); % v
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
                J = J - log(norm(a_cor-b_cor) - Di);
            end
            
            % outside collision avoidance
            R_col_mat = cell2mat(R_col_set)
            R_len = size(R_col_mat,1);
            for i = 1:R_len
                if size(R_col_set{i,1},2) ~= 0
                    a_index = R_col_set{i,1}(1);
                    b_index = R_col_set{i,1}(2);
                    a_cor = path_gene(ee_set(a_index,:),Ls,w{a_index,2*cnt+1}(1),20);
                    b_cor = path_gene(ee_set(b_index,:),Ls,w{b_index,2*cnt+1}(1),20);
                    J = J - log(norm(a_cor-b_cor) - Do);
                end
            end
            
            % in and out collision avoidance
            IR_len = size(IR_col_set,1);
            for i = 1:IR_len
                a_index = IR_col_set(i,1);
                b_index = IR_col_set(i,2);
                a_cor = path_gene(ee_set(a_index,:),Ls,w{a_index,2*cnt+1}(1),20);
                b_cor = path_gene(ee_set(b_index,:),Ls,w{b_index,2*cnt+1}(1),20);
                J = J - log(norm(a_cor-b_cor) - Di);
            end       
        end
        
        % create an NLP solver
        prob = struct('f',J,'x',vertcat(w{:}),'g',vertcat(g{:}));
        setting.ipopt.print_level = 5;
        solver = nlpsol('solver','ipopt',prob,setting);
        
        % solve the NLP
        sol = solver('x0',vertcat(w0(:)),'lbx',vertcat(lbw(:)),'ubx',vertcat(ubw(:)),'lbg',vertcat(lbg(:)),'ubg',vertcat(ubg(:)));
        w_opt = full(sol.x);
        
        result = 1;
        trj_pred = 1;
        trj_appl = 1;
        
        x_new = 1;
        x_last = 1;
        
        trj = 1;
        
    end
    
end
