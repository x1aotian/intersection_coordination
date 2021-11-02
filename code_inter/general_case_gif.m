function output = general_case_gif(Ls,Lb,Do,Di,Nv,V,Bdy,T,N,NK)

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
    % diff = 0.4452;

    
    %% Model
    addpath('D:\Control_System\casadi')
    import casadi.*
    
    % dynamic: x_k+1 = A*x_k + B*u_k
    Ts = T/N;
    A = [1,Ts;0,1];
    B = [0.5*Ts^2;Ts];
    
    % car's shape
    car_ls = 2.3;
    car_ss = 1.5;
    L = 50;
    v_ref = 10;
        
    % priority set and inside/outside index set
    ee_set = V(:,1:2);
    
    % implement the io_gene function
    io_cors = zeros(Nv,2); %(x,y) of in and out points
    for i = 1:Nv
        io_cors(i,:) = io_gene(ee_set(i,:),Ls,Lb);
    end
    
    % boundary
    v_min = Bdy(1,1);
    v_max = Bdy(1,2);
    a_min = Bdy(2,1);
    a_max = Bdy(2,2);
    
    %% Cell of coordination generate functions
    F_cell = cell(Nv,3);
    for i = 1:Nv
        temp_data = path_gene(ee_set(i,:),Ls,20);
        px = casadi.interpolant('PX','bspline',{temp_data(1,:)},temp_data(2,:));
        py = casadi.interpolant('PY','bspline',{temp_data(1,:)},temp_data(3,:));
        pthe = casadi.interpolant('PThe','bspline',{temp_data(1,:)},temp_data(4,:));
        F_cell{i,1} = px;
        F_cell{i,2} = py;
        F_cell{i,3} = pthe;
    end
    
    %% Initial
    % the intial states
    x_ini = V(:,3:4);
    x_last = reshape(x_ini',2*Nv,1);
    x_left = zeros(3*N*Nv-3*Nv,1);
    x_left_copy = zeros(3*Nv,1);
    
    % reserve some space for trajectories
    xux_all = cell(Nv,2*NK+1);
    x_trj = zeros(2*Nv,NK+1);
    appl_trj = zeros(2*Nv,1);
    pxy_trj = zeros(2*Nv,NK+1);
    
    px_trj = zeros(Nv,NK+1);
    py_trj = zeros(Nv,NK+1);
    p_trj = zeros(Nv,NK+1);
    p_trj(:,1) = x_ini(:,1);
    v_trj = zeros(Nv,NK+1);
    v_trj(:,1) = x_ini(:,2);
    a_trj = zeros(Nv,NK);
    
    for i = 1:Nv
        pv_temp = x_ini(i,:)';
        xux_all{i,1} = pv_temp;
        xy_temp = full([F_cell{i,1}(pv_temp(1)),F_cell{i,2}(pv_temp(1))]);
        pxy_trj(2*i-1:2*i,1) = xy_temp';
        px_trj(i,1) = full(F_cell{i,1}(p_trj(i,1)));
        py_trj(i,1) = full(F_cell{i,2}(p_trj(i,1)));
    end
    x_trj(:,1) = x_last;
    
    %% Plot: Ready
    export_gif = true;
    export_pdf = false;
    gif_filename = 'MPC_trj.gif';
    
    color_set = ['r';'g';'b';'c';'m';'y';'k';];
    color_pred_set = [':.r';':.g';':.b';':.c';':.m';':.y';':.k'];
    
    gatherplt = figure('Name','Trajectory Evolution', 'NumberTitle','off', 'PaperPositionMode', 'auto', 'PaperUnits', 'points', 'PaperSize', [825 255], 'DefaultAxesFontSize',14);

    trjplt = figure('Name','Trajectory', 'NumberTitle','off', 'PaperPositionMode', 'auto', 'PaperUnits', 'points', 'PaperSize', [260 260],'DefaultAxesFontSize',14);
    set(gca, 'FontName', 'Times New Roman', 'XTick', 10:10:20, 'YTick', 10:10:20);

    hold on
    axis([-0.5*L 0.5*L -0.5*L 0.5*L])
    axis square
    
    road_x = [0.5*L; 0.5*Ls;0.5*Ls;-0.5*Ls;-0.5*Ls;-0.5*L; -0.5*L;-0.5*Ls;-0.5*Ls;0.5*Ls; 0.5*Ls;0.5*L];
    road_y = [0.5*Ls;0.5*Ls;0.5*L ; 0.5*L ; 0.5*Ls;0.5*Ls;-0.5*Ls;-0.5*Ls; -0.5*L;-0.5*L;-0.5*Ls;-0.5*Ls];
    fill(road_x,road_y,[0.5 0.5 0.5],'FaceAlpha', 0.2, 'LineStyle', 'none');

    ylabel('$$p_y \; \textrm{[m]}$$','Interpreter','latex')
    xlabel('$$p_x \; \textrm{[m]}$$','Interpreter','latex')
    
    cartrj_cell = cell(Nv,1);
    recs_cell = cell(Nv,1);
    predtrj_cell = cell(Nv,1);
    
    for i = 1:Nv
        cartrj_cell{i,1} = animatedline('Color',color_set(i,:),'LineWidth',1.5);
        [rec_x,rec_y] = ellip_gen(pxy_trj(2*i-1,1),pxy_trj(2*i,1),0,car_ls,car_ss);
        recs_cell{i,1} = fill(rec_x,rec_y,color_set(i,:),'FaceAlpha', 0.5, 'LineStyle', 'none');
        addpoints(cartrj_cell{i,1},pxy_trj(2*i-1,1),pxy_trj(2*i,1));
    end
     
    
    %% MPC loop
    for K = 1:NK
        
        % reserve some space for w and w0
        w = cell(Nv,2*N+1);
        lbw = zeros(3*Nv*N+2*Nv,1);
        ubw = zeros(3*Nv*N+2*Nv,1);
        g = {};
        lbg = [];
        ubg = [];
        J = 0;
        
        % initial conditions       
        w(:,1) = MXcell(1,Nv,['X' int2str(K-1) '_'],2);
        w0 = [x_last;x_left;x_left_copy];
        % w0 = [x_last;zeros(3*N*Nv,1)];
        % w0 = zeros(3*N*Nv+2*Nv,1);
        lbw(1:2*Nv) = x_last;
        ubw(1:2*Nv) = x_last;
        
        % construct lbw and ubw
        for k = 1:N
            lbw(3*Nv*k-Nv+1:3*Nv*k) = a_min*ones(Nv,1);
            ubw(3*Nv*k-Nv+1:3*Nv*k) = a_max*ones(Nv,1);
            lbw(3*Nv*k+1:2:3*Nv*k+2*Nv-1) = -inf;
            ubw(3*Nv*k+1:2:3*Nv*k+2*Nv-1) = +inf;
            lbw(3*Nv*k+2:2:3*Nv*k+2*Nv) = v_min;
            ubw(3*Nv*k+2:2:3*Nv*k+2*Nv) = v_max;
        end
            
        %% collision avoidance priciples             
        
        I_set = [];  % set of inside vehicles' index
        R_set = cell(8,1);  % set of outside vehicles' index, cell size = 8*1
        margin_set = [];
        I_col_set = [];
        R_col_set = cell(8,1);
        
        % set of vehicles Inside
        % I_col_set
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
        
        if size(I_set,2) > 1
            I_col_set = nchoosek(I_set,2);
        end
        
        % set of vehicels in the margin
        % R_col_set
        % IR_col_set
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
        
        %% Predict horizon
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
                % Objective of reference
                J = J + 0.2*w{h,2*cnt}^2 + (w{h,2*cnt+1}(2)-v_ref)^2;
            end
            
            % inside collision avoidance
            I_len = size(I_col_set,1);
            for i = 1:I_len
                a_index = I_col_set(i,1);
                b_index = I_col_set(i,2);
                a_p = w{a_index,2*cnt+1}(1);
                b_p = w{b_index,2*cnt+1}(1);
                a_cor = [F_cell{a_index,1}(a_p);F_cell{a_index,2}(a_p)];
                b_cor = [F_cell{b_index,1}(b_p);F_cell{b_index,2}(b_p)];
                g = [g {norm(a_cor-b_cor)}];
                lbg = [lbg; Di];
                ubg = [ubg; +inf];
            end
            
            % outside collision avoidance
            R_col_mat = cell2mat(R_col_set);
            R_len = size(R_col_mat,1);
            if k >= K+4
                for i = 1:R_len
                    a_index = R_col_mat(i,1);
                    b_index = R_col_mat(i,2);
                    % g = [g { w{a_index,2*cnt+1}(1)-w{b_index,2*cnt-1}(1) }];
                    a_p = w{a_index,2*cnt+1}(1);
                    b_p = w{b_index,2*cnt+1}(1);
                    a_cor = [F_cell{a_index,1}(a_p);F_cell{a_index,2}(a_p)];
                    b_cor = [F_cell{b_index,1}(b_p);F_cell{b_index,2}(b_p)];
                    g = [g {norm(a_cor-b_cor)}];
                    lbg = [lbg; Do];
                    ubg = [ubg; +inf];
                end
            end
            
            % in and out collision avoidance
            IR_len = size(IR_col_set,1);
            for i = 1:IR_len
                a_index = IR_col_set(i,1);
                b_index = IR_col_set(i,2);
                a_p = w{a_index,2*cnt+1}(1);
                b_p = w{b_index,2*cnt+1}(1);
                a_cor = [F_cell{a_index,1}(a_p);F_cell{a_index,2}(a_p)];
                b_cor = [F_cell{b_index,1}(b_p);F_cell{b_index,2}(b_p)];
                g = [g {norm(a_cor-b_cor)}];
                lbg = [lbg; Di];
                ubg = [ubg; +inf];
            end       
        end
        
        % create an NLP solver
        prob = struct('f',J,'x',vertcat(w{:}),'g',vertcat(g{:}));
        setting.ipopt.print_level = 5;
        solver = nlpsol('solver', 'ipopt', prob, setting);
        
        % solve the NLP
        sol = solver('x0',w0,'lbx',lbw,'ubx',ubw,'lbg',lbg,'ubg',ubg);
        K
        w_opt = full(sol.x);
        

        for i = 1:Nv
            xux_all{i,2*K} = w_opt(3*Nv-Nv+i); %u0 ?
            xux_all{i,2*K+1} = w_opt(3*Nv+2*i-1: 3*Nv+2*i); %x1
            appl_trj(2*i-1:2*i,1) = w_opt(3*Nv+2*i-1: 3*Nv+2*i);
        end
        
        x_left = w_opt(5*Nv+1:end);
        x_left_copy = x_left(end-3*Nv+1:end);
        x_last = appl_trj; 

        p_trj(:,K+1) = x_last(1:2:end);
        v_trj(:,K+1) = x_last(2:2:end);
        a_trj(:,K) = w_opt(2*Nv+1:3*Nv);
        for i = 1:Nv
            px_trj(i,K+1) = full(F_cell{i,1}(p_trj(i,K+1)));
            py_trj(i,K+1) = full(F_cell{i,2}(p_trj(i,K+1)));
        end
            
        for i = 1:Nv
            x_trj(:,K+1) = x_last;
            xy_temp = full([F_cell{i,1}(x_last(2*i-1)),F_cell{i,2}(x_last(2*i-1))]);
            pxy_trj(2*i-1:2*i,K+1) = xy_temp';  
        end
        
        pxy_pred = zeros(2*Nv,N+1);
        for i = 1:N+1
            for j = 1:Nv
                 pxy_pred(2*j-1,i) = full(F_cell{j,1}( w_opt( 3*Nv*(i-1)+2*j-1 ,1) ));
                 pxy_pred(2*j,i) = full(F_cell{j,2}( w_opt( 3*Nv*(i-1)+2*j-1 ,1) ));
                
            end
        end
        
        %% Plot: Step by step
        for i = 1:Nv
            angle = atan( (pxy_trj(2*i,K+1)-pxy_trj(2*i,K))/(pxy_trj(2*i-1,K+1)-pxy_trj(2*i-1,K)) );
            [rec_x, rec_y] = ellip_gen(pxy_trj(2*i-1,K+1),pxy_trj(2*i,K+1),angle,0.8*car_ls,0.8*car_ss);
            set(recs_cell{i,1},'XData',rec_x,'YData',rec_y);
            addpoints(cartrj_cell{i,1},pxy_trj(2*i-1,K+1),pxy_trj(2*i,K+1));
        end
        
        if K == 1
            for i = 1:Nv
                predtrj_cell{i,1} = plot(pxy_pred(2*i-1,:), pxy_pred(2*i,:), color_pred_set(i,:));
            end
            
            figcaptext = text(8, 8, sprintf('MPC loop %d', K), 'Interpreter', 'latex');
        else
            
            for i = 1:Nv
                set(predtrj_cell{i,1},'XData',pxy_pred(2*i-1,:),'YData',pxy_pred(2*i,:));
            end
            
            set(figcaptext, 'String', sprintf('MPC loop %d', K));
        end
            
        drawnow
            
        if export_gif
            frame = getframe(trjplt);
            im = frame2im(frame);
            [imind,cm] = rgb2ind(im,256);
            % Write to the GIF File 
            if K == 1 
              imwrite(imind,cm,gif_filename,'gif','Loopcount',inf); 
            else
              imwrite(imind,cm,gif_filename,'gif','WriteMode','append'); 
            end
        end
        
        
    end

    
    % plot
    v_p = x_trj(2:2:end,:); % 7* NK+1
    a_p = xux_all(:,2:2:end); a_p = vertcat(a_p{:}); a_p = reshape(a_p,Nv,NK); % 7*NK
    t_v = [0:T/N:T/N*NK];
    t_a = [0:T/N:T/N*(NK-1)];
    style_set = ['-.r';'-.g';'-.b';'-.c';'-.m';'-.y';'-.k';];
    
    figure(3);
    for i = 1:Nv
        plot(t_v,v_p(i,:),style_set(i));
        hold on;
    end
    title('velocity');
    ylabel('v-m/s');
    xlabel('t-s');
    
    figure(4);
    for i = 1:Nv
        plot(t_a,a_p(i,:),style_set(i));
        hold on;
    end
    title('acceleration');
    ylabel('a-m/s^2');
    
    xlabel('t-s');

    save('trjs','p_trj','px_trj','py_trj', 'v_trj', 'a_trj');
    output = {p_trj, v_trj, a_trj};
end
    
