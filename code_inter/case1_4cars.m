close all
clear
clc

addpath('D:\Control_System\casadi')
import casadi.*

% % description
% % state: x = (p,v)
% % control: u = a

% ODE
x = SX.sym('x',2);
u = SX.sym('u',1);
dx = [x(2);u];

% parameter
T = 0.5; N = 40;
NK = 1; Ts  = T/N;

x1trj = []; x2trj = []; x3trj = []; x4trj = [];
predicted1trj = cell(NK,1); predicted2trj = cell(NK,1);
predicted3trj = cell(NK,1); predicted4trj = cell(NK,1);
applied1trj = []; applied2trj = []; 
applied3trj = []; applied4trj = [];

dae = struct('x',x,'p',u,'ode',dx);
opts = struct('tf',Ts);
F = casadi.integrator('F','cvodes',dae,opts);

% modeling
x1 = SX.sym('x1',2); x2 = SX.sym('x2',2); 
x3 = SX.sym('x3',2); x4 = SX.sym('x4',2);

% boundary
lbx = [0;1]; ubx = [inf;20];
lbu = -3; ubu = 3;

% initial state p=0 v=v_ref=13
initial_x1 = [0;13]; initial_x2 = [0;13];
initial_x3 = [0;13]; initial_x4 = [0;13];

lastX1 = initial_x1; lastX2 = initial_x2;
lastX3 = initial_x3; lastX4 = initial_x4;

% % MPC
for K = 1:NK
    % Start with an empty NLP
    w1 = {}; w01 = []; lbw1 = []; ubw1 = []; 
    g1 = {}; lbg1 = []; ubg1 = [];
    w2 = {}; w02 = []; lbw2 = []; ubw2 = []; 
    g2 = {}; lbg2 = []; ubg2 = [];
    w3 = {}; w03 = []; lbw3 = []; ubw3 = []; 
    g3 = {}; lbg3 = []; ubg3 = [];
    w4 = {}; w04 = []; lbw4 = []; ubw4 = []; 
    g4 = {}; lbg4 = []; ubg4 = []; 
    J = 0;
    
    % Lift initial conditions
    X10 = MX.sym(['X1_' int2str(K-1)], 2);
    w1 = [w1 {X10}];
    w01 = [w01; lastX1];
    lbw1 = [lbw1; lastX1];
    ubw1 = [ubw1; lastX1];
    
    X20 = MX.sym(['X2_' int2str(K-1)], 2);
    w2 = [w2 {X20}];
    w02 = [w02; lastX2];
    lbw2 = [lbw2; lastX2];
    ubw2 = [ubw2; lastX2];
    
    X30 = MX.sym(['X3_' int2str(K-1)], 2);
    w3 = [w3 {X30}];
    w03 = [w03; lastX3];
    lbw3 = [lbw3; lastX3];
    ubw3 = [ubw3; lastX3];
    
    X40 = MX.sym(['X4_' int2str(K-1)], 2);
    w4 = [w4 {X40}];
    w04 = [w04; lastX4];
    lbw4 = [lbw4; lastX4];
    ubw4 = [ubw4; lastX4];
    
    % Formulate the NLP
    X1k = X10; X2k = X20; X3k = X30; X4k = X40;
    
    % Do some judge and collision avoidance constraint    
    
    % % Construct new horizon prediction
    for k = 0+K-1:N+K-2
        
        % New NLP variable for the control
        U1k = MX.sym(['U1_' num2str(k)],1);
        w1 = [w1 {U1k}]; w01 = [w01;0];
        lbw1 = [lbw1;lbu]; ubw1 = [ubw1;ubu];
        
        U2k = MX.sym(['U2_' num2str(k)],1);
        w2 = [w2 {U2k}]; w02 = [w02;0];
        lbw2 = [lbw2;lbu]; ubw2 = [ubw2;ubu];
        
        U3k = MX.sym(['U3_' num2str(k)],1);
        w3 = [w3 {U3k}]; w03 = [w03;0];
        lbw3 = [lbw3;lbu]; ubw3 = [ubw3;ubu];
        
        U4k = MX.sym(['U4_' num2str(k)],1);
        w4 = [w4 {U4k}]; w04 = [w04;0];
        lbw4 = [lbw4;lbu]; ubw4 = [ubw4;ubu];        
        
        J = J + U1k'*U1k + U2k'*U2k + U3k'*U3k + U4k'*U4k;
        
        % Integrate till the end of the interval
        F1k = F('x0',X1k,'p',U1k);
        X1k_end = F1k.xf;
        
        F2k = F('x0',X2k,'p',U2k);
        X2k_end = F2k.xf;
        
        F3k = F('x0',X3k,'p',U3k);
        X3k_end = F3k.xf;
        
        F4k = F('x0',X4k,'p',U4k);
        X4k_end = F4k.xf;        
        
        X1k = MX.sym(['X1_' num2str(k+1)], 2);
        w1 = [w1 {X1k}];
        w01 = [w01; 10*(Ts*(k+1)); 10];
        lbw1 = [lbw1; lbx];
        ubw1 = [ubw1; ubx];        

        X2k = MX.sym(['X2_' num2str(k+1)], 2);
        w2 = [w2 {X2k}];
        w02 = [w02; 10*(Ts*(k+1)); 10];
        lbw2 = [lbw2; lbx];
        ubw2 = [ubw2; ubx];
        
        X3k = MX.sym(['X3_' num2str(k+1)], 2);
        w3 = [w3 {X3k}];
        w03 = [w03; 10*(Ts*(k+1)); 10];
        lbw3 = [lbw3; lbx];
        ubw3 = [ubw3; ubx]; 
        
        X4k = MX.sym(['X4_' num2str(k+1)], 2);
        w4 = [w4 {X4k}];
        w04 = [w04; 10*(Ts*(k+1)); 10];
        lbw4 = [lbw4; lbx];
        ubw4 = [ubw4; ubx]; 
        
        pos1k = [2;4-X1k(1)]; pos2k = [-3+X2k(1);-2];
%         pos3k = [-7+X3k(1);-2]; pos4k = [-10+X4k(1);-2];
        cons1k = (pos1k-pos2k)'*(pos1k-pos2k)-0.2;
%         cons2k = norm(pos1k-pos3k) - 1.5;
%         cons3k = norm(pos2k-pos3k) - 1.5;
%         cons4k = norm(pos3k-pos4k) - 1.5;
        J = J - 1.0e-5 * (log(cons1k));... + log(cons2k) + log(cons3k) + log(cons4k));
        
        % Add equality constraints to make x continous
        g1 = [g1 {X1k_end-X1k}];
        lbg1 = [lbg1; zeros(2,1)];
        ubg1 = [ubg1; zeros(2,1)]; 
        
        g2 = [g2 {X2k_end-X2k}];
        lbg2 = [lbg2; zeros(2,1)];
        ubg2 = [ubg2; zeros(2,1)]; 
        
        g3 = [g3 {X3k_end-X3k}];
        lbg3 = [lbg3; zeros(2,1)];
        ubg3 = [ubg3; zeros(2,1)]; 
        
        g4 = [g4 {X4k_end-X4k}];
        lbg4 = [lbg4; zeros(2,1)];
        ubg4 = [ubg4; zeros(2,1)];         
        
    end
    
    % % Create an NLP solver
    prob = struct('f', J, 'x', vertcat(w1{:},w2{:},w3{:},w4{:}),...
                  'g', vertcat(g1{:},g2{:},g3{:},g4{:}));
    setting.ipopt.print_level = 5;
    solver = nlpsol('solver', 'ipopt', prob,setting);
    
    % Solve the NLP
    sol = solver('x0', [w01;w02;w03;w04],...
                 'lbx', [lbw1;lbw2;lbw3;lbw4], 'ubx', [ubw1;ubw2;ubw3;ubw4],...
                 'lbg', [lbg1;lbg2;lbg3;lbg4], 'ubg', [ubg1;ubg2;ubg3;ubg4]);
    w_opt = full(sol.x);
    
    % save optimized predicted trajectory
    result1 = w_opt(1:2+3*N);
    result2 = w_opt(3+3*N:4+2*3*N);
    result3 = w_opt(5+2*3*N:6+3*3*N);
    result4 = w_opt(7+3*3*N:end);
    
    state_control1 = reshape(result1(3:end),3,N);
    predicted1trj{K} = state_control1;
    state_control2 = reshape(result2(3:end),3,N);
    predicted2trj{K} = state_control2;
    state_control3 = reshape(result3(3:end),3,N);
    predicted3trj{K} = state_control3;
    state_control4 = reshape(result4(3:end),3,N);
    predicted4trj{K} = state_control4;
    
    %%
    p1_trj = state_control1(2,:);
    px1_trj = 2*ones(1,N); py1_trj = -p1_trj+4;
    plot(px1_trj,py1_trj,'--rs','Linewidth',1,'MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',3);
    hold on;

    p2_trj = state_control2(2,:);
    px2_trj = p2_trj-4; py2_trj = -2*ones(1,N);
    plot(px2_trj,py2_trj,'--gs','Linewidth',1,'MarkerEdgeColor','g','MarkerFaceColor','g','MarkerSize',3);
    hold on;

%     p3_trj = state_control3(2,:);
%     px3_trj = p3_trj-7; py3_trj = -2*ones(1,50);
%     plot(px3_trj,py3_trj,'--bs','Linewidth',1,'MarkerEdgeColor','b','MarkerFaceColor','b','MarkerSize',3);
%     hold on;
% 
%     p4_trj = state_control4(2,:);
%     px4_trj = p4_trj-10; py4_trj = -2*ones(1,50);
%     plot(px4_trj,py4_trj,'--ks','Linewidth',1,'MarkerEdgeColor','k','MarkerFaceColor','k','MarkerSize',3);
    %%
    newX1 = state_control1(2:3,1);
    newX2 = state_control2(2:3,1);
    newX3 = state_control3(2:3,1);        
    newX4 = state_control4(2:3,1);
    
    % updates
    
    lastX1 = newX1;
    lastX2 = newX2;
    lastX3 = newX3;
    lastX4 = newX4;
    
    applied1trj = [applied1trj state_control1(:,1)];
    applied2trj = [applied2trj state_control2(:,1)];
    applied3trj = [applied3trj state_control3(:,1)];
    applied4trj = [applied4trj state_control4(:,1)];
    
    x1trj = [x1trj lastX1];
    x2trj = [x2trj lastX2];
    x3trj = [x3trj lastX3];
    x4trj = [x4trj lastX4];
    
end

