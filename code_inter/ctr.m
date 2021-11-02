close all
clear
clc

addpath('D:/Control_System/casadi')
import casadi.*

% ODE
x = SX.sym('x',4);
u = SX.sym('u',2);
dx = [x(3)*cos(x(4));...
      x(3)*sin(x(4));...
      u(1);u(2)];
%para
T = 1;
N = 20;
NK = 64;
Ts = T/N;
 
x1trj = []; x2trj = []; x3trj = []; x4trj = [];
predicted1trj = cell(NK,1); predicted2trj = cell(NK,1);
predicted3trj = cell(NK,1); predicted4trj = cell(NK,1);
applied1trj = []; applied2trj = []; 
applied3trj = []; applied4trj = [];

dae = struct('x',x,'p',u,'ode',dx);
opts = struct('tf',Ts);
F = casadi.integrator('F','cvodes',dae,opts);

%% modeling
x1 = SX.sym('x1',4); x2 = SX.sym('x2',4); x3 = SX.sym('x3',4); x4 = SX.sym('x4',4);

% constraints
cons1 = norm(x1(1:2)-x2(1:2));
cons2 = norm(x1(1:2)-x3(1:2));
cons3 = norm(x2(1:2)-x3(1:2));
cons4 = norm(x3(1:2)-x4(1:2));
G = Function('col_avoid',{x1,x2,x3,x4},{[cons1;conos2;cons3;cons4]});

% initial
initial_x1 = ; initial_x2 = ;
initial_x3 = ; initial_x4 = ;

% boundary
lbu = [-5;-5];
ubu = [5;5];

lbx = [-inf;-inf;5;-inf];
ubx = [inf;inf;20;inf];

% initial conditions
initial_x1 = []; initial_x2 = [];
initial_x3 = []; initial_x4 = [];

lastX1 = initial_x1; lastX2 = initial_x2;
lastX3 = initial_x3; lastX4 = initial_x4;
lastPhi1 = 0; lastPhi2 = 0;
lastPhi3 = 0; lastPhi4 = 0;

% MPC
for k=1:NK
    
    w1 = {}; w01 = []; lbw1 = []; ubw1 = []; 
    g1 = {}; lbg1 = []; ubg1 = [];
    w2 = {}; w02 = []; lbw2 = []; ubw2 = []; 
    g2 = {}; lbg2 = []; ubg2 = [];
    w3 = {}; w03 = []; lbw3 = []; ubw3 = []; 
    g3 = {}; lbg3 = []; ubg3 = [];
    w4 = {}; w04 = []; lbw4 = []; ubw4 = []; 
    g4 = {}; lbg4 = []; ubg4 = [];
    w = {}; w0 = []; lbw = []; ubw = []; 
    g = {}; lbg = []; ubg = [];
    J = 0;
    
    % Lift initial conditions
    X10 = MX.sym(['X1_' int2str(K-1)], 4);
    w1 = [w1 {X10}];
    w01 = [w01; lastX1];
    lbw1 = [lbw1; lastX1];
    ubw1 = [ubw1; lastX1];
    
    X20 = MX.sym(['X2_' int2str(K-1)], 4);
    w2 = [w2 {X20}];
    w02 = [w02; lastX2];
    lbw2 = [lbw2; lastX2];
    ubw2 = [ubw2; lastX2];
    
    X30 = MX.sym(['X3_' int2str(K-1)], 4);
    w3 = [w3 {X30}];
    w03 = [w03; lastX3];
    lbw3 = [lbw3; lastX3];
    ubw3 = [ubw3; lastX3];
    
    X40 = MX.sym(['X4_' int2str(K-1)], 4);
    w4 = [w4 {X40}];
    w04 = [w04; lastX4];
    lbw4 = [lbw4; lastX4];
    ubw4 = [ubw4; lastX4];
    
    % Formulate the NLP
    X1k = X10; X2k = X20; X3k = X30; X4k = X40;
    phi1 = lastPhi1; phi2 = lastPhi2;
    phi3 = lastPhi3; phi4 = lastPhi4;
    
    for k = 0+K-1:N+K-2
        U1k = MX.sym(['U1_' num2str(k)], 2);
        w1 = [w1 {U1k}];
        w01 = [w01;zeros(2,1)];
        lbw1 = [lbw1; lbu];
        ubw1 = [ubw1; ubu];
        J = J + U1k'*R*U1k; 
        
        U2k = MX.sym(['U2_' num2str(k)], 2);
        w2 = [w2 {U2k}];
        w02 = [w02;zeros(2,1)];
        lbw2 = [lbw2; lbu];
        ubw2 = [ubw2; ubu];
        J = J + U1k'*R*U1k;
        
        U3k = MX.sym(['U3_' num2str(k)], 2);
        w3 = [w3 {U3k}];
        w03 = [w03;zeros(2,1)];
        lbw3 = [lbw3; lbu];
        ubw3 = [ubw3; ubu];
        J = J + U3k'*R*U3k;
        
        U4k = MX.sym(['U1_' num2str(k)], 2);
        w4 = [w4 {U4k}];
        w04 = [w04;zeros(2,1)];
        lbw4 = [lbw4; lbu];
        ubw4 = [ubw4; ubu];
        J = J + U4k'*R*U4k;
        
        % Integrate till the end of the interval
        F1k = F('x0', X1k, 'p', U1k);
        X1k_end = F1k.xf;
        F2k = F('x0', X2k, 'p', U2k);
        X2k_end = F2k.xf;
        F3k = F('x0', X3k, 'p', U3k);
        X3k_end = F3k.xf;
        F4k = F('x0', X4k, 'p', U4k);
        X4k_end = F4k.xf;        
        
        % New NLP variable for state at the end of interval
        phi1 = phi1 + X1k(3)*Ts; phi2 = phi2 + X2k(3)*Ts;
        phi3 = phi3 + X1k(3)*Ts; phi4 = phi4 + X4k(3)*Ts;
        % ???
        
    end
end    