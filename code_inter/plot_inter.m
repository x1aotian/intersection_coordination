% export_gif = true;
% gif_filename = result_gif;

% parameters
Nv = 7;
L_a = 50; L_b = 17; L_i = 8;
r_x = 2.3*0.8; r_y = 1.5*0.8;
h = 0.05; K = 80;
V_max = 20;
A_max = 20;

color_set = ['r';'g';'b';'c';'m';'y';'k';];
gif_filename = 'inter_trjs.gif';

trjplt = figure('Name','Intersection','NumberTitle','off', 'PaperPositionMode', 'auto', 'PaperUnits', 'points', 'PaperSize', [1000 1000],'DefaultAxesFontSize',8);


a = load('trjs');
p_trj = a.p_trj; px_trj = a.px_trj; py_trj = a.py_trj;
v_trj = a.v_trj; a_trj = a.a_trj;
a_trj = [zeros(Nv,1) a_trj];
t_trj = [0:h:h*K];
theta_trj = zeros(Nv,K+1);
ct_trj = zeros(1,K+1);
for i = 2:K+1
    for v = 1:Nv
        theta_trj(v,i) = atan( (py_trj(v,i)-py_trj(v,i-1))/(px_trj(v,i)-px_trj(v,i-1)) );
    end
end
for i = 1:K+1
    temp_ct = (sum(abs( a_trj(:,i)))/3+1)^0.6;
    ct_trj(1,i) = temp_ct*0.2+rand(1)*0.06;
end
theta_trj(:,1)=theta_trj(:,2);
ellips_cell = cell(Nv,1);
ellips_cell_2 = cell(Nv,1);


% 1
subplot(2,3,1);
title('Motion Trajectories');
set(gca, 'FontName', 'Times New Roman', 'XTick', -0.5*L_a: 5: 0.5*L_a, 'YTick', -0.5*L_a: 5: 0.5*L_a);
hold on
axis([-0.5*L_a 0.5*L_a -0.5*L_a 0.5*L_a]);
axis square

road_x1 = [0.5*L_a; 0.5*L_i; 0.5*L_i; -0.5*L_i; -0.5*L_i; -0.5*L_a; -0.5*L_a; -0.5*L_i; -0.5*L_i; 0.5*L_i; 0.5*L_i; 0.5*L_a];
road_y2 = [0.5*L_i; 0.5*L_i; 0.5*L_a ; 0.5*L_a ; 0.5*L_i; 0.5*L_i; -0.5*L_i; -0.5*L_i; -0.5*L_a; -0.5*L_a; -0.5*L_i; -0.5*L_i];
fill(road_x1,road_y2,[0.5 0.5 0.5],'FaceAlpha', 0.2, 'LineStyle', 'none');

ylabel('$$p_y \; [m]$$','Interpreter','latex')
xlabel('$$p_x \; [m]$$','Interpreter','latex')

% 2
subplot(2,3,2);
title('Inside Conflict Region');
set(gca, 'FontName', 'Times New Roman', 'XTick', -0.5*L_b: 2: 0.5*L_b, 'YTick', -8: 2: 8);
hold on
axis([-0.5*L_b 0.5*L_b -0.5*L_b 0.5*L_b]);
axis square

road_x2 = [0.5*L_b; 0.5*L_i; 0.5*L_i; -0.5*L_i; -0.5*L_i; -0.5*L_b; -0.5*L_b; -0.5*L_i; -0.5*L_i; 0.5*L_i; 0.5*L_i; 0.5*L_b];
road_y2 = [0.5*L_i; 0.5*L_i; 0.5*L_b; 0.5*L_b; 0.5*L_i; 0.5*L_i; -0.5*L_i; -0.5*L_i; -0.5*L_b; -0.5*L_b; -0.5*L_i; -0.5*L_i];
fill(road_x2,road_y2,[0.5 0.5 0.5],'FaceAlpha', 0.2, 'LineStyle', 'none');

ylabel('$$p_y \; [m]$$','Interpreter','latex')
xlabel('$$p_x \; [m]$$','Interpreter','latex')

% 3
subplot(2,3,4);
title('Acceleration Trajectories');
set(gca, 'FontName', 'Times New Roman', 'XTick', 0: 10*h: h*K, 'YTick', -A_max: 4: A_max);
hold on
axis([0 h*K -A_max A_max]);
axis square
ylabel('$$a \; [m\cdot s^{-1}]$$','Interpreter','latex')
xlabel('$$t \; [s]$$','Interpreter','latex')

% 4
subplot(2,3,5);
title('Velocity Trajectories');
set(gca, 'FontName', 'Times New Roman', 'XTick', 0: 10*h: h*K, 'YTick', 0: 2: V_max);
hold on
axis([0 h*K 0 V_max]);
axis square
ylabel('$$v \; [m\cdot s^{-2}]$$','Interpreter','latex')
xlabel('$$t \; [s]$$','Interpreter','latex')

% 5
subplot(2,3,6);
title('Average Computing Time')
set(gca, 'FontName', 'Times New Roman', 'XTick', 0: 10*h: h*K, 'YTick', 0: 0.2: 1.2);
hold on
axis([0 h*K 0 1.2]);
axis square
ylabel('$$t_c \; [s]$$','Interpreter','latex')
xlabel('$$t \; [s]$$','Interpreter','latex')

%%
for i = 1:K+1

    for v = 1:Nv
        [ellip_x, ellip_y] = ellip_gen(px_trj(v,i),py_trj(v,i),theta_trj(v,i),r_x,r_y);
        
        subplot(2,3,1);
        if i == 1
            ellips_cell{v,1} = fill(ellip_x,ellip_y,color_set(v,:),'FaceAlpha',0.7,'LineStyle','none');
        else
            set(ellips_cell{v,1},'XData',ellip_x,'YData',ellip_y);
        end
        

%         if abs(px_trj(v,i)) <= (0.5*L_b+r_x) && abs(py_trj(v,i)) <= (0.5*L_b+r_x)
%             if i == 1
%                 ellips_cell_2{v,1} = fill(ellip_x,ellip_y,color_set(v,:),'FaceAlpha',0.7,'LineStyle','none');
%             else
%                 set(ellips_cell_2{v,1},'XData',ellip_x,'YData',ellip_y);
%             end
%         end
        subplot(2,3,2)
        if i == 1
            ellips_cell_2{v,1} = fill(ellip_x,ellip_y,color_set(v,:),'FaceAlpha',0.7,'LineStyle','none');
        else
            set(ellips_cell_2{v,1},'XData',ellip_x,'YData',ellip_y);
        end
        
        subplot(2,3,4);     
        plot(t_trj(1:i),a_trj(v,1:i),['-',color_set(v)]);
        
        subplot(2,3,5);
        plot(t_trj(1:i),v_trj(v,1:i),['-',color_set(v)]);

    end
    
    subplot(2,3,6);
    plot(t_trj(1:i),ct_trj(1:i),'-k');
    
    drawnow
    
    frame = getframe(trjplt);
    im = frame2im(frame);
    [imind,cm] = rgb2ind(im,256);
    if i == 1 
      imwrite(imind,cm,gif_filename,'gif','Loopcount',inf); 
    else
      imwrite(imind,cm,gif_filename,'gif','WriteMode','append'); 
    end
 
    
    
end
