
function output = ellip_cons(rx,ry,cor1,the1,cor2,the2,delta)
    
    Q1_11 = (cos(the1)^2)*(rx^2) + (sin(the1)^2)*(ry^2);
    Q1_12 = 0.5*sin(2*the1)*(rx^2-ry^2);
    Q1_21 = Q1_12;
    Q1_22 = (sin(the1)^2)*(rx^2) + (cos(the1)^2)*(ry^2);
    Q1 = [Q1_11,Q1_12;Q1_21,Q1_22];
    
    Q2_11 = (cos(the2)^2)*(rx^2) + (sin(the2)^2)*(ry^2);
    Q2_12 = 0.5*sin(2*the2)*(rx^2-ry^2);
    Q2_21 = Q2_12;
    Q2_22 = (sin(the2)^2)*(rx^2) + (cos(the2)^2)*(ry^2);
    Q2 = [Q2_11,Q2_12;Q2_21,Q2_22];
    
    xi_1 = ((eye(2)+delta(2)*Q2)*(eye(2)+delta(1)*Q1)-eye(2))\ ...   
       (delta(1)*(eye(2)+delta(2)*Q2)*Q1*cor1+delta(2)*Q2*cor2);
   
    xi_2 = ((eye(2)+delta(1)*Q1)*(eye(2)+delta(2)*Q2)-eye(2))\ ...
       (delta(2)*(eye(2)+delta(1)*Q1)*Q2*cor2+delta(1)*Q1*cor1);

    Llam1 = (xi_1-cor1)'* Q1 * (xi_1-cor1) - 1;
    Llam2 = (xi_2-cor2)'* Q2 * (xi_2-cor2) - 1;

    cons1 = norm(xi_1 - xi_2); % >ds_in
    cons2 = [Llam1;Llam2]; % = [0;0]
    
    output = [cons1;cons2];

end