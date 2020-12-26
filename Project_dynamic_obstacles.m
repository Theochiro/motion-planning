%% Notes about the code
%Author: Theodor Chakhachiro
%Convex Optimization Project
%% Comment
%User is able to add/remove as many obstacles as desired by changing the 
%Obs matrix below. Also, the starting position and the goal pose can be
%modified. Finally, the number of sides of the polygon can also be changed
%to further constrict/relax the constraints.
%Code will also create a gif of the actual robot trajectory.
%% Innitializing
clear all
close all
clc
filename = 'testAnimated.gif';
%% Required Constants
h=10;              %Length of Control Horizon 
r_e=0.04;          %Radius embracing the perimiter of the robot
v_max=0.05;        %Maximum Allowed Velocity in m/s
s=6;               %Polygon dimention to convefiy the obstacle
M=2e2;             %From Big-M Method constant
tol=1e-1;          %Tolerance
z_s=[0 0];         %Coordinates of initial position
z_g=[1 1];         %Coordinates of final position
Weights=[1 10 1];  %Weights of the objective function (2 10 1)
Z_S=z_s;      
Total_Length=0;
X_aux=[];
Y_aux=[];
OBS=[];
%% Adding Obstacles
%Obs(i,:) = [x_obs y_obs radius_obs] in m
% Obs=[0.15 0.25 0.1;0.6 0.5 0.15];                         %2 Obstacles
% Obs=[0.6 0.5 0.15];                                       %1 Obstacle
% Obs=[0.1 0.2 0.1;0.4 1 0.15; 0.6 0.5 0.2;1 0.2 0.15];     %4 Obtacles
%% Main Code
T=0;
theta_e0=0;
dz=norm(z_s-z_g);
tic;
nnn=1;
while dz>tol
    %% Dynamic Obstacle
Obs=[-0.25*cos(0.15*T-2)+0.7,0.25*sin(0.15*T-2)+0.4,0.05];
OBS=[OBS;Obs];
[r,c]=size(Obs);
%% Obstacle circumscription & Setting up the plot
for fig=1:2
    figure(fig);
    hold on
    f=1:s;
    b=[];
    A=[];
    for j=1:r
        for i=1:s
            phi(i)=2*pi*i/s;
            a_t(i)=tan(-phi(i)+pi/2-pi/s);
            x_aux(i)=(Obs(j,3)/cos(pi/s))*cos(-phi(i))+Obs(j,1);
            y_aux(i)=(Obs(j,3)/cos(pi/s))*sin(-phi(i))+Obs(j,2);
            b_t(i)=y_aux(i)-a_t(i)*x_aux(i);                
        end
        X_aux=[x_aux x_aux(1)]';
        Y_aux=[y_aux y_aux(1)]';
        A=[A;[(diff(Y_aux)./diff(X_aux)) ones(length(x_aux),1)]];
        b=[b;y_aux'-A((j-1)*s+1:s*j,1).*x_aux'];
        circle(Obs(j,1), Obs(j,2), Obs(j,3), 'r');
        v1 = [x_aux;y_aux]';
        patch('Faces',f,'Vertices',v1,...
            'EdgeColor','black','FaceColor','none','LineWidth',2);
    end
    axis([-0.2 1.2 -0.2 1.2])
    xlabel("$$x \ (m)$$",'interpreter', 'latex')
    ylabel("$$y \ (m)$$",'interpreter', 'latex')
    if fig==1
        title("\textbf{Actual Robot Trajectory}",'interpreter', 'latex')
    else
        title("\textbf{Motion Planning}",'interpreter', 'latex')
    end
    hold off
end
A(A(:,1)>0&A(:,1)<1e-3)=0;
for n=1:r
    for m=1+6*(n-1):s*n
        if Obs(n,2)>(A(m,1)*Obs(n,1)+b(m))
            A(m,2)=-A(m,2);
            b(m)=-b(m);
        elseif Obs(n,2)<(A(m,1)*Obs(n,1)+b(m))
            A(m,1)=-A(m,1);
        end
    end
end
    [r1,c1]=size(Z_S);
    figure(1);
    hold on
    for i=1:r1
        circle(Z_S(i,1), Z_S(i,2), r_e, [1-i/r1 1-i/r1 1]);
    end
    axis([-0.2 1.2 -0.2 1.2])
    cvx_begin  
        cvx_solver mosek
        variables z_e(h,2) dt(h-1,1)
        variable v(r*s,1) binary
        summation=0;
        v_sum=[];
        max_const=[];
        obs_right=[];
        obs_left=[];
        for j=1:h-1
            summation=summation+pow_pos(norm(z_e(j+1,:)-z_e(j,:)),2);
            max_const=[max_const;norm(z_e(j+1,:)-z_e(j,:))];
            obs_right=[obs_right;A(:,1)*(z_e(j,1)')+A(:,2)*(z_e(j,2))'];
            obs_left=[obs_left;b+(v-1)*M];
        end
        for k=1:r
            v_sum=[v_sum;sum(v((k-1)*s+1:s*k))];
        end
        minimize(Weights(1)*sum(dt)+Weights(2)*summation+Weights(3)*pow_pos(norm(z_e(h,:)-z_g),2))
        subject to
            obs_right >= obs_left ;
            v_sum == ones(r,1);
            z_e(1,:)==z_s;
            max_const <= v_max*dt;
            dt>=2;

    cvx_end
    d_ze=diff(z_e);
    if d_ze(1,1)>=0
        theta_e=atan2(d_ze(1,2),d_ze(1,1));
    else
        theta_e=atan2(d_ze(1,2),d_ze(1,1))+pi;
    end
    vel=norm(d_ze(1,:))/dt(1,1);
    w=(theta_e-theta_e0)/dt(1);
    z_s=z_s+[vel*cos(theta_e)*dt(1,1) vel*sin(theta_e)*dt(1,1)];
    theta_e0=theta_e0+w*dt(1,1);
    dz=norm(z_s-z_g);
    Z_S=[Z_S;z_s];
    T=T+dt(1,1);
    [r_obs,c_obs]=size(OBS);
    hh=figure(1);
    hold off
    figure(1);
    hold on
    for i=1:r_obs
        circle(OBS(i,1), OBS(i,2), OBS(i,3), [1 1-i/r_obs 1-i/r_obs]);
    end
    axis([-0.2 1.2 -0.2 1.2])
              frame = getframe(hh); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if nnn == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end
    hold off
    figure(2);
    hold on
    for i=1:h
        circle(z_e(i,1), z_e(i,2), r_e, [1-i/h 1-i/h 1]);
    end
    for i=1:r_obs
        circle(OBS(i,1), OBS(i,2), OBS(i,3), [1 1-i/r_obs 1-i/r_obs]);
    end
    axis([-0.2 1.2 -0.2 1.2])
    hold off
    nnn=nnn+1;
end
Solving_Time=toc;
Total_Length=sum(sqrt(sum(diff(Z_S).^2,2)));
fprintf("\nSolving Time: %f seconds\nTime Taken: %f seconds\nTrajectory Length: %f m \nFinal Robot Pose: [%f,%f]\n", ...
Solving_Time,T,Total_Length,z_s(1,1),z_s(1,2));
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
