close all
clear all

addpath('../');
show_3D=0;

path = './';
V = zeros(256,440,400);


centre = [128 220 30];

R = 130;

r = 36;

for i=1:256
    for j=1:440
        for k=1:400
            if (((sqrt((k-centre(3))^2 + (j-centre(2))^2) - R)^2 + (i-centre(1))^2) <=   r^2) 
                V(i,j,k) = 1.;
            end
        end
    end
end

% right brachiocephalic artery


centre1 = [0 170 180];

R = 110;

r = 12;

for i=1:256
    for j=1:440
        for k=160:400
            if (((sqrt((k-centre1(3))^2 + (i-centre1(1))^2) - R)^2 + (j-centre1(2))^2) <=   r^2) 
                V(i,j,k) = 1.;
            end
        end
    end
end

% left common carotic artery


centre2 = [485 220 180];

R = 350;

r = 9;

for i=1:256
    for j=1:440
        for k=160:400
            if (((sqrt((k-centre2(3))^2 + (i-centre2(1))^2) - R)^2 + (j-centre2(2))^2) <=   r^2) 
                V(i,j,k) = 1.;
            end
        end
    end
end

% left subclavian artery 

% centre3 = [128 420 170];
% 
% R = 130;
% 
% r = 12;
% 
% for i=1:256
%     for j=1:440
%         for k=160:400
%             if (((sqrt((k-centre3(3))^2 + (j-centre3(2))^2) - R)^2 + (i-centre3(1))^2) <=   r^2) 
%                 V(i,j,k) = 1.;
%             end
%         end
%     end
% end


centre1 = [256 300 180];

R = 110;

r = 12;

for i=1:256
    for j=1:440
        for k=160:400
            if (((sqrt((k-centre1(3))^2 + (i-centre1(1))^2) - R)^2 + (j-centre1(2))^2) <=   r^2) 
                V(i,j,k) = 1.;
            end
        end
    end
end

%vtkwrite('CT_3D.vtk', 'structured_points', 'indicator', (((V(1:330,150:end,k-329:k)))));


% for CT_3D_5.vtk
%V =  smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3...
%    (smooth3(permute(V,[2 1 3])))))))))));


V =  smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3...
    (smooth3(permute(V,[1 2 3])))))))))));


%original bounds
%vtkwrite('CT_3D_4.vtk', 'structured_points', 'indicator', (((V(150:end,1:330,k-329:k)))));
%vtkwrite('CT_3D_5.vtk', 'structured_points', 'indicator', (((V(100:end,1:330,k-329:k)))));
vtkwrite('CT_3D_6.vtk', 'structured_points', 'indicator', (((V))));
%vtkwrite('CT_3D_lowRes.vtk', 'structured_points', 'indicator', (((V(150:4:end,1:4:330,k-329:4:k)))));


% [m,n,k] = size(V);
% x = (150:m)/(m-150+1);
% y = (1:330)/331;
% z = (k-328:k)/329;
% [X,Y,Z] = meshgrid(x,y,z);
% X =  permute(X,[2 1 3]);
% Y =  permute(Y,[2 1 3]);
% Z =  permute(Z,[2 1 3]);

%vtkwrite('CT_3D_3.vtk', 'structured_grid',X,Y,Z, 'scalars', 'indicator', (((V(150:end,1:330,k-328:k)))));



 

