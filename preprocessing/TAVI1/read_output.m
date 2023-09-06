close all
clear all

addpath('../');
show_3D=0;

path = './';
% was 1_1 for the first simulation
name = 'Segmentation.nii';
filename = strcat(path,name)
V = niftiread(filename);
refImage = V(:,:,:);

[m,n,k] = size(refImage);

% perfrom upsampling

MU = 2; % upsampling factor

B = zeros([size(V,1)*MU size(V,2)*MU size(V,3)*MU]);
B(1:MU:end,1:MU:end,1:MU:end) = V;

sz = MU^3;
H = fspecial3('average',[sz sz sz]);
C = convn(B,H,'same');
[m,n,k] = size(C);


% be careful of the 0.05  residual for the smoothing  (it is usually around max(C)/2.)
for ii =1:m
    for jj=1:n
        for kk=1:k
            if (C(ii,jj,kk) >= (0.05))
                C(ii,jj,kk) =1;
            else 
                C(ii,jj,kk) =0;
            end
                
        end
    end
end


if show_3D

    

    
x = (1:m);
y = (1:n);
z = (155:490);
%z = (1:k);

[X,Y,Z] = meshgrid(y,x,z);



isosurface(X,Y,Z, ((C(:,:,155:490))),0.99)
%isosurface(X,Y,Z, ((C(:,:,1:k))),0.99)

axis equal
colormap hot
h = get(gca,'DataAspectRatio') 
% if h(3)==1
%       set(gca,'DataAspectRatio',[1 1 1*max(h(1:2))])
% else
%       set(gca,'DataAspectRatio',[1 1 h(3)])
% end

else



x = (1:376);
y = (24:488);

%vtkwrite('CT_3D.vtk', 'structured_points', 'indicator', (((V(1:330,150:end,k-329:k)))));


C =  smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3...
    (smooth3(permute(C,[1 2 3])))))))))));

C = flip(C, 2);
C = flip(C, 1);

%original bounds
%vtkwrite('CT_3D_4.vtk', 'structured_points', 'indicator', (((V(150:end,1:330,k-329:k)))));
vtkwrite('CT_3D_TAVI001_3.vtk', 'structured_points', 'indicator', (((C(:,:,155:490)))));
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
end


 

