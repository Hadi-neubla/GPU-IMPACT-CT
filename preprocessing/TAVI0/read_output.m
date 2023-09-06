close all
clear all

addpath('../');
show_3D=0;

path = './';
% was 1_1 for the first simulation
name = 'Segmentation_1_1_1_1.nii';
filename = strcat(path,name)
V = niftiread(filename);
refImage = V(:,:,:);

[m,n,k] = size(refImage);

if show_3D

x = (1:m);
y = (1:n);
z = (k-340:k);

[X,Y,Z] = meshgrid(x,y,z);

isosurface(X,Y,Z, medfilt3((refImage(:,:,k-340:k))),0.99)
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


% for CT_3D_5.vtk
%V =  smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3...
%    (smooth3(permute(V,[2 1 3])))))))))));


V =  smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3(smooth3...
    (smooth3(permute(V,[1 2 3])))))))))));


%original bounds
%vtkwrite('CT_3D_4.vtk', 'structured_points', 'indicator', (((V(150:end,1:330,k-329:k)))));
%vtkwrite('CT_3D_5.vtk', 'structured_points', 'indicator', (((V(100:end,1:330,k-329:k)))));
vtkwrite('CT_3D_6.vtk', 'structured_points', 'indicator', (((V(1:330,100:end,k-329:k)))));
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


 

