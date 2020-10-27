%% chenyu: question d-e

%set inputs
I = double(imread('0001_CC_Con.png') > 0);
J = double(imread('0003_CC_Alz.png') > 0);

epsilon = 0.05; 
nIter = 4000;

%start computing!
[IDchenyu,bxchenyu,bychenyu] = splineImagechenyu(I,J,epsilon,nIter)
initcost = 0.5*sum(sum((IDchenyu - J).^2))

%visualization
figure(2);
imagesc(I);
title('Original(x)');

figure(3);
imagesc(J);
title('Target(x)');

figure(4);
imagesc(IDchenyu);
title('Deformed(x)');

figure(5);
imagesc(I-J);
title('the template image minus the target image')
colorbar

figure(6);
imagesc(IDchenyu-J);
title('the translated template image minus the target image')
colorbar
%% spline: question f

%set inputs
I = IDchenyu; %use my deformed image instead
J = double(imread('0003_CC_Alz.png') > 0);

epsilon = 0.005; %make it much smaller when using slineImage
nIter = 4000;
sigma = 0.01;
alpha = 20;

%start computing!
[IDsp,vxsp,vysp] = splineImage(I,J,alpha,sigma,epsilon,nIter)
finalcost = 0.5*sum(sum((IDsp - J).^2))

%visualization
figure(2);
imagesc(IDchenyu);
title('First Time Deformed(x)');

figure(3);
imagesc(J);
title('Target(x)');

figure(4);
imagesc(IDsp);
title('Second Time Deformed(x)');

figure(5);
imagesc(IDchenyu-J);
title('the original image minus the target image')
colorbar

figure(6);
imagesc(IDsp-J);
title('the deformed image minus the target image')
colorbar
%% implement splineImage directly on the template
%set inputs
I = double(imread('0001_CC_Con.png') > 0);
J = double(imread('0003_CC_Alz.png') > 0);

epsilon = 0.005; 
nIter = 4000;
sigma = 0.01;
alpha = 20;

%start computing!
[ID,vx,vy] = splineImage(I,J,alpha,sigma,epsilon,nIter)

%visualization
figure(2);
imagesc(I);
title('Original(x)');

figure(3);
imagesc(J);
title('Target(x)');

figure(4);
imagesc(ID);
title('Deformed(x)');

figure(5);
imagesc(I-J);
title('the template image minus the target image')
colorbar

figure(6);
imagesc(ID-J);
title('the translated template image minus the target image')
colorbar

%% Jacobian Calculations
[gradvx_x,gradvx_y] = gradient(vx);
[gradvy_x,gradvy_y] = gradient(vy);
deter = gradvx_x + gradvy_y + gradvx_x.*gradvy_y - gradvy_x .* gradvx_y +1;

%figure
figure(7)
imagesc(deter)
title('the determinant of the Jacobian of the transformation')
colorbar

