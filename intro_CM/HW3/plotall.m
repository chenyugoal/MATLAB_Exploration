%% chenyu: question d-e
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

%% implement splineImage directfigure(2);
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
colorbarly on the template
