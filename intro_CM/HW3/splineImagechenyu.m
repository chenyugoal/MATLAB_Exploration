%I is the origial image
%J is the Target image, 
%epsilon is the step, 
%nIter is the iteration times (should be large enough)

function [ID,bx,by] = splineImagechenyu(I,J,epsilon,nIter)
% initialize velocity field
bx = zeros(size(I));
by = zeros(size(I)); 
for i = 1 : nIter
    % deform image
    ID = applyVToImage(I,bx,by);
    [gradIx,gradIy] = gradient(ID);
    gradCostx = (J-ID).*gradIx;  %J is the Target image
    gradCosty = (J-ID).*gradIy;  %ID is the deformed I
    
    %gradCostx = sum((J-ID).*gradIx); wierd
    %gradCosty = sum((J-ID).*gradIy); wierd
    
    %gradCostx = sum(sum((J-ID)*gradIx)); wierd
    %gradCosty = sum(sum((J-ID)*gradIy)); wierd
    
    bx = bx - gradCostx*epsilon;
    by = by - gradCosty*epsilon;
    
    %compute cost function after each iteration
    Cost(i) = 0.5*sum(sum((ID - J).^2));

%plot cost function to check whether it's decreasing
figure(1)
plot(Cost,'m')
title('Cost Function(iteration)');
xlabel('Iteration Time'); 
ylabel('Value of Cost Function');

end


function ID = applyVToImage(I,bx,by)
% deform image by composing with x-b (interpolating at specified points)
[X,Y] = meshgrid(1:size(I,2),1:size(I,1));
ID = interp2(I,X-bx,Y-by,'linear',0);