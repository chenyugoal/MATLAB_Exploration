%% Generate the image

% The location of each pixel
nX = 200;
nY = 200;
% Create a blank image
I = zeros(nX, nY);
cx = 75;
cy = 75;
rx = 30;
ry = 45;
for i = 1:nX
    for j = 1:nY
        if ((i < (cx +rx) && i > (cx - rx)) && (j<(cy+ry) && j > (cy-ry)))
            I(j,i) =100;
        end
    end
end

%Display
figure;
imagesc(I);
axis image;
title('I(x)');
set(gca,'ydir', 'normal'); % Put origin at bottom left

%% Record landmarks and Find scale factor

X = [50 100 75 50 100; 35 35 75 115 115];
Y = [75 150 115 75 150;50 50 115 175 175];


syms s
part = (X.*s-Y).^2
E = sum(part(:))
scale = solve(diff(E)==0)

%% Transform the landmarks X

sX = X.*scale;
figure;
scatter(sX(1,:),sX(2,:),'b');
hold on;
scatter(X(1,:),X(2,:),'c');
scatter(Y(1,:),Y(2,:),'r');
grid on;

%% Transform the image "naively"

% initialize an image of all zeros
ITransformed = zeros(size(I));
for i = 1 : nY % loop through each row
    for j = 1 : nX % loop through each column
        % we are looking for the value to assign to Ibx(j,i)
        % find the position to look at in the image J
        iLook = i/scale;% you should implement this line
        jLook = j/scale;% you should implement this line
        % round them to the nearest integer
        iLookRound = round(iLook);
        jLookRound = round(jLook);
        % check if we're out of bounds,
        if iLookRound < 1 || iLookRound > nY || jLookRound < 1 || jLookRound > nX
        % if so, fill the image with the value zero
        ITransformed(j,i) = 0;
        else
        % otherwise, assign the value in our image at this point
        ITransformed(j,i) = I(jLookRound,iLookRound);
        end
        % don't forget to index your images by (row,column) and not (x,y) !
    end
end
%% 

%Display
figure;
imagesc(ITransformed);
axis image;
title('J(x)');
set(gca,'ydir', 'normal'); % Put origin at bottom left
grid on
