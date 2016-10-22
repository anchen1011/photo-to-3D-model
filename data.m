load('X.dat')
load('Y.dat')
load('Z.dat')
load('I.dat')
img = zeros(256,373,3);
img(:,:,1) = I(1:256,:);
img(:,:,2) = I(257:512,:);
img(:,:,3) = I(513:768,:);

figure
surf(X(2:end-1, 2:end-1),Z(2:end-1, 2:end-1),Y(2:end-1, 2:end-1)),double(img(2:end-1, 2:end-1, :))/256)
axis('equal')
shading interp