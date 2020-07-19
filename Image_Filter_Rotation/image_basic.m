% 1
% load image
A=imread('Lenna.png');
% convert image to grayscale.
A1 = rgb2gray(A);
% map image to its negative image
A2 = imcomplement(A1);
subplot(3,2,1),imshow(A1);
subplot(3,2,2),imshow(A2);

% 2
% Flip the image vertically.
A3 = flipud(A2);
subplot(3,2,3),imshow(A3);

% 3
% load colour image
A4=A;
%$ swap the red and blue colour channels.
tmp= A4(:,:,1);
A4(:,:,1)=A4(:,:,3);
A4(:,:,3)=tmp;
subplot(3,2,4),imshow(A4);


% 4
% average the input image with its vertically flipped images.
A5=(A+A3)/2;
subplot(3,2,5),imshow(A5);

% 5
A6=A1;
% read shape of image
[rows, columns] = size(A6);
% create a matrix with random number between 0-255
x1 = uint8(randi([0, 255], [rows,columns]));
% add random value to image
A6(:,:)=A6(:,:)+x1;
% clip the new image
A6(A6>255)=255;
subplot(3,2,6),imshow(A6);
