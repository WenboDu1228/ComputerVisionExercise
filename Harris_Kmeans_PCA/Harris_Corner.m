% load images
image_1=rgb2gray(imread("Harris_1.jpg"));
image_2=imread("Harris_2.jpg");
image_3=rgb2gray(imread("Harris_3.jpg"));
image_4=rgb2gray(imread("Harris_4.jpg"));

% parameters
sigma = 2; 
thresh = 0.01;
% (adjust in different image)
% threshold = 100000000;

% image 1
[i,j]=Harris_Cornerness(image_1,thresh,sigma,50000000);
default_corner=corner(image_1);

% matlab build-in corner detector
subplot(2,4,1),imshow(image_1);
hold on
plot(default_corner(:,1),default_corner(:,2),'g*');

% my corner detector
subplot(2,4,5),imshow(image_1);
hold on
plot(j,i,'g*');

% image 2
[i,j]=Harris_Cornerness(image_2,thresh,sigma,100000000);
default_corner=corner(image_2);

subplot(2,4,2),imshow(image_2);
hold on
plot(default_corner(:,1),default_corner(:,2),'g*');

subplot(2,4,6),imshow(image_2);
hold on
plot(j,i,'g*');

% image 3
[i,j]=Harris_Cornerness(image_3,thresh,sigma,150000000);
default_corner=corner(image_3);

subplot(2,4,3),imshow(image_3);
hold on
plot(default_corner(:,1),default_corner(:,2),'g*');

subplot(2,4,7),imshow(image_3);
hold on
plot(j,i,'g*');

% image 4
[i,j]=Harris_Cornerness(image_4,thresh,sigma,30000000);
default_corner=corner(image_4);

subplot(2,4,4),imshow(image_4);
hold on
plot(default_corner(:,1),default_corner(:,2),'g*');

subplot(2,4,8),imshow(image_4);
hold on
plot(j,i,'g*');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task: Compute the Harris Cornerness                   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [i,j]=Harris_Cornerness(input_image,thresh,sigma,threshold)
    [nrow,ncol] = size(input_image);
    % derivative masks
    dx = [-1 0 1;-1 0 1; -1 0 1];
    dy = dx';
    % compute x and y derivative of image.
    Ix = conv2(input_image,dx,'same');
    Iy = conv2(input_image,dy,'same');
    % Gaussian windows function
    g = fspecial('gaussian',max(1,fix(3*sigma)*2+1),sigma);
    % compute the sum of the product of the derivative at each pixel
    % x and x
    Sx2 = conv2(Ix.^2,g,'same');
    % y and y
    Sy2 = conv2(Iy.^2,g,'same');
    % x and y
    Sxy = conv2(Ix.*Iy,g,'same');
    % init matrix of cornerness value of each pixel.
    R=zeros(nrow,ncol);
    for i=1:nrow
        for j=1:ncol
            H=[Sx2(i,j) Sxy(i,j); Sxy(i,j) Sy2(i,j)];
            % compute cornerness of the pixel cell.
            R(i,j) = det(H)-thresh*(trace(H)*trace(H));
        end
    end
    % perform non-maximum supression and return the index
    [i,j]=NMS(R,5,threshold);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task: Perform non-maximum suppression and             %
%       thresholding, return the N corner points        %
%       as an Nx2 matrix of x and y coordinates         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [i,j] = NMS(R,window_size,threshold)
    [nrow,ncol] = size(R);
    mx=zeros([nrow,ncol]);
    % zero padding matrix with window size
    R_l=zeros([nrow+2*window_size,ncol+2*window_size]);
    R_l(window_size:nrow+window_size-1,window_size:ncol+window_size-1)=R;
    for i=1+window_size:nrow+window_size
        for j=1+window_size:ncol+window_size
            % assign all cells in the windows to the maximum value in the
            % windows.
            tmp=max(max(R_l(i-window_size:i+window_size,j-window_size:j+window_size)));
            mx(i-window_size,j-window_size)=tmp;
        end
    end
    % find all the index with local maximum and the local maximum is
    % greater then threshold
    index = find((R == mx) &(R > threshold));
    % convert index format to [rows,cols]
    [i,j]=ind2sub(size(R),index);
end