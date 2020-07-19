% Q4
% 1
% load image
face_image=imread("me2.jpg");
% crop image
face_image_cropped=imcrop(face_image,[350,650,100+256,100+256]);

% convert to grayscale image
face_iamge_cropped_gray=rgb2gray(face_image_cropped);
subplot(1,2,1),imshow(face_image_cropped);
subplot(1,2,2),imshow(face_iamge_cropped_gray);


%2
% add Gaussian Noise
image_gaussian=imnoise(face_iamge_cropped_gray,'gaussian',0,(15/256)^2);
% filter the iamge range within 0-255.
image_gaussian(image_gaussian>255)=255;
image_gaussian(image_gaussian<0)=0;

% display images(uncomment to see)
% subplot(1,2,1),imshow(face_iamge_cropped_gray);
% subplot(1,2,2),imshow(image_gaussian);


% display histograms(uncomment to see)
% subplot(1,2,1),imhist(face_iamge_cropped_gray);
% subplot(1,2,2),imhist(image_gaussian);


% 4
% load Gaussian kernel and return the filtered image
kernel_1 = fspecial('gaussian',[5 5],0.1);
filtered_image_1=my_Gauss_filter(image_gaussian,kernel_1);
kernel_2 = fspecial('gaussian',[5 5],0.5);
filtered_image_2=my_Gauss_filter(image_gaussian,kernel_2);
kernel_3 = fspecial('gaussian',[5 5],1);
filtered_image_3=my_Gauss_filter(image_gaussian,kernel_3);
kernel_4 = fspecial('gaussian',[5 5],5);
filtered_image_4=my_Gauss_filter(image_gaussian,kernel_4);
kernel_5 = fspecial('gaussian',[5 5],10);
filtered_image_5=my_Gauss_filter(image_gaussian,kernel_5);

% print and compare origin image and filtered images.(uncomment to see)
% subplot(2,3,1),imshow(image_gaussian);
% subplot(2,3,2),imshow(filtered_image_1);
% subplot(2,3,3),imshow(filtered_image_2);
% subplot(2,3,4),imshow(filtered_image_3);
% subplot(2,3,5),imshow(filtered_image_4);
% subplot(2,3,6),imshow(filtered_image_5);


%5
% compare with build-in filter. (uncomment to see)
% subplot(3,2,1),imshow(filtered_image_1);
% subplot(3,2,2),imshow(imfilter(image_gaussian,kernel_1));
% subplot(3,2,3),imshow(filtered_image_3);
% subplot(3,2,4),imshow(imfilter(image_gaussian,kernel_3));
% subplot(3,2,5),imshow(filtered_image_5);
% subplot(3,2,6),imshow(imfilter(image_gaussian,kernel_5));

function output_image=my_Gauss_filter(noisy_image,kernel)
    [row,col]=size(noisy_image);
    [ker_row,ker_col]=size(kernel);
    % init output image
    output_image=zeros(row,col);
    mid_row=(ker_row+1)/2;
    mid_col=(ker_col+1)/2;
    for i=1:row-ker_row
        for j=1:col-ker_col
            % sliding sub-matrix windows from noisy image
            tmp=noisy_image(i:i+ker_row-1,j:j+ker_col-1);
            % element-wise multiple sliced windows and kernel then sum up.
            tmp=sum(sum(double(tmp).*kernel));
            % add gussian noise to the midium cell.
            output_image(i+mid_row,j+mid_col)=tmp;
        end
    end
    output_image=uint8(output_image);
end

