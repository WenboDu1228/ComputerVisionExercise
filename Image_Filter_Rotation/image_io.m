%Q3 
%1
me1=imread('me1.jpg');
me2=imread('me2.jpg');
me3=imread('me3.jpg');
% subplot(1,3,1),imshow(me1);
% subplot(1,3,2),imshow(me2);
% subplot(1,3,3),imshow(me3);

%2
% resize the image
face01u6361796=imresize(me1,[720 1024]);
face02u6361796=imresize(me2,[720 1024]);
face03u6361796=imresize(me3,[720 1024]);
imwrite(face01u6361796,"face01u6361796.jpg");
imwrite(face02u6361796,"face02u6361796.jpg");
imwrite(face03u6361796,"face03u6361796.jpg");

subplot(1,3,1),imshow(face01u6361796);
subplot(1,3,2),imshow(face02u6361796);
subplot(1,3,3),imshow(face03u6361796);

%3 
% a
% I choosed 3rd image
the_image=imread('me3.jpg');
% resize the image
the_image=imresize(the_image,[512 768]);
%imshow(the_image);
% b
% display each of the three channel grayscale images separately.
the_image_R=the_image(:,:,1);
the_image_G=the_image(:,:,2);
the_image_B=the_image(:,:,1);
% display the three channel images(uncomment to see)
% subplot(1,3,1),imshow(the_image_R);
% subplot(1,3,2),imshow(the_image_G);
% subplot(1,3,3),imshow(the_image_B);

% c
% display the histgrams.(uncomment to see)
% subplot(1,3,1),imhist(the_image_R);
% subplot(1,3,2),imhist(the_image_G);
% subplot(1,3,3),imhist(the_image_B);

% d
% apply hist equalisation
the_image_eq=histeq(the_image);
the_image_R_eq=histeq(the_image_R);
the_image_G_eq=histeq(the_image_G);
the_image_B_eq=histeq(the_image_B);

% display equalized images(uncomment to see)
% subplot(2,2,1),imshow(the_image_eq);
% subplot(2,2,2),imshow(the_image_R_eq);
% subplot(2,2,3),imshow(the_image_G_eq);
% subplot(2,2,4),imshow(the_image_B_eq);
% subplot(2,2,1),imhist(the_image_eq);

% display histogram of equalized images(uncomment to see)
% subplot(2,2,2),imhist(the_image_R_eq);
% subplot(2,2,3),imhist(the_image_G_eq);
% subplot(2,2,4),imhist(the_image_B_eq);

