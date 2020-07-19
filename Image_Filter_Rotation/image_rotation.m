
% NOTE this script require a lot computation thus is a bit slow.
origin_image=imread('me3.jpg');
% resize image
origin_image=imresize(origin_image,[512 512]);
% default interpolation method
method='linear';
% method='nearest';
% method='natural';

% Q1 and Q2
% rotate with forward, backward and matlab build-in rotation funtion afor comparsion.

% -90
rotate_negative_90_forward=my_rotation_forward(origin_image,-90);
rotate_negative_90_backward=my_rotation_backward(origin_image,-90,method);
rotate_negative_90=imrotate(origin_image,90);

% -45
rotate_negative_45_forward=my_rotation_forward(origin_image,-45);
rotate_negative_45_backward=my_rotation_backward(origin_image,-45,method);
rotate_negative_45=imrotate(origin_image,45);

%  -15
rotate_negative_15_forward=my_rotation_forward(origin_image,-15);
rotate_negative_15_backward=my_rotation_backward(origin_image,-15,method);
rotate_negative_15=imrotate(origin_image,15);

% 45
rotate_45_forward=my_rotation_forward(origin_image,45);
rotate_45_backward=my_rotation_backward(origin_image,45,'linear');
rotate_45=imrotate(origin_image,-45);

% 90
rotate_90_forward=my_rotation_forward(origin_image,90);
rotate_90_backward=my_rotation_backward(origin_image,90,'linear');
rotate_90=imrotate(origin_image,-90);

% plot all these images in order.
subplot(5,3,1),imshow(rotate_negative_90_forward);
subplot(5,3,2),imshow(rotate_negative_90_backward);
subplot(5,3,3),imshow(rotate_negative_90);

subplot(5,3,4),imshow(rotate_negative_45_forward);
subplot(5,3,5),imshow(rotate_negative_45_backward);
subplot(5,3,6),imshow(rotate_negative_45);

subplot(5,3,7),imshow(rotate_negative_15_forward);
subplot(5,3,8),imshow(rotate_negative_15_backward);
subplot(5,3,9),imshow(rotate_negative_15);

subplot(5,3,10),imshow(rotate_45_forward);
subplot(5,3,11),imshow(rotate_45_backward);
subplot(5,3,12),imshow(rotate_45);


subplot(5,3,13),imshow(rotate_90_forward);
subplot(5,3,14),imshow(rotate_90_backward);
subplot(5,3,15),imshow(rotate_90);

% Q3

% 15
rotate_linear_15=my_rotation_backward(origin_image,15,'linear');
rotate_nearest_15=my_rotation_backward(origin_image,15,'nearest');
rotate_natural_15=my_rotation_backward(origin_image,15,'natural');

% 45
rotate_linear_45=my_rotation_backward(origin_image,45,'linear');
rotate_nearest_45=my_rotation_backward(origin_image,45,'nearest');
rotate_natural_45=my_rotation_backward(origin_image,45,'natural');

% 90
rotate_linear_90=my_rotation_backward(origin_image,90,'linear');
rotate_nearest_90=my_rotation_backward(origin_image,90,'nearest');
rotate_natural_90=my_rotation_backward(origin_image,90,'natural');

% subplot side by side(uncomment to see)
% subplot(3,3,1),imshow(rotate_linear_15);
% subplot(3,3,2),imshow(rotate_nearest_15);
% subplot(3,3,3),imshow(rotate_natural_15);
% 
% subplot(3,3,4),imshow(rotate_linear_45);
% subplot(3,3,5),imshow(rotate_nearest_45);
% subplot(3,3,6),imshow(rotate_natural_45);
% 
% subplot(3,3,7),imshow(rotate_linear_90);
% subplot(3,3,8),imshow(rotate_nearest_90);
% subplot(3,3,9),imshow(rotate_natural_90);



function output_image=my_rotation_forward(image,degree)
    % matlab default positive is anti-clockwise thus we need to rotation in
    % opposite direction.
    degree=-degree;
    [row,col,channel]=size(image);
    
    % coordinate of the bottom right corner of the image.
    br_x=int16(cosd(degree)*(col-1)+sind(degree)*(row-1));
    br_y=int16(-sind(degree)*(col-1)+cosd(degree)*(row-1));
    
    % coordinate of the bottom left corner of the image.
    bl_x=int16(cosd(degree)*(col-1)+sind(degree)*0);
    bl_y=int16(-sind(degree)*(col-1)+cosd(degree)*0);
    
    % coordinate of the top right corner of the image.
    tr_x=int16(cosd(degree)*0+sind(degree)*(row-1));
    tr_y=int16(-sind(degree)*0+cosd(degree)*(row-1));
    
    % find the minimum x, minimum y the rotated image.
    x_min=min([br_x,bl_x,tr_x,0]);
    y_min=min([br_y,bl_y,tr_y,0]);
    
    % the image matrix will be extended automatically if needed. 
    output_image=zeros(1,1,channel);
    for a=1:row
        for b=1:col
            % for each pixel in origin image, find its new position in rotated image.
            t_i=int16(cosd(degree)*(b-1)+sind(degree)*(a-1));
            t_j=int16(-sind(degree)*(b-1)+cosd(degree)*(a-1));
            % move the R,G and B value to the new position.
            output_image(t_j-y_min+1,t_i-x_min+1,1)=image(a,b,1);
            output_image(t_j-y_min+1,t_i-x_min+1,2)=image(a,b,2);
            output_image(t_j-y_min+1,t_i-x_min+1,3)=image(a,b,3);
        end
    end
    % convert matrix to uint8 format.
    output_image=uint8(output_image);
end

function output_image=my_rotation_backward(image,degree,method)
    degree=-degree;
    [row,col,channel]=size(image);
    
    % coordinate of the bottom right corner of the image.
    br_x=int16(cosd(degree)*(col-1)+sind(degree)*(row-1));
    br_y=int16(-sind(degree)*(col-1)+cosd(degree)*(row-1));
    
    % coordinate of the top right corner of the image.
    bl_x=int16(cosd(degree)*(col-1)+sind(degree)*0);
    bl_y=int16(-sind(degree)*(col-1)+cosd(degree)*0);
    
    % coordinate of the top right corner of the image.
    tr_x=int16(cosd(degree)*0+sind(degree)*(row-1));
    tr_y=int16(-sind(degree)*0+cosd(degree)*(row-1));
    
    % find the range of x,y axis of rotated image.  
    x_min=min([br_x,bl_x,tr_x,0]);
    x_max=max([br_x,bl_x,tr_x,0]);
    y_min=min([br_y,bl_y,tr_y,0]);
    y_max=max([br_y,bl_y,tr_y,0]);
    n_col=x_max-x_min;
    n_row=y_max-y_min;
    
    % init rotated image.
    output_image=zeros(n_row,n_col,channel);
    
    % init interpolation input vector
    F_x=zeros(1,row*col);
    F_y=zeros(1,row*col);
    F_v_R=zeros(1,row*col);
    F_v_G=zeros(1,row*col);
    F_v_B=zeros(1,row*col);
    % iteratively feed in coordinate and intensity value to input vectors.
    index=1;
    for i=1:row
        for j=1:col
              F_x(index)=i;
              F_y(index)=j;
              F_v_R(index)=image(i,j,1);
              F_v_G(index)=image(i,j,2);
              F_v_B(index)=image(i,j,3);
              index=index+1;
        end
    end
    % define interpolation function with 3 color channels.
    F_R=scatteredInterpolant(F_y(:),F_x(:),F_v_R(:),method);
    F_G=scatteredInterpolant(F_y(:),F_x(:),F_v_G(:),method);
    F_B=scatteredInterpolant(F_y(:),F_x(:),F_v_B(:),method);
    
    % 
    for a=1:n_row
        for b=1:n_col 
            % each pixel in the rotated imagem, find its origin image
            % we rotate from center if a image.
            t_x=(cosd(-degree)*(b-1-n_col/2)+sind(-degree)*(a-1-n_row/2));
            t_y=(-sind(-degree)*(b-1-n_col/2)+cosd(-degree)*(a-1-n_row/2));
            if (t_y>=-row/2) && (t_y<row/2) && (t_x>=-col/2) && (t_x<col/2)
                % copy RGB value from origin image to rotated image if the pixel is in the origin image.
                output_image(a,b,1)=F_R(double(t_x+row/2+1),double(t_y+col/2+1));
                output_image(a,b,2)=F_G(double(t_x+row/2+1),double(t_y+col/2+1));
                output_image(a,b,3)=F_B(double(t_x+row/2+1),double(t_y+col/2+1));
            else
                % assign pixel as black hole if the pixel position is out
                % of range.
                output_image(a,b,1)=0;
                output_image(a,b,2)=0;
                output_image(a,b,3)=0;
            end
        end
    end
    % convert matrix to uint8 format.
    output_image=uint8(output_image);
end

