%Q5
% load origin image
origin_image=imread('me1.jpg');

% convert to grayscale
origin_image=rgb2gray(origin_image);

% filter with my sobel fucntion and offical sobel function 
my_sobel_image=my_sobel_filter(origin_image);
sobel_image=edge(origin_image,'sobel');

% compare
subplot(1,2,1),imshow(my_sobel_image);
subplot(1,2,2),imshow(sobel_image);

function output_image=my_sobel_filter(original_image)
    %  X axis filter 
    X=double([-1,0,1;-2,0,2;-1,0,1]);
    %  Y axis filter
    Y=double([1,2,1;0,0,0;-1,-2,-1]);
    
    %  horizontal derivative approximations 
    G_x=my_filter(original_image,X);
    %  vertical derivative approximations 
    G_y=my_filter(original_image,Y);
    
    % combine two gradient approximations 
    output_image=uint8(sqrt(double(G_x).^2+double(G_y).^2));
end

function output_image=my_filter(input_image,kernel)
    [row,col]=size(input_image);
    [ker_row,ker_col]=size(kernel);
    
    % init output image
    output_image=zeros(row,col);
    mid_row=(ker_row+1)/2;
    mid_col=(ker_col+1)/2;
    
    for i=1:row-ker_row
        for j=1:col-ker_col
            % sliding sub-matrix windows from noisy image
            tmp=input_image(i:i+ker_row-1,j:j+ker_col-1);
            % element-wise multiple sliced windows and kernel then sum up.
            tmp=sum(sum(double(tmp).*kernel));
            % add gussian noise to the midium cell.
            output_image(i+mid_row,j+mid_col)=tmp;
        end
    end
    output_image=uint8(output_image);
end