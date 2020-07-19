% load two images
image1=imread("mandm.png");
image2=imread("peppers.png");
image2=im2uint8(image2);


%---------------different clusters number comparsion test-------------------

% image1, 3 coordinates (3,5,7) clusters.
% subplot(2,2,1);
% imshow(image1);
% subplot(2,2,2);
% imshow(image_segmentation(image1,3,3, false));
% subplot(2,2,3);
% imshow(image_segmentation(image1,3,5, false));
% subplot(2,2,4);
% imshow(image_segmentation(image1,3,7, false));

% image1, 5 coordinates (3,5,7) clusters.
% subplot(2,2,1);
% imshow(image1);
% subplot(2,2,2);
% imshow(image_segmentation(image1,5,3, false));
% subplot(2,2,3);
% imshow(image_segmentation(image1,5,5, false));
% subplot(2,2,4);
% imshow(image_segmentation(image1,5,7, false));

% image2, 3 coordinates (3,5,7) clusters.
% subplot(2,2,1);
% imshow(image2);
% subplot(2,2,2);
% imshow(image_segmentation(image2,3,3, false));
% subplot(2,2,3);
% imshow(image_segmentation(image2,3,5, false));
% subplot(2,2,4);
% imshow(image_segmentation(image2,3,7, false));

% image2, 5 coordinates (3,5,7) clusters.
% subplot(2,2,1);
% imshow(image2);
% subplot(2,2,2);
% imshow(image_segmentation(image2,5,3, false));
% subplot(2,2,3);
% imshow(image_segmentation(image2,5,5, false));
% subplot(2,2,4);
% imshow(image_segmentation(image2,5,7, false));
% ----------------------------------------------------------------------



%-----------with or without corrdinates comparsion test-------------------
% subplot(1,3,1);
% imshow(image1);
% subplot(1,3,2);
% imshow(image_segmentation(image1,3,5, false));
% subplot(1,3,3);
% imshow(image_segmentation(image1,5,5, false));

% subplot(1,3,1);
% imshow(image2);
% subplot(1,3,2);
% imshow(image_segmentation(image2,3,5, false));
% subplot(1,3,3);
% imshow(image_segmentation(image2,5,5, false));
% ----------------------------------------------------------------------

%-----------kmeans++ and random initialization test-------------------
% subplot(1,3,1);
% imshow(image1);
% subplot(1,3,2);
% imshow(image_segmentation(image1,5,3, false));
% subplot(1,3,3);
% imshow(image_segmentation(image1,5,3, true));
% 
% subplot(1,3,1);
% imshow(image1);
% subplot(1,3,2);
% imshow(image_segmentation(image1,5,5, false));
% subplot(1,3,3);
% imshow(image_segmentation(image1,5,5, true));
% 
% subplot(1,3,1);
% imshow(image1);
% subplot(1,3,2);
% imshow(image_segmentation(image1,5,7, false));
% subplot(1,3,3);
% imshow(image_segmentation(image1,5,7, true));
%---------------------------------------------------------------------


%-----------kmeans++ and random initialization test-------------------
subplot(1,3,1);
imshow(image2);
subplot(1,3,2);
imshow(image_segmentation(image2,5,3, false));
subplot(1,3,3);
imshow(image_segmentation(image2,5,3, true));
% 
% subplot(1,3,1);
% imshow(image2);
% subplot(1,3,2);
% imshow(image_segmentation(image2,5,5, false));
% subplot(1,3,3);
% imshow(image_segmentation(image2,5,5, true));
% % 
% subplot(1,3,1);
% imshow(image2);
% subplot(1,3,2);
% imshow(image_segmentation(image2,5,7, false));
% subplot(1,3,3);
% imshow(image_segmentation(image2,5,7, true));
%---------------------------------------------------------------------


function segamented_image=image_segmentation(input_image,dimension,N, is_kmeans_pp)
    [nrows,ncols,tmp]=size(input_image);
    
    % convert pixel value to high dimension points.
    points=image2points(input_image,dimension);
    
    % feed points to kmean(with or without kmeans++ init method)
    % return a vector contain mean value and center of each point. 
    [means,centers]=my_kmeans(points, N, is_kmeans_pp);
    
    % do color segamentation using result from kmeans algorithms.
    segamented_image=zeros(nrows,ncols,3);
    for y=1:nrows
        for x=1:ncols
            segamented_image(y,x,:)=means(centers((y-1)*ncols+x,:),1:3);
        end
    end
    
    % convert image back to rgb values.
    segamented_image=lab2rgb(segamented_image);
end

function points=image2points(input_image,dimension)
     % numbers of rows and columns
    [nrows,ncols,tmp]=size(input_image);
    % convert rgb to lab
    input_image=rgb2lab(input_image);
    index=1;
    % each cells in vectors is the encode of a pixel.
    points=zeros(nrows*ncols,dimension);
    for y=1:nrows
        for x=1:ncols
            % 5d encoding
            if dimension==5
                points(index,:)=[input_image(y,x,1) input_image(y,x,2) input_image(y,x,3) x y ];
            % 3d encoding
            else
                points(index,:)=[input_image(y,x,1) input_image(y,x,2) input_image(y,x,3)];
            end
            index = index + 1;
        end
    end
end

function [means,centers]=my_kmeans(points, N, is_kmeans_pp)
    [length,dimension]=size(points);
    if is_kmeans_pp
        % init center of clusters using k means ++.
        init_means=K_means_pp_init(points,N,dimension);
    else
        % init center of clusters randomly from points.
        init_means=zeros(N,dimension);
        for i=1:N
            init_means(i,:)= points(randi(length),:);
        end
    end
    means=init_means;
    last_means=init_means;
    while 1
        centers=zeros(length,1);
        % assign each pixel to its closest centers.
        for a=1:length
            pixel=points(a,:);
            min_dis=norm(pixel-means(1,:));
            min_b=1;
            for b=2:N
                dis=norm(pixel-means(b,:));
                if (dis<min_dis)
                    min_dis=dis;
                    min_b=b;
                end
            end
            % record the index of centers.
            centers(a,:)=min_b;
        end
        % re initlize mean
        means=zeros(N,dimension);
        counts=ones(N,dimension);
        
        % for all points in a cluster, computer new center.
        for k=1:length
            means(centers(k),:)=means(centers(k),:)+points(k,:);
            counts(centers(k),:)=counts(centers(k),:)+1;
        end
        means=means./counts;
        
        % converaged, escape the loop
        if norm(last_means-means)==0
            break
        end
        last_means=means;
    end
end

function means=K_means_pp_init(vectors,N,coor)
    [length,tmp]=size(vectors);
    means=zeros(N,coor);
    % uniform select the first means.
    means(1,:)=vectors(randi(length),:);
    % record the distance of each point to it closest center.
    distances=zeros(length,1);
    mean_id=2;
    % break when we find enough centers.
    while (mean_id<=N)
        % compute the distance for each point to its closest center.
        for i=1:length
            min_dis=norm(vectors(i,:)-means(1,:));
            for j=2:N
                if ~(means(N,:)==zeros(1,coor))
                    dis=norm(vectors(i,:)-means(N,:));
                    if (dis<min_dis)
                        min_dis=dis;
                    end
                end
            end
            distances(i,:)=min_dis*min_dis;
        end
        % choose a center based on the weighted disturbution.
        % proportional to squared distance.
        means(mean_id,:)=vectors(randsample(1:length,1,true,distances),:);
        mean_id=mean_id+1;
    end
end



