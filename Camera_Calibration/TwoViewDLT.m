% CLAB3 




% read images
left=imread('Left.jpg');
right=imread('Right.jpg');

% get the imput uv points
% imshow(left);
% uv_left = ginput(6);
% imshow(right);
% uv_right = ginput(6);
% save data_c uv_left uv_right

% load the saved uv coordinates(load data_C is for comparsion);
load data_b


% plot the selcted points
% imshow(left);
% hold on;
% plot(uv_left(:,1),uv_left(:,2),'r*','LineWidth', 2, 'MarkerSize', 10);
% imshow(right);
% hold on;
% plot(uv_right(:,1),uv_right(:,2),'r*','LineWidth', 2, 'MarkerSize', 10);


% find the matrix H

H = homography(uv_right(:,1), uv_right(:,2), uv_left(:,1), uv_left(:,2));
[nrows ,ncols ,tmp]=size(left);
disp(H);

% find the boundary of warped image
index=1;
rec=zeros([nrows*ncols,2]);
for i=1:nrows
    for j=1:ncols
        x=H*[j;i;1];
        rec(index,:)=int16([x(1)/x(3),x(2)/x(3)]);
        if x(1)/x(3)>32700
            disp([x(1)/x(3),x(2)/x(3)]);
        end
        index=index+1;
    end
end
t_min=min(rec,[],1);
t_max=max(rec,[],1);
y_min=t_min(2);
x_min=t_min(1);
y_max=t_max(2);
x_max=t_max(1);

% warp the left image. 
wrapped=zeros([int16(y_max-y_min), int16(x_max-x_min)]);
for i=1:nrows
    for j=1:ncols
        x=H*[j;i;1];
        a=int16([x(1)/x(3),x(2)/x(3)]);
        wrapped((a(2)-y_min)+1,(a(1)-x_min)+1)=left(i,j);
    end
end
imshow(uint8(wrapped));



function H = homography(u2Trans, v2Trans, uBase, vBase)
        
        % combine uv to a matrix
        uv_base=[uBase, vBase];
        uv_trans=[u2Trans, v2Trans];
        
        % convert to homogenous system
        homo=ones([6,1]);
        uv_base_homo=horzcat(uv_base,homo);
        uv_trans_homo=horzcat(uv_trans,homo);
        
        % construct matrix A
        A=zeros([12 9]);
        for i=1:6
            l=2*i;
            A(l-1,:)=horzcat([0 0 0],horzcat(-uv_trans_homo(i,3)*uv_base_homo(i,:),uv_trans_homo(i,2)*uv_base_homo(i,:)));
            A(l,:)=horzcat(uv_trans_homo(i,3)*uv_base_homo(i,:),horzcat([0 0 0],-uv_trans_homo(i,1)*uv_base_homo(i,:)));
        end
        [U,S,V] = svd(A);
        
        % find h and reshape to H
        h=V(:,end);
        h=h/norm(h);
        H=reshape(h,[3,3])';
        % this piece can be used to varification H is correct
%         t=(H*uv_base_homo')'
%         for i=1:6
%             new_t(i,:)=[t(i,1)/t(i,3),t(i,2)/t(i,3)] 
%         end
%         disp(new_t);
%         disp(uv_trans_homo);
end




