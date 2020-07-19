% CLAB3 


%   
I = imread('stereo2012a.jpg');

% imshow(I);
% 
% display('click mouse for 12 features...');
% uv = ginput(12);    % Graphical user interface to get 12 points
% display(uv);
% 
% % Pre-defined XYZ coordinates 
% XYZ=[0,0,0;
%     0,42,21;
%     0,42,7;
%     7,35,0;
%     21,35,0;
%     0,28,21;
%     0,28,7;
%     7,21,0;
%     21,21,0;
%     14,0,7;
%     14,0,21;
%     21,0,21];
% 
% % save uv coordinates and XYZ coordinates
% save data_a uv XYZ

% load uv and XYZ
load data_a

C= calibrate(I,XYZ,uv);
[K, R, t]=vgg_KR_from_P(C,0);
K=K/K(3,3);
theta_y=atan2(-R(3,1),sqrt(R(3,2)^2+R(3,3)^2));
theta_y=(theta_y/3.14195)*180;

%display C, K, R, t, theta_y
disp(C);
disp(K);
disp(R);
disp(t);
disp(theta_y);




function P= calibrate(im,XYZ,uv)
    homo=ones([12,1]);
    % convert to homogenous system
    uv_homo=horzcat(uv,homo);
    XYZ_homo=double(horzcat(XYZ,homo));
    % construct matrix A
    A=zeros([2*12,12]);
    for i=1:12
        l=i*2;
        A(l-1,:)=horzcat([0 0 0 0],horzcat((-uv_homo(i,3))*XYZ_homo(i,:),uv_homo(i,2)*XYZ_homo(i,:)));
        A(l,:)=horzcat(uv_homo(i,3)*XYZ_homo(i,:),horzcat([0 0 0 0], (-uv_homo(i,1))*XYZ_homo(i,:))); 
    end
    
    % computer P
    [U,S,V] = svd(A);
    p=V(:,12);
    p=p/norm(p);
    P=reshape(p,[4,3])';
    
    % plot the uv coordinates
    imshow(im);
    hold on;
    plot(uv(:,1),uv(:,2),'b*','LineWidth', 2, 'MarkerSize', 10);
    
    % project XYZ coordinates back into images.
    MSE=0;
    uv_check=(P*(XYZ_homo'))';
    for j=1:12
     uv_check_norm(j,:)=[uv_check(j,1)/uv_check(j,3) uv_check(j,2)/uv_check(j,3)];
     MSE=MSE+ norm(uv_check_norm(j,:)-uv(j,:));
    end
    MSE=MSE/12;
    % print mean squared error
    disp(MSE);
    hold on;
    plot(uv_check_norm(:,1),uv_check_norm(:,2),'r+','LineWidth', 2, 'MarkerSize', 10);
    
    
    % find the vinishing point and draw lines
    % dash line: parallel lines
    % solid line: line from origin to vanishing point

    % X axis as magenta line
    find_vanishing_point(uv,uv_homo,4,5,8,9,'m--','m-');
    % Y axis as yellow line 
    find_vanishing_point(uv,uv_homo,2,6,4,8,'y--','y-');
    % Z axis as white line
    find_vanishing_point(uv,uv_homo,2,3,6,7,'w--','w-');
end
function find_vanishing_point(uv,uv_homo,p1,p2,p3,p4,c1,c2)
    % 2 line equations 
    s1=cross(uv_homo(p1,:),uv_homo(p2,:));
    s2=cross(uv_homo(p3,:),uv_homo(p4,:));
    
    % the vanishing point
    C=cross(s1,s2);
    % normalised vanishing point
    C_norm=[C(1)/C(3),C(2)/C(3)];
    
    % plot dash lines(2 line segments) and solid line(origin to vanishing point)
    hold on;
    plot([uv(p1,1) C_norm(1)],[uv(p1,2) C_norm(2)], c1,'LineWidth', 2, 'MarkerSize', 10);
    plot([uv(p3,1) C_norm(1)],[uv(p3,2) C_norm(2)], c1,'LineWidth', 2, 'MarkerSize', 10);
    plot([uv(1,1) C_norm(1)],[uv(1,2) C_norm(2)], c2,'LineWidth', 2, 'MarkerSize', 10);
end

