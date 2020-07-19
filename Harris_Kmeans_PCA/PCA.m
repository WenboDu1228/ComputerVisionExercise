% load images
train_path = 'Yale-FaceA/trainingset';
test_path = 'Yale-FaceA/testset';
test = dir(fullfile(test_path,'*.png'));
train = dir(fullfile(train_path,'*.png')); 

% top k eigenface
top_k=10;

% without my face
%train_sample_size=135;
% with my face
train_sample_size=144;

% without my face
% test_sample_size=10;
% with my face
test_sample_size=11;

% find the eigenface of the image dataset
% return the mean_face, eigen vectors(U), and A.
[U,mean,mean_face,B,A]=eigenface(top_k,train_sample_size,train_path,train);

% find the dimension of image
[m,n]=size(mean_face);
% represent each image as a single data point in a high dimensional space

% show mean face
% imshow(uint8(mean_face));

% show time k eigen face.
for i=1:top_k
    x=U(:,i)+mean;
    res=reshape(x,[],n);
%     subplot(4,3,i);
%     imshow(uint8(res));
end

k_nearest=3;

% compute the projection of test faces 
[test_A,test_B]=test_PCA(mean,test_path,test,test_sample_size,U);

% find the k neighbor of test images
 for i=1:test_sample_size
        dists=zeros(0,0);
        % find the vector have smallest distance with test image.
        for j =1:train_sample_size
            dists(j)=norm(test_B(:,i)-B(:,j));
        end
        [~,idx]=sort(dists);
        idx=idx(1:k_nearest);
        % show origin image and top k neighbors side by side 
        subplot(test_sample_size,k_nearest+1,1+(i-1)*(k_nearest+1));
        a=reshape(test_A(:,i)+mean,[],n);
        imshow(uint8(a));
        for r=1: k_nearest
            a=reshape(A(:,idx(r))+mean,[],n);
            subplot(test_sample_size,k_nearest+1,1+(i-1)*(k_nearest+1)+r);
            imshow(uint8(a));
        end
 end
    
function [test_A,test_B]=test_PCA(mean,test_path,test,test_sample_size,U)
    % read test dataset as a big matrix A
    test_A=zeros(0,0);
    for i=1:test_sample_size
        F = fullfile(test_path,test(i).name);
        F=imread(F);
        I = reshape(F,[],1);
        test_A(:,i)=I;
    end
    test_A=double(test_A)-mean;
    % compute the project onto U of A
    test_B=transpose(U)*test_A;
end




% function to find the eigenface
function [U,mean,mean_face,B,A]=eigenface(top_k,sample_size,train_path,train)
    mean_face=0.0
    A=zeros(0,0);
    % read images to a big matrix A
    for k = 1:sample_size
        F = fullfile(train_path,train(k).name);
        I = reshape(imread(F),[],1);
        mean_face=mean_face+uint16(imread(F));
        A(:,k)=I;
    end
    % substract mean from A
    mean_face=double(mean_face)/sample_size;
    mean=reshape(mean_face,[],1);
    A=double(A)-mean;
    C=(transpose(A)*A);
    % get the sorted top k eigen vectors.
    [V,D]=eigs(C,top_k);
    U=A*V;
    B=transpose(U)*A;
end



