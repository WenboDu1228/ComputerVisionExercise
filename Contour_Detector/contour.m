

clc;clearvars;close all;
addpath('./pybsds');
GT_DIR = './contour-data/groundTruth';
IMAGE_DIR = './contour-data/images';
N_THRESHOLDS = 99;

imset = 'val';
imlist = get_imlist(imset);
output_dir = './contour-output/demo/';
% fn = @(x)compute_edges_dxdy(x); %% function for computing edges 
if ~exist(output_dir)
       mkdir(output_dir);
end

fprintf('Running detector:\n');
% detect_edges(imlist, fn, output_dir) % comment this once you implemented
%your own edge detector

% please implement your better edge detection function, and call it

% <c> Non-maximum-supression
%fn_yours = @(x)my_compute_edges_dxdy(x);

% <d> Color information detector
fn_yours = @(x)my_color_compute_edges_dxdy(x);
detect_edges(imlist, fn_yours, output_dir)
% 


%%% for evaluation of the results %%%%
load_pred_x = @(x)load_pred(output_dir, x);
load_gt_boundaries_x = @(x)load_gt_boundaries(x);
fprintf('Evaluating:\n');
[sample_results, threshold_results, overall_result] = pr_evaluation(int32(N_THRESHOLDS), imlist, load_gt_boundaries_x, load_pred_x);
file_name = sprintf('%s_out.txt',output_dir(1:end-1));
f = figure('visible','off');
display_results(file_name, sample_results, threshold_results, overall_result);
saveas(f,sprintf('%s_pr.pdf',output_dir(1:end-1)));


function imlist = get_imlist(name)
    fileID = fopen(sprintf('contour-data/%s.imlist',name),'r');
    imlist = fscanf(fileID,'%d');
end

function mag = compute_edges_dxdy(I)
    % """Returns the norm of dx and dy as the edge response function."""
    I = double(I)/255;
    dx = conv2(I, [-1, 0, 1],'same');
    dy = conv2(I, [-1; 0; 1], 'same');
    mag = (dx .^ 2 + dy .^ 2).^(1/2);
    mag = mag / max(mag,[],'all');
    mag = mag * 255;
    mag(mag<0) = 0;
    mag(mag>255) = 255;
    mag = uint8(mag);
end

function mag = my_compute_edges_dxdy(I)
    I = double(I)/255;
    
%    origin code    
%    dx = conv2(I, [-1, 0, 1],'same');
%    dy = conv2(I, [-1, 0, 1],'same');

% %   <a> reduce artifcts near boundarys.
%     dx = conv2(I, [-1, 0, 1],'valid');
%     dy = conv2(I, [-1; 0; 1], 'valid');
%     dx = padarray(dx,[0 1]);
%     dy = padarray(dy,[1 0]);

%   <b> Gussian filter and its derivative
    sigma=2.5;
    I = imgaussfilt(I,sigma);
    [dx, dy]= gradient(I);

%   <c> Non-maximum supression 
    mag = (dx .^ 2 + dy .^ 2).^(1/2);
    mag=non_maximum_supression(mag,dx,dy);
    mag = mag / max(mag,[],'all');
    
%   <d> weak edge detector, not used since performance is not good
%   mag=weak_edges_hy(mag);

    mag = mag * 255;
    mag(mag<0) = 0;
    mag(mag>255) = 255;
    mag = uint8(mag);
end

% <d> color information edge detector.
function mag = my_color_compute_edges_dxdy(I)
    % convert to LAB 
    LAB = rgb2lab(I);
    L=LAB(:,:,1);
    A=LAB(:,:,2);
    B=LAB(:,:,3);
    
    % performance edge detection on 3 channels.
    magL= my_compute_edges_dxdy(L);
    magA=my_compute_edges_dxdy(A);
    magB=my_compute_edges_dxdy(B);
    
    % perform element-wise max operations.
    mag=max(magL,magA);
    mag=max(mag,magB);
end

% double thresholding and Edge Tracking by Hysteresis
function mag=weak_edges_hy(mag)
    [nrows,ncols]=size(mag);
    
    % track the index of strong edges and weak edges
    strongRow = zeros(1,nrows*ncols);
    strongCol = zeros(1,nrows*ncols); 
    weakRow = zeros(1,nrows*ncols);  
    weakCol = zeros(1,nrows*ncols);  
    strong = 1;
    weak = 1;
    
    
    % set the threholdings. 
    high= max(max(mag))*0.25;
    low = high*0.51;
    
    % label each pixel as strong, weak or non-edge.
    for i=2:nrows-1 % row
        for j=2:ncols-1 % col
            if mag(i,j) > high
                mag(i,j) = 1;
                strongRow(strong) = i;
                strongCol(strong) = j;
                strong = strong + 1;
            elseif mag(i,j) > low 
                weakRow(weak) = i;
                weakCol(weak) = j;
                weak = weak + 1;
            else      
                mag(i,j) = 0;
            end
        end
    end
    recordMatrix=zeros(nrows,ncols);
    
    % find weak edge connected to strong edge.
    for i=1:strongIndex-1
        mag = FindConnectedWeakEdges(mag, strongRow(i),...
            strongCol(i),recordMatrix);
        break;
    end 
end



function detect_edges(imlist, fn, out_dir)
    IMAGE_DIR = './contour-data/images/';
    lis_len = length(imlist);
    tic
    for i = 1:lis_len
        imname = imlist(i);
        I = imread(sprintf('%s%s.jpg',IMAGE_DIR, string(imname)));
        %gray = rgb2gray(I);
        %mag = fn(gray);  
        mag = fn(I);
        out_file_name = sprintf('%s%s.png',out_dir, string(imname));
        imwrite(mag,out_file_name);
    end
    timeElapsed = toc;
    disp(timeElapsed);
end

function new_mag=non_maximum_supression(mag,dx,dy)
    angle = atan2(dy,dx)*180/pi;
    [nrows, ncols] = size(mag);
    new_mag=zeros(size(mag));
    for i=2:nrows-1
        for j=2:ncols-1
            % perform non maximum supression on each pixel
            new_mag(i,j)=non_maximum_supression_helpper(mag,angle,i,j,dx,dy);
        end
    end
end

function pixel_mag=non_maximum_supression_helpper(mag,angle,i,j,dx,dy)
    pixel_mag=0;
    % 4 cases of norm direction
    if (angle(i,j)>=0 &&angle(i,j)<=45 ||angle(i,j)>=-180 &&angle(i,j)<=-135)
        right_bottom=[mag(i,j+1), mag(i+1,j+1)];
        left_top=[mag(i,j-1), mag(i-1,j-1)];
    elseif (angle(i,j)>45 &&angle(i,j)<=90 ||angle(i,j)>=-135 &&angle(i,j)<=-90)
        right_bottom=[mag(i+1,j), mag(i+1,j+1)];
        left_top=[mag(i-1,j-1), mag(i-1,j)];
    elseif (angle(i,j)>=90 &&angle(i,j)<=135 ||angle(i,j)>=-90 &&angle(i,j)<=-45)
        right_bottom=[mag(i+1,j), mag(i+1,j-1)];
        left_top=[mag(i-1,j), mag(i-1,j+1)];
    elseif (angle(i,j)>=135 &&angle(i,j)<=180 ||angle(i,j)>=-45 &&angle(i,j)<=0)
        right_bottom=[mag(i,j-1), mag(i+1,j-1)];
        left_top=[mag(i,j+1), mag(i-1,j+1)];
    end
    % test if we elinimate the pixel because it is not local-maximum
    delta=abs(dy(i,j)/mag(i,j));
    if (mag(i,j) >= ((left_top(2)-left_top(1))*delta+left_top(1))) && ...
    mag(i,j) >= ((right_bottom(2)-right_bottom(1))*delta+right_bottom(1))
        pixel_mag= mag(i,j);
    end
end

% recursive find week edges
% reference: http://justin-liang.com/tutorials/canny/#grayscale
function[Gmag] = FindConnectedWeakEdges(Gmag, row, col,recordMatrix)
    for i = -3:1:3
        for j = -3:1:3
            if (i==0) && (j==0)
                continue;
            end
            if recordMatrix(row,col)==0
                if (row+i > 0) && (col+j > 0) && (row+i < size(Gmag,1)) && ...
                        (col+j < size(Gmag,2))
                    if (Gmag(row+i,col+j) > 0) && (Gmag(row+i,col+j) < 1)
                        Gmag(row+i,col+j) = 1;
                        recordMatrix(row,col)=1;
                        Gmag = FindConnectedWeakEdges(Gmag, row+i, col+j,recordMatrix);
                    end
                end
            end
        end
    end
end

function boundary = load_gt_boundaries(imname)
    GT_DIR = './contour-data/groundTruth/';
    gt_path = sprintf('%s%s.mat',GT_DIR, string(imname));
    bd = bsds_dataset;
    boundary = bd.load_boundaries(gt_path);
end


function img = load_pred(output_dir, imname)
    pred_path = sprintf('%s%s.png',output_dir,string(imname));
    img = double(imread(pred_path))/255.0;
end

function display_results(f, im_results, threshold_results, overall_result)
    % out_keys = ["threshold", "f1", "best_f1", "area_pr"];
    % out_name = ["threshold", "overall max F1 score", "average max F1 score", "area_pr"];
    
    overall_result
    fileID = fopen(f,'w');
    fprintf(fileID,'%s %10.6f\n','threshold',overall_result.threshold);
    fprintf(fileID,'%s %10.6f\n','overall max F1 score',overall_result.f1);
    fprintf(fileID,'%s %10.6f\n','average max F1 score',overall_result.best_f1);
    fprintf(fileID,'%s %10.6f\n','area_pr',overall_result.area_pr);
    fclose(fileID);
    
    res = reshape(struct2array(threshold_results),[],length(threshold_results))';
    recall = res(:, 2);
    precision = res(recall > 0.01, 3);
    recall = recall(recall > 0.01);
    label_str = sprintf('%0.2f, %0.2f, %0.2f', overall_result.f1, overall_result.best_f1, overall_result.area_pr);

    p = plot(recall, precision, 'Color','r');
    p.LineWidth = 2;
    grid on;
    xlim([0, 1]);
    ylim([0, 1]);
    legend([label_str]);
    xlabel('Recall');
    ylabel('Precision');
end
    