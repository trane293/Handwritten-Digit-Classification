%=========================================================================%
%                 HANDWRITTEN DIGIT CLASSIFICATION USING IMAGE MOMENTS                        %
%=========================================================================%
%           Author - Anmol Sharma
%      Affiliation - DAV Institute of Engineering & Technology
%      Description - The code is used to calculate Zernike Moments that are
%                    invariant to translation, scaling and rotation of a
%                    mass in mammograms. The input images are segmented
%                    using Chan and Vese segmentation algorithm. After
%                    that, the moments are computed using a fast method
%                    called Kintner's method. Rotation invariance using
%                    Phase cancellation is utilized to preserve moment
%                    complexity. Translation invariance is achieved by
%                    taking the centroid of mass as origin. Scale
%                    invariance is achieved using normalizing the moments
%                    themselves by taking m00 as normalizing moment. 
%        Copyright - Under GNU GPL license. Please visit GNU GPL Homepage
%                    to read and accept the terms. In short, personal use,
%                    distribution, modification, editing of any sort is
%                    permitted as long as original credits are included in
%                    the file. You may not remove the original authors
%                    name. 
%        Last Edit - 12/07/2014 1.15 am
%=========================================================================%
%% Variables Used

%       database_image = Cell array 1x57 containing 57 mass images obtained
%                        after running Chan and Vese segmentation algorithm
%  database_edge_image = Cell array 1x57 containing margin images of the
%                        mass image
%                thres = Threshold value. Used to limit the value of
%                        moments. If moment value becomes lesser than thres
%                        it is regarded as zero.
%             typecomp = 1 means Kintner's Method for computation of
%                        Zernike Moments are used. 
%              typeinv = 1 means Rotation invariance by Phase Cancellation
%                        is utilized as proposed by Flusser et al.
%           orderShape = Order upto which moments from shape image are
%                        calcualated. (5 for lower order, 10 for higher
%                        order).
%          orderMargin = Order upto which moments from margin image are
%                        calculated. (5 for lower order,10 for higher order). 
%% Read the images. 
clear all
close all
tic
currentFolder = pwd;
t = 1;
for k=100:115
    try
    filename = strcat(num2str(k),'.png');
        % Using the original masks output by chan vese segmentation algorithm.
    database_image{t} = not(logical(rgb2gray(imread([currentFolder, '\Images\', filename]))));
    database_edge_image{t} = edge(database_image{t}, 'sobel');
    t = t+1;
    catch
    end
end
for k=10:95
    try
    filename = strcat(num2str(k),'.png');
        % Using the original masks output by chan vese segmentation algorithm.
    database_image{t} = not(logical(rgb2gray(imread([currentFolder, '\Images\', filename]))));
    database_edge_image{t} = edge(database_image{t}, 'sobel');
    t = t+1;
    catch
    end
end
j = 1;
labels = zeros(50, 10);
for i = 1:50
    labels(i, j) = 1;
    if (mod(i, 5) == 0)
        j = j + 1;
    end
end
%% For debugging
% for i = 1:50
%     imshow(database_image{i});
%     disp(i);
%     pause;
% end
%% Find margin images. 
% for k=1:3
%     database_edge_image{k} = edge(database_image{k}, 'sobel');
% end
%% Do not touch these
thres = 1e-3;
% Best results with 1, 0
typecomp = 0; %0 Use Kintner's method 1 Using Geomteric central moments
typeinv = 0; % 1 Use rotation invariance by phase cancellation 0 Rotation normalization

%  
 
% database_image{4} = not(logical(rgb2gray(imread('C:\Users\Anmol\Desktop\im4.png'))));

%% Set orders of moments to be computed
orderShape = 4;
orderMargin = 4;
%% Calculate Zernike moments from shape images.
for i=1:size(database_image, 2)
    [invts{i},invmat{i},ind{i},normom{i}]=zermi(database_image{i},orderShape,thres,typecomp,typeinv);
    invariant_values_shape(i,:) = invts{1,i}(1,:);
    %invts{i} = rotmi(database_image{i},orderShape);
end

for i = 1:11 % wrap it up. Convert cell to matrix.
    max_im = max(invariant_values_shape(:, i));
    invariant_values_shape(:, i) = invariant_values_shape(:, i)./max_im;
end

%% Calculate Zernike moments from margin images.
for i=1:size(invariant_values_shape, 2)
    [invts{i},invmat{i},ind{i},normom{i}]=zermi(database_edge_image{i},orderMargin,thres,typecomp,typeinv);
    invariant_values_margin(i,:) = invts{1,i}(1,:);
    max_im = max(invariant_values_margin(:, i));
    invariant_values_margin(:, i) = invariant_values_margin(:, i)./max_im;
end


for i = 1:size(invariant_values_shape, 2)
    standardDev(1, i) = std(invariant_values_shape(:, i));
end

%% Combine moments from shape and margin (For complex feature subsets)

invariant_values = [invariant_values_shape];%, invariant_values_margin];

%% Load Target matrix for ANN based Classification (DO NOT TOUCH) (KEEP COMMENTED)
%load('C:\Users\Anmol\Desktop\Project @ NIT Silchar\Implementation\zernike_cancer_targets.mat');

%% Final householding chores
% Best results 3, 4
col1 = 7;
col2 = 8;
c = linspace(1,100,length(invariant_values(:,1)));
figure,
hold on
scatter(invariant_values(1:5, col1), invariant_values(1:5, col2),100, 'p'), 
scatter(invariant_values(6:10, col1), invariant_values(6:10, col2),200,  'o'),
scatter(invariant_values(11:15, col1), invariant_values(11:15, col2),300,  'x'),
scatter(invariant_values(16:20, col1), invariant_values(16:20, col2),300,  'd'),
scatter(invariant_values(21:25, col1), invariant_values(21:25, col2),300,  '<'),
scatter(invariant_values(26:30, col1), invariant_values(26:30, col2),300,  '>'),
scatter(invariant_values(31:35, col1), invariant_values(31:35, col2),300,  '^'),
scatter(invariant_values(36:40, col1), invariant_values(36:40, col2),300,  '+'),
scatter(invariant_values(41:45, col1), invariant_values(41:45, col2),300,  's'),
scatter(invariant_values(46:50, col1), invariant_values(46:50, col2),100, '*'),
legend('0', '1', '2', '3', '4', '5', '6', '7', '8', '9');
axis([-1 1 -1 1]);

%% For 3D
% scatter(invariant_values(46:50, col1), invariant_values(46:50, col2),invariant_values(46:50, col3),300, '*'), 
% scatter(invariant_values(1:5, col1), invariant_values(1:5, col2),invariant_values(1:5, col3),300, 'p'), 
% hold
% scatter(invariant_values(6:10, col1), invariant_values(6:10, col2), invariant_values(6:10, col3),300,  'o'),
% scatter(invariant_values(11:15, col1), invariant_values(11:15, col2),invariant_values(11:15, col3),300,  'x'),
% scatter(invariant_values(16:20, col1), invariant_values(16:20, col2),invariant_values(16:20, col3), 300,  'd'),
% scatter(invariant_values(21:25, col1), invariant_values(21:25, col2),invariant_values(21:25, col3), 300,  '<'),
% scatter(invariant_values(26:30, col1), invariant_values(26:30, col2),invariant_values(26:30, col3), 300,  '>'),
% scatter(invariant_values(31:35, col1), invariant_values(31:35, col2),invariant_values(31:35, col3), 300,  '^'),
% scatter(invariant_values(36:40, col1), invariant_values(36:40, col2),invariant_values(36:40, col3), 300,  '+'),
% scatter(invariant_values(41:45, col1), invariant_values(41:45, col2),invariant_values(41:45, col3), 300,  's'),
%% 
elasped_time = toc;
warning('on', 'all'); %Turns on warnings, to prevent this script to break MATLAB
disp(['Calculation is complete now, ' 'Total time taken is --> ' num2str(elasped_time) ' seconds']);