%% clear the workspace and select data 
clear; clc; close all;
% set(0,'DefaultFigureVisible','off')

%% choose data
run ./Sources/CNMF_E/cnmfe_setup.m
neuron = Sources2D();

% nams = '/media/s2v2/Arash/TestData/c142.tif';
% nam = neuron.select_data(nams);

foldername = '/home/arash/Desktop';
filetype = 'tif';
files = subdir(fullfile(foldername,['*.',filetype]));
nams = files(1).name;
nam = neuron.select_data(nams);

%% parameters  
% -------------------------    COMPUTATION    -------------------------  %
pars_envs = struct('memory_size_to_use', 120, ...   % GB, memory space you allow to use in MATLAB
    'memory_size_per_patch', 0.5, ...   % GB, space for loading data within one patch
    'patch_dims', [30, 40]);  %GB, patch size

% -------------------------      SPATIAL      -------------------------  %
gSig = 0.5;       % pixel, gaussian width of a gaussian kernel for filtering the data. 0 means no filtering
gSiz = 30;        % pixel, neuron diameter
ssub = 1;         % spatial downsampling factor
with_dendrites = false;   % with dendrites or not
if with_dendrites
    % determine the search locations by dilating the current neuron shapes
    updateA_search_method = 'dilate';
    updateA_bSiz = 20;
    updateA_dist = neuron.options.dist;
else
    % determine the search locations by selecting a round area
    updateA_search_method = 'ellipse';
    updateA_dist = 5;
    updateA_bSiz = neuron.options.dist;
end
spatial_constraints = struct('connected', true, 'circular', false);  % you can include following constraints: 'circular'

% -------------------------      TEMPORAL     -------------------------  %
Fs = 30;             % frame rate
tsub = 1;           % temporal downsampling factor
deconv_options = struct('type', 'ar1', ... % model of the calcium traces. {'ar1', 'ar2'}
    'method', 'foopsi', ... % method for running deconvolution {'foopsi', 'constrained', 'thresholded'}
    'smin', -3, ...         % minimum spike size. When the value is negative, the actual threshold is abs(smin)*noise level 
    'optimize_pars', true, ...  % optimize AR coefficients
    'optimize_b', true);% optimize the baseline); 
nk = 1;             % detrending the slow fluctuation. usually 1 is fine (no detrending)
                    % when changed, try some integers smaller than total_frame/(Fs*30) 
detrend_method = 'local_min';  % compute the local minimum as an estimation of trend. 

% -------------------------     BACKGROUND    -------------------------  %
bg_model = 'svd';  % model of the background {'ring', 'svd'(default), 'nmf'}
nb = 1;             % number of background sources for each patch (only be used in SVD and NMF model)
ring_radius = 13;  % when the ring model used, it is the radius of the ring used in the background model. 
                    %otherwise, it's just the width of the overlapping area 

% -------------------------      MERGING      -------------------------  %
show_merge = false;  % if true, manually verify the merging step 
merge_thr = 0.7;     % thresholds for merging neurons; [spatial overlap ratio, temporal correlation of calcium traces, spike correlation]
method_dist = 'max';   % method for computing neuron distances {'mean', 'max'}
dmin = 10;       % minimum distances between two neurons. it is used together with merge_thr
dmin_only = 2;  % merge neurons if their distances are smaller than dmin_only. 

% -------------------------  INITIALIZATION   -------------------------  %
K = [];             % maximum number of neurons per patch. when K=[], take as many as possible 
min_corr = 0.7;     % minimum local correlation for a seeding pixel
min_pnr =10;       % minimum peak-to-noise ratio for a seeding pixel
min_pixel = 2^2;      % minimum number of nonzero pixels for each neuron
bd = 0;             % number of rows/columns to be ignored in the boundary (mainly for motion corrected data)
frame_range = [];   % when [], uses all frames 
save_initialization = false;    % save the initialization procedure as a video. 
use_parallel = true;    % use parallel computation for parallel computing 
show_init = true;   % show initialization results 
choose_params = false; % manually choose parameters 
center_psf = false;  % set the value as true when the background fluctuation is large (usually 1p data) 
                    % set the value as false when the background fluctuation is small (2p)

% ----------------------   MANUAL INTERVENTION  --------------------  %
with_manual_intervention = false; 

% -------------------------  FINAL RESULTS   -------------------------  %
save_demixed = true;    % save the demixed file or not 
kt = 1;                 % frame intervals 

% -------------------------    UPDATE ALL    -------------------------  %
neuron.updateParams('gSig', gSig, ...       % -------- spatial -------- 
    'gSiz', gSiz, ...
    'ring_radius', ring_radius, ...
    'ssub', ssub, ...
    'search_method', updateA_search_method, ...
    'bSiz', updateA_bSiz, ...
    'dist', updateA_bSiz, ...
    'spatial_constraints', spatial_constraints, ...
    'tsub', tsub, ...                       % -------- temporal -------- 
    'deconv_options', deconv_options, ...
    'nk', nk, ...
    'detrend_method', detrend_method, ...
    'background_model', bg_model, ...       % -------- background -------- 
    'nb', nb, ...
    'ring_radius', ring_radius, ...
    'merge_thr', merge_thr, ...             % -------- merging ---------
    'dmin', dmin, ...
    'method_dist', method_dist, ...
    'min_corr', min_corr, ...               % ----- initialization ----- 
    'min_pnr', min_pnr, ...
    'min_pixel', min_pixel, ...
    'bd', bd, ...
    'center_psf', center_psf);
neuron.Fs = Fs; 

%% distribute data and be ready to run source extraction 
neuron.getReady(pars_envs); 
fprintf('data distribution is done\n')

%% initialize neurons from the video data within a selected temporal range 
if choose_params
    % change parameters for optimized initialization
    [gSig, gSiz, ring_radius, min_corr, min_pnr] = neuron.set_parameters();
end

[center, Cn, PNR] = neuron.initComponents_parallel(K, frame_range, save_initialization, use_parallel, false); 
if show_init
    figure;
    imagesc(Cn, [0, 1]); colormap gray;
    hold on;
    plot(center(:, 2), center(:, 1), '.r', 'markersize', 10);
end
neuron_init = neuron.copy(); 
fprintf('initialize neuron ... is done\n')

%% estimation of background components 
neuron.update_background_parallel(use_parallel);
% neuron.update_background_batch(use_parallel); 
fprintf('estimation of background components is done\n')

%% pick neurons from the residual 
[center_res, Cn_res, PNR_res] =neuron.initComponents_residual_parallel(K, save_initialization, use_parallel);
neuron_init_res = neuron.copy(); 
fprintf('pick neurons from the residual is done\n')

%% estimate the temporal components 
neuron.update_temporal_parallel(use_parallel);
% neuron.update_temporal_batch(use_parallel);  
fprintf('estimate the temporal components is done\n')

%% estimate the spatial components 
neuron.update_spatial_parallel(use_parallel); 
% neuron.update_spatial_batch(use_parallel); 
fprintf('estimate the spatial components is done\n')

%% deal with bad neurons 
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.delete(tags>0); 
fprintf('deal with bad neurons is done\n')

%% merge neurons
neuron.merge_neurons_dist_corr(show_merge); 
neuron.merge_close_neighbors(show_merge, dmin_only); 
fprintf('merge neurons is done\n')

%% add a manual intervention and run the whole procedure for a second time 
neuron.orderROIs('snr');   % order neurons in different ways {'snr', 'decay_time', 'mean', 'circularity'}
if with_manual_intervention
    neuron.viewNeurons([], neuron.C_raw);
end
neuron.update_background_parallel(use_parallel); 
neuron.update_temporal_parallel(use_parallel); 
neuron.update_spatial_parallel(use_parallel); 

K = size(neuron.A,2); 
tags = neuron.tag_neurons_parallel();  % find neurons with fewer nonzero pixels than min_pixel and silent calcium transients
neuron.delete(tags>0);
neuron.merge_neurons_dist_corr(show_merge);
neuron.merge_close_neighbors(show_merge, dmin_only);
if K~=size(neuron.A,2)
    neuron.update_temporal_parallel(use_parallel);
    neuron.update_spatial_parallel(use_parallel);
    neuron.update_background_parallel(use_parallel);
end
fprintf('add a manual intervention and run the whole procedure for a second time is done\n')

%% save the workspace for future analysis 
neuron.save_workspace(); 
% neuron.save_workspace_batch(); 
fprintf('save the workspace for future analysis is done\n')

%% get the correlation image and PNR image for all neurons 
% neuron.correlation_pnr_parallel();
% neuron.correlation_pnr_batch(); 

% %% concatenate temporal components 
% neuron.concatenate_temporal_batch(); 
% neuron.viewNeurons([], neuron.C_raw);

%% show neuron contours  
Coor = neuron.show_contours(); 
fprintf('show neuron contours is done\n')

%% create a video for displaying 
% amp_ac = 5000; 
% range_ac = [0, 5000]; 
% range_Y = [0, 10000]; 
% neuron.show_demixed_video(save_demixed, kt, [], amp_ac, range_ac, range_Y); 
% fprintf('create a video for displaying is done\n')

%% save neurons shapes 
neuron.save_neurons(); 
fprintf('save neurons shapes is done\n')


