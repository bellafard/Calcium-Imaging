clear all; clc; close all;

addpath(genpath('Sources/NoRMCorre'));
addpath(genpath('utilities'));
addpath(genpath('deconvolution'));

gcp;

foldername = '/home/arash/s3/C143/2p/'; % folder where all the files are located.

filetype = 'raw';
files = subdir(fullfile(foldername,['*.',filetype]));
% FOV = size(read_file(files(1).name,1,1));
FOV = [512,512];
numFiles = length(files);

non_rigid = false;      % flag for non-rigid motion correction
output_type = 'tif';    % format to save registered files

if non_rigid; append = '_nr'; else; append = '_rig'; end    % use this to save motion corrected files

options_mc = NoRMCorreSetParms('d1',FOV(1),'d2',FOV(2),'grid_size',[128,128],'init_batch',200,...
                'overlap_pre',32,'mot_uf',4,'bin_width',200,'max_shift',24,'max_dev',8,'us_fac',50,...
                'output_type',output_type);

template = [];
col_shift = [];
for i = 1:numFiles
    fullname = files(i).name;
    [folder_name,file_name,ext] = fileparts(fullname);
    output_filename = fullfile(folder_name,[file_name,append,'.',output_type]);
    options_mc = NoRMCorreSetParms(options_mc,'output_filename',output_filename,'h5_filename','','tiff_filename',''); % update output file name
    
    [M,shifts,template,options_mc,col_shift] = normcorre_batch_even(fullname,options_mc,template);
    save(fullfile(folder_name,[file_name,'_shifts',append,'.mat']),'shifts','-v7.3'); % save shifts of each file at the respective folder        

    fprintf(['\n\n', file_name,' motion correction is done.\n\n'])
end

