function tomography_parameters_02;
% This function creates a structure containing relevant parameter
% information for processing tomographic piv data.  The data created by
% this function is inputted into 'reprojection_tomography_06' to perform
% reconstructions.  This version is similar to version 01 except that this
% version creates a structure that contains all information for running
% reprojection_tomography_06 including the calibration method and
% reconstruction method.

% This is the directory where the parameter structure is written
structure_pathname='/mnt/current_storage/Projects2/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Test_03/MatrixProduct_MLOS_Recon_01/';
% This is the parameter structure filename
structure_filename='tomography_parameters.mat';
% This is the domain of the volume to reconstruct in world coordinates

% xmin=-83.1966;
% xmax=82.6777;
% ymin=-89.9141;
% ymax=59.4309;
% zmin=-24.5;
% zmax=1.5;

xmin=-47.868291802704334;
xmax=47.116555027663708;
ymin=-74.607152275741100;
ymax=16.652798600494862;
zmin=-24.444600000000001;
zmax=1.513240000000000;

% This is the resolution of the reconstruction volume

% xres_volume=1426;
% yres_volume=1284;
% zres_volume=224;

xres_volume=817;
yres_volume=785;
zres_volume=224;

% This is the number of voxels to reconstruct at one time
voxel_recon_num=2e5;
% This is the range of frames to process
frame_min=1;
frame_max=10;
% This is the calibration method
calibration_method='cubic';
% This is the reconstruction method
reconstruction_method='multiplicative';
% These are the reconstruction parameters
reconstruction_parameters=1;
% This is a boolean value that tells MATLAB whether or not to save a *.mat
% version of the reconstruction (which may have issues on computers with
% less RAM)
mat_file_save=true;
% This is the directory which contains the world to image coordinate
% correspondences from the calibration data for each camera
calibration_directory='/mnt/current_storage/Projects2/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Test_03/MLOS_Reconstruction/Calibration/';
% This is the calibration directory prefix
calibration_dir_prefix='cam_';
% This is the directory containing the images to process for each camera
image_directory='/mnt/current_storage/Projects2/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Test_03/Preprocessed_Images/';
% This is the image directory prefix
image_dir_prefix='cam_';
% This is the directory to which the reconstructed volumes will be written
reconstruction_directory='/mnt/current_storage/Projects2/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Test_03/MatrixProduct_MLOS_Recon_01/Reconstruction/';
% This is the list of calibration directories
calibration_dir_list=dir([calibration_directory,calibration_dir_prefix,'*']);
% This is the list of image directories
image_dir_list=dir([calibration_directory,image_dir_prefix,'*']);
% This checks that there are the same number of calibration directories as
% there are image directories
if length(calibration_dir_list)==length(image_dir_list);
    camera_num=length(calibration_dir_list);
else;
    error('The number of calibration directories does not match the number of image directories . . . ');
end;
% This initializes the structure
parameter_data=struct;
% This iterates through the cameras filling in the structure
for ii=1:camera_num;
    % This is the directory path for the iith calibration directory
    parameter_data(ii).calibration_dir=[calibration_directory,calibration_dir_list(ii).name,'/'];
    % This is the directory path for the iith image directory
    parameter_data(ii).image_dir=[image_directory,image_dir_list(ii).name,'/'];
    % This creates the reconstruction directories
    if ii==1;
        mkdir(reconstruction_directory);
        mkdir([reconstruction_directory,image_dir_list(ii).name]);
        parameter_data(ii).reconstruction_dir=[reconstruction_directory,image_dir_list(ii).name,'/'];
    else;
        mkdir([reconstruction_directory,image_dir_list(ii).name]);
        parameter_data(ii).reconstruction_dir=[reconstruction_directory,image_dir_list(ii).name,'/'];
    end;
    % This is the full reconstruction directory
    mkdir([reconstruction_directory,'Full_Reconstruction/']);
    parameter_data(ii).full_reconstruction_dir=[reconstruction_directory,'Full_Reconstruction/'];
    % This is the list of images in the current image directory
    image_list=dir([parameter_data(ii).image_dir,'*.tif']);
    % The saves the number of images to the structure
    parameter_data(ii).image_num=length(image_list);
    % This loads an image to extract the image resolution
    I=imread([parameter_data(ii).image_dir,image_list(1).name]);
    % This extracts the resolution information
    [yres_cam,xres_cam]=size(I);
    % This saves the resolution information to the structure
    parameter_data(ii).image_res=[xres_cam,yres_cam];
    % This saves the volume domain to the structure
    parameter_data(ii).vol_domain=[xmin,xmax,ymin,ymax,zmin,zmax];
    % This saves the volume resoltuion to the structure
    parameter_data(ii).vol_res=[xres_volume,yres_volume,zres_volume];
    % This save the number of voxels to reconstruct in each slice of the
    % volume to the structure
    parameter_data(ii).voxel_recon_num=voxel_recon_num;
    % This saves the total number of cameras to the structure
    parameter_data(ii).camera_num=camera_num;
    % This saves the frame range to be processed
    parameter_data(ii).frame_domain=[frame_min,frame_max];
    % This saves the calibration method
    parameter_data(ii).cal_method=calibration_method;
    % This saves the reconstruction method
    parameter_data(ii).recon_method=reconstruction_method;
    % This saves the reconstruction method parameters vector (this is
    % currently only used for the exponent of the multiplicative method)
    parameter_data(ii).recon_param=reconstruction_parameters;
    % This saves the boolean value of whether or not to save a mat file
    % (the option of 'true' may not work well on computers with little RAM)
    parameter_data(ii).mat_file_save=mat_file_save;
end;
% This saves the tomography parameters structure to the drive
save([structure_pathname,structure_filename],'parameter_data');