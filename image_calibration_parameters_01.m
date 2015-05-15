function image_calibration_parameters_01;
% This function creates a structure containing all revelant data for
% simulating PIV images from multiple cameras including the camera
% calibration, the camera resolution, the volume to image, the resolution
% of the volume reconstruction, and the directories from which to read and
% to which to write all data.

% This is the directory to save the structures to
structure_pathname='/home/rlafoy/Documents/AEThER/Tomo_PIV/Function_Testing/';
% This is the filename of the structure to save
structure_filename='image_parameters.mat';
% This creates the image_parameters structure
image_parameters=create_image_structure;
% This saves the parameters structure
save([structure_pathname,structure_filename],'image_parameters');



function image_parameters=create_image_structure;
% This function creates the image structure used by the rest of this
% program to create simulated particle fields.

% This is the volume to be imaged
wrld_xmin=-1;
wrld_xmax=1;
wrld_ymin=-1;
wrld_ymax=1;
wrld_zmin=11.75;
wrld_zmax=12.25;
% This is a vector containing the volume coodinate range
v_domain=[wrld_xmin,wrld_xmax,wrld_ymin,wrld_ymax,wrld_zmin,wrld_zmax];
% This is the resolution of the reconstruction volume
xres_volume=512;
yres_volume=512;
zres_volume=128;
% This is a boolean value determining whether or not to save the 3D
% intensity field
intensity_field_save=true;
% This is a boolean value determining whether or not to save the 2D
% simulated camera images
camera_images_save=true;
% This is the simulated image resolution
xres_cam=512;
yres_cam=512;
% These are the numbers of cameras in the x and y directions
x_camera_number=2;
y_camera_number=2;
% This is the angle of the outer-most camera from the vertical axis (ie in
% a 4 x 4 grid this is the angle of cameras (1,1), (1,4), (4,1), and (4,4)
% from the center of the camera grid)
x_camera_theta=20;
y_camera_theta=20;
% This is a boolean value determining whether or not to display the camera
% calibration gui
display_calibration_gui=true;
% This is the standard deviation of the gaussian distribution at the
% reference plane
sigma_mean=1;
% This is the number of particles to simulate in the images (this number of
% particles will be used regardless of the size of the simulation domain)
prt_num=5e3;
% This is the number to seed the particle generator with
seed=0;
% This is the directory to write the image data to
image_directory='/home/rlafoy/Documents/AEThER/Tomo_PIV/Function_Testing/Images/';
% This is the directory to write the 3D intensity field to
intensity_directory='/home/rlafoy/Documents/AEThER/Tomo_PIV/Function_Testing/Intensity_Field/';
% This is the data directory to read the particle positions from
data_directory='/home/rlafoy/Documents/AEThER/Tomo_PIV/Field_Simulation/Data/Stokes_Vortex_Ring/Particle_Positions/';
% This is the calibration data directory
calibration_directory='/home/rlafoy/Documents/AEThER/Tomo_PIV/Function_Testing/Calibration/';
% This is the minimum frame to load
frame_min=1;
% This the maximum frame to load
frame_max=10;
% These are the frame step size
frame_step=1;

% This initializes the structure to save the data in
image_parameters=struct;
% This saves the volumne domain to the parameter structure
image_parameters.vol_domain=v_domain;
% This saves the reconstruction volume resolution to the structure
image_parameters.vol_res=[xres_volume,yres_volume,zres_volume];
% This save the boolean value determinging whether or not to save the 3D
% intensity field
image_parameters.intensity_save=intensity_field_save;
% This saves the boolean value determining whether or not to save the 2D
% simulated camera images
image_parameters.images_save=camera_images_save;
% This saves the image resolution to the structure
image_parameters.image_res=[xres_cam,yres_cam];
% This saves the camera numbers
image_parameters.camera_numbers=[x_camera_number,y_camera_number];
% This saves the camera angles
image_parameters.camera_angles=[x_camera_theta,y_camera_theta];
% This saves the camera calibration gui display boolean value
image_parameters.display_gui=display_calibration_gui;
% This saves the particle size to the structure
image_parameters.particle_diameter=sigma_mean;
% This saves the particle number to the structure
image_parameters.particle_number=prt_num;
% This is the random number seed for extracting particles
image_parameters.particle_seed=seed;
% This saves the image writing directory path to the structure
image_parameters.image_directory=image_directory;
% This saves the 3D intensity field directory to the structure
image_parameters.intensity_directory=intensity_directory;
% This saves the directory containging the particle positions to the
% structure
image_parameters.particle_directory=data_directory;
% This saves the directory containging the calibration data to the
% structure
image_parameters.calibration_directory=calibration_directory;
% This saves the vector containing the frames to create to the structure
image_parameters.frame_domain=[frame_min,frame_step,frame_max];



