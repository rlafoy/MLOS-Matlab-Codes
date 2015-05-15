function reprojection_tomography_08(parameter_data);
% This function is for generating tomographic reconstruction of volumes
% from multiple calibrated cameras.  This algorithm works similarly to the
% homographic tomography code except that this code is generallized by
% allowing the mapping function (x,y)=f(X,Y,Z) relating the world to image 
% coordinates to be arbitrary instead of a linear homography.  This
% function works by calculating the location in the image to which each 
% point in the tomographic volume maps to and then interpolates that point
% back onto the volume.  Initially an interpolation matrix must be
% calculated for each camera, but once this is completed, the reprojection
% step becomes a simple matrix mulitplication (followed by a reshape
% command).  This version is the same as version 05 except that it is more
% functionalized, all data is passed in as an argument to the overall
% function for automation purposes, and the reconstruction process is
% divided into the sections as defined by the voxel_recon_number (for
% computers with small amounts of RAM).
%
% For linear calibrations, only the first point_max=1000 points are
% randomly selected from the data set to avoid ou-of-memory errors.
%
% Things to work on: change scaling in multiplicative reconstruction so
% that the scaling is not done on each voxel section, change frame indexing
% so that the step size can be arbitrary (ie
% frame_min:frame_step:frame_max), parallelize the code, . . .
%
% Version 08 includes the option to perform a cubic calibration of the
% camera data.

% This initializes a timer
tic;

% This displays that the camera calibration is being generated
disp('Generating camera calibration interpolation file . . . ');
% % This calculates the linear calibration from the data
cal_struct=camera_calibration(parameter_data);

% This displays that the cameras are being reprojected
disp('Reprojecting images from each camera . . . ');
% This reprojects each camera
camera_reprojection(parameter_data,cal_struct);

% This displays that the full reconstruction is being generated
disp('Reconstructing full intensity field . . . ');
% This combines the reconstructions from each camera
full_reconstruction(parameter_data);

% If specified by parameter_data.mat_file_save, this saves the
% reconstructions as '*.mat' files
if parameter_data(1).mat_file_save;
    % This displays that the binary files are being converted
    disp('Converting binary reconstruction files to ''mat'' files . . . ');
    % This reads the reconstruction binary files and saves them as mat
    % files
    save_mat_file(parameter_data);
end;

% This returns the time of the timer
t=toc;
disp(['The total calculation time is ',num2str(t),' seconds . . . ']);



function cal_struct=camera_calibration(parameter_data);

% This performs the camera calibration

% This is a general calibration function that returns the calibration
% structure for later evaluation.  The type may be 'linear' or 'kriging'
% or whatever the user wishes to add.

% This is the calibration method
type=parameter_data(1).cal_method;
% This checks the format to generate the calibration structure
if strcmp(type,'linear');
    % This creates a linear calibration function
    cal_struct=linear_calibration(parameter_data);
elseif strcmp(type,'cubic');
    % This creates the cubic calibration function
    cal_struct=cubic_calibration(parameter_data);
elseif strcmp(type,'kriging');
    % This creates the kriging calibration function
    cal_struct=kriging_calibration(parameter_data);
end;



function camera_reprojection(parameter_data,cal_struct);

% This maps the reconstructed voxel world coordinates to image coordinates

% This loads the volume to reconstruct from the parameter data
volume_domain=parameter_data(1).vol_domain;
% This saves the volume extents in world coordinates to the domain
% variables
xmin=volume_domain(1);
xmax=volume_domain(2);
ymin=volume_domain(3);
ymax=volume_domain(4);
zmin=volume_domain(5);
zmax=volume_domain(6);
% This loads the volume resolution from the paramter data
volume_res=parameter_data(1).vol_res;
% This saves the resolutions to the x,y,z resolution variables
xres_volume=volume_res(1);
yres_volume=volume_res(2);
zres_volume=volume_res(3);
% This is the number of voxels to reconstruct at one time
voxel_recon_num=parameter_data(1).voxel_recon_num;
% This is the range of frames to be processed
frame_domain=parameter_data(1).frame_domain;
frame_min=frame_domain(1);
frame_max=frame_domain(2);
% This is the number of cameras
camera_num=parameter_data(1).camera_num;
% These are vectors of the volume coordinates
x_world_vect=linspace(xmin,xmax,xres_volume)';
y_world_vect=linspace(ymin,ymax,yres_volume)';
z_world_vect=linspace(zmin,zmax,zres_volume)';
% This is the total number of voxels
voxel_num=xres_volume*yres_volume*zres_volume;
% The total volume will be divided into voxel_recon_num long sections of
% voxel (which may or may not correspond with the planes) this is a vector
% of voxel indicies of each of these sections
voxel_index_vector=1:voxel_recon_num:voxel_num;
if voxel_index_vector(end)~=voxel_num;
    voxel_index_vector=[voxel_index_vector,voxel_num+1];
end;
% This iterates through the cameras
for n=1:camera_num;
    
    % This displays the current camera
    disp(['Processing camera number ',num2str(n),' . . . ']);
    
    % This is the filename of the camera calibration
    filename_read=[parameter_data(n).calibration_dir,'calibration_struct.mat'];
    % This loads the calibration structures into MATLAB
    load(filename_read);
    % This extracts the current camera's resolution
    image_res=parameter_data(n).image_res;
    xres_cam=image_res(1);
    yres_cam=image_res(2);
    % This is the list of images for the current camera
    image_list=dir([parameter_data(n).image_dir,'*.tif']);
    % If there are no images in the list for some silly reason (like a
    % crappy harddrive died for example), then the current camera is
    % skipped and the loop moves to the next camera
    if isempty(image_list);
        continue;
    end;
    % This creates a string that is passed into sprintf as an argument to
    % append the correct number of zeros onto the filename string
    sprintf_string = sprintf('%%s%%s%%0%0.0fd%%s',floor(log10(parameter_data(n).image_num))+1);
    % This iterates through the sections of voxels that are voxel_recon_num long
    for ii=1:length(voxel_index_vector)-1;
        
        if mod(ii,10)==0;
            % This displays the current voxel volume being processed
            disp(['Processing voxel section ',num2str(ii),' of ',num2str(length(voxel_index_vector)-1),' sections . . . ']);
        end;
        
        % This is the minimum index
        index_min=voxel_index_vector(ii);
        % This is the maximum index
        index_max=voxel_index_vector(ii+1)-1;
        % This extracts the vectors of voxel world coordinates
        [x_world_slice,y_world_slice,z_world_slice]=meshgrid_indexing_01(x_world_vect,y_world_vect,z_world_vect,index_min:index_max);
        % This returns the image coordinates of the current world coordinates
        [x_image,y_image]=calibration_evaluate(cal_struct(n),x_world_slice,y_world_slice,z_world_slice);
        % This creates the interpolation matrix
        M=interp_matrix_01(x_image,y_image,xres_cam,yres_cam);
        % This iterates through the frames to be processed
        for jj=frame_min:frame_max;
            % This is the filename of the image to be loaded
            filename_read=[parameter_data(n).image_dir,image_list(jj).name];
            % This loads the current image
            I=imread(filename_read);
            % This performs the reprojection interpolation
            V_Reproj=M*double(I(:));
            % This is the filename for writing using the magical sprintf
            % function
            filename_write=sprintf(sprintf_string,parameter_data(n).reconstruction_dir,'frame_',jj,'.dat');
            % This checks whether a copy of the reconstruction data already
            % exists for the current frame and deletes it if so (this is so
            % that the data is not simply appended onto the end of an
            % already existing file)
            if (ii==1)&&(exist(filename_write,'file'));
                delete(filename_write);
            end;
            % This opens the file for writing
            fid=fopen(filename_write,'a');
            % This writes the data
            fwrite(fid,V_Reproj,'double');
            % This closes the file
            fclose(fid);
        end;
    end;
end;



function full_reconstruction(parameter_data);

% This section loads and combines the reprojections

% Once the volumes have been reprojected from each camera they are combined
% together in some magical way to make the full 3D intensity field; this
% iterates through the frames

% This loads the volume resolution from the paramter data
volume_res=parameter_data(1).vol_res;
% This saves the resolutions to the x,y,z resolution variables
xres_volume=volume_res(1);
yres_volume=volume_res(2);
zres_volume=volume_res(3);
% This is the reconstruction method
reconstruction_method=parameter_data(1).recon_method;
% These are the reconstruction parameters (which are currently only used
% for the exponent of the multiplicative method)
reconstruction_parameters=parameter_data(1).recon_param;
% This is the number of voxels to reconstruct at one time
voxel_recon_num=parameter_data(1).voxel_recon_num;
% This is the range of frames to be processed
frame_domain=parameter_data(1).frame_domain;
frame_min=frame_domain(1);
frame_max=frame_domain(2);
% This is the number of cameras
camera_num=parameter_data(1).camera_num;
% This is the total number of voxels
voxel_num=xres_volume*yres_volume*zres_volume;
% The total volume will be divided into voxel_recon_num long sections of
% voxels (which may or may not correspond with the planes) this is a vector
% of voxel indicies of each of these sections
voxel_index_vector=1:voxel_recon_num:voxel_num;
if voxel_index_vector(end)~=voxel_num;
    voxel_index_vector=[voxel_index_vector,voxel_num+1];
end;
% This creates a string that is passed into sprintf as an argument to
% append the correct number of zeros onto the filename string
sprintf_string = sprintf('%%s%%s%%0%0.0fd%%s',floor(log10(parameter_data(1).image_num))+1);
% This iterates through the frames reconstructing the volumes
for jj=frame_min:frame_max;
    % This iterates through the slices for each frame
    for ii=1:length(voxel_index_vector)-1;
        % This is the minimum index
        index_min=voxel_index_vector(ii);
        % This is the maximum index
        index_max=voxel_index_vector(ii+1)-1;
        % This calculates the offset to use for the fseek command (from the
        % beginning of the file in bytes; since double precision is being
        % used to store data, there are 8 bytes for every data vector
        % component)
        fileread_offset=8*(index_min-1);
        % This is the numer of elements to read in at one time
        fileread_size=index_max-index_min+1;
        % This iterates through each camera loading the current cameras
        % reprojection volume
        for n=1:camera_num;
            % This initializes the reconstruction vector
            if n==1;
                V_Initial='initialize';
            end;
            % This is the filename for writing using the magical sprintf
            % function
            filename_read=sprintf(sprintf_string,parameter_data(n).reconstruction_dir,'frame_',jj,'.dat');
            % If there are no reconstructions in the list for some silly reason (like a
            % crappy harddrive died for example), then the current camera is
            % skipped and the loop moves to the next camera
            if not(exist(filename_read,'file'));
                continue;
            end;
            % This opens the file for reading
            fid=fopen(filename_read);
            % This moves the file pointer to the start of the current voxel
            % section counting from the beginning of file
            fseek(fid,fileread_offset,'bof');
            % This loads the current slice of the reprojection volume
            V_Reproj=fread(fid,fileread_size,'double');            
            % This closes the file
            fclose(fid);
            % This recombines the reprojection with the total
            % reconstruction vector
            V_Initial=volume_reconstruct(V_Initial,V_Reproj,camera_num,reconstruction_method,reconstruction_parameters);
        end;
        % These are the leading zeros to append to the filename
        zero_string=repmat('0',1,floor(log10(parameter_data(1).image_num))-floor(log10(jj)));
        % This is the filename to write the reconstruction data to
        filename_write=[parameter_data(1).full_reconstruction_dir,'frame_',zero_string,num2str(jj),'.dat'];
        % This checks whether a copy of the reconstruction data already
        % exists for the current frame and deletes it if so (this is so
        % that the data is not simply appended onto the end of an
        % already existing file)
        if (ii==1)&&(exist(filename_write,'file'));
            delete(filename_write);
        end;
        % This opens the file for writing
        fid=fopen(filename_write,'a');
        % This writes the data
        fwrite(fid,V_Initial,'double');
        % This closes the file
        fclose(fid);
    end;
end;



function save_mat_file(parameter_data);
% This function loads the binary data file created by full_reconstruction, 
% converts the file to a properly shaped variable 'I', and then saves the
% data to a '*.mat' file in the same directory as 'Full_Reconstruction'.
% Since the whole file may be very large, this function may cause memory
% issues.

% This loads the volume resolution from the paramter data
volume_res=parameter_data(1).vol_res;
% This saves the resolutions to the x,y,z resolution variables
xres_volume=volume_res(1);
yres_volume=volume_res(2);
zres_volume=volume_res(3);
% This is the range of frames to be processed
frame_domain=parameter_data(1).frame_domain;
frame_min=frame_domain(1);
frame_max=frame_domain(2);
% This iterates through the frames saving the reconstructions as *.mat
% files
for ii=frame_min:frame_max;
    % These are the leading zeros to append to the filename
    zero_string=repmat('0',1,floor(log10(parameter_data(1).image_num))-floor(log10(ii)));
    % This is the filename to read the reconstruction data from
    filename_read=[parameter_data(1).full_reconstruction_dir,'frame_',zero_string,num2str(ii),'.dat'];
    % This is the filename to write the reconstruction data to
    filename_write=[parameter_data(1).full_reconstruction_dir,'frame_',zero_string,num2str(ii),'.mat'];
    % This opens the file for reading
    fid=fopen(filename_read,'r');
    % This writes the data
    I=fread(fid,'double');
    % This closes the file
    fclose(fid);
    % This reshapes the matrix to the appropriate size
    I=reshape(I,[yres_volume,xres_volume,zres_volume]);
    % This saves the data as a *.mat file format
    save(filename_write,'I','-v7.3');
end;



function V_Recon=volume_reconstruct(V_Initial,V_Reproj,camera_num,method,method_parameters);
% This function combines the reprojection volumes from each camera to
% produce the total 3D intensity field.  The variable 'V_Initial' is the 
% current reconstruction, the variable 'V_Reproj' is a vector containing 
% the current reprojection, 'camera_num' is the number of cameras, and 
% 'method' is a string giving the reconstruction method.  The variable 
% 'method' may be
%
%   'additive'          This method recombines the images by averaging the
%                       intensities together.
%   'multiplicative'    This method multiplies the images together raised
%                       to the 1/camera_num power.
%   'minimization'      This method takes the final intensity as the
%                       minimum intensity of all cameras.
%
% To iteratively perform the reconstruction, the output of the previous
% function call is input into the current function call.

% This implements the additive reconstruction method
if strcmp(method,'additive')==1;
    % This checks whether to initialize the reconstruction
    if strcmp(V_Initial,'initialize');
        V_Initial=zeros(size(V_Reproj));
    end;
    % This adds the current reprojection onto the total reprojection
    V_Recon=V_Initial+V_Reproj/camera_num;
% This implements the multiplicative reconstruction method    
elseif strcmp(method,'multiplicative')==1;
    % This is the exponent for the reconstruction
    gamma=method_parameters(1);
    % This checks whether to initialize the reconstruction
    if strcmp(V_Initial,'initialize');
        V_Initial=ones(size(V_Reproj));
    end;
    % This scales the current reprojection to lie within the range of [0,1]
    %V_Reproj=(V_Reproj-min(V_Reproj(:)))/(max(V_Reproj(:))-min(V_Reproj(:)));
    
    V_Reproj=(V_Reproj-0)/(2^16-1-0);
    
    % This multiplies the current reconstruction to the total
    % reconstruction
    V_Recon=V_Initial.*(V_Reproj.^(gamma));
% This implments the minimization reconstruction method
elseif strcmp(method,'minimization')==1;
    % This checks whether to initialize the reconstruction
    if strcmp(V_Initial,'initialize');
        V_Initial=inf(size(V_Reproj));
    end;
    % This takes the minimum value of the current reprojection and the
    % total reprojection
    V_Recon=min(V_Initial,V_Reproj);
end;
    


function x=kriging_calibration_eval(C,X,Y,Z);
% This function calculates the coordinate transformation x_i=f(X,Y,Z) for a
% set of world coordinates [X,Y,Z] and a the kriging camera calibration
% structure C using the kriging camera calibration interpolation.  Only the
% variable 'x' is used for the camera coordinates in this function, but it
% is called to generate both the x and y camera coordiantes by the function
% 'calibration_evaluate'.

% This first tries the function to see if there is sufficient memory
try;
    % This evaluates the coordinate transformation
    x=predictor([X,Y,Z],C);
catch ErrorID;
    % This checks if the error was 'Out of memory.'
    if strcmp(ErrorID.identifier,'MATLAB:nomem')||strcmp(ErrorID.identifier,'MATLAB:pmaxsize');
        % This is the index at roughly the half the length of the vectors
        half_index=round(length(X)/2);
        % This breaks the vectors in half to reduce the memory requirements
        X1=X(1:half_index);
        X2=X(half_index+1:end);
        Y1=Y(1:half_index);
        Y2=Y(half_index+1:end);
        Z1=Z(1:half_index);
        Z2=Z(half_index+1:end);
        % This attempts to calculate the interpolations with vectors half of
        % the current length
        x1=kriging_calibration_eval(C,X1,Y1,Z1);
        x2=kriging_calibration_eval(C,X2,Y2,Z2);
        % This recombines the vectors for output
        x=[x1;x2];
    else;
        % If the error was something other than running out of memory, the
        % function is halted and the error is printed to the terminal
        error(ErrorID.message);
    end;
end;



function [x,y]=calibration_evaluate(cal_struct,X,Y,Z);
% This function evaluates the mapping of the world coordinates (X,Y,Z) to 
% the camera coordinates (x,y) using the calibration function described in
% cal_struct.

% This checks the format of the cabration structure
if strcmp(cal_struct(1).type,'linear');
    % This evaluates the world coordinate mapping
    [x,y]=linear_calibration_eval(cal_struct.x_cal,X,Y,Z);
elseif strcmp(cal_struct(1).type,'cubic');
    % This extracts the cubic calibration function structure
    calibration_function=cal_struct.calibration_function;
    % This evaluates the world coordinate mapping
    [x,y]=cubic_calibration_eval_01(calibration_function,X,Y,Z);
elseif strcmp(cal_struct(1).type,'kriging');
    % This evaluates the world coordinate mapping
    x=kriging_calibration_eval(cal_struct.x_cal,X,Y,Z);
    y=kriging_calibration_eval(cal_struct.y_cal,X,Y,Z);
end;



function cal_struct=kriging_calibration(parameter_data);
% This function returns a calibration structure using the kriging camera
% calibration model.  The form of the structure is
%
%  cal_struct.type='kriging'
%  cal_struct.x_cal=C_x;
%  cal_struct.y_cal=C_y;
%
% The index of the structure corresponds to the camera index.

% This initializes the calibration structure
cal_struct=struct([]);
% This iterates through the cameras calculating the calibrations
for n=1:parameter_data(1).camera_num;
    % This is the calibration pathname to read into MATLAB
    pathname_read=parameter_data(n).calibration_dir;
    % This is the calibration filename to read into MATLAB
    filename_read=dir([pathname_read,'*_calibration.mat']);
    % This loads the nth camera's point correspondence vectors which are of
    % the form [x,y]=F([X,Y,Z])
    load([pathname_read,filename_read.name]);
    % These are the initial values of theta
    theta=[10,10,10];
    % This is the lower bound on theta
    lower_bound=[1e-1,1e-1,1e-1];
    % This is the upper bound on theta
    upper_bound=[30,30,30];
    % This calculates the Kriggin model for the x-coordinate transformation
    [C_x,perf_x]=dacefit([X,Y,Z],x,@regpoly2,@corrgauss,theta,lower_bound,upper_bound);
    % This calculates the Kriggin model for the x-coordinate transformation
    [C_y,perf_y]=dacefit([X,Y,Z],y,@regpoly2,@corrgauss,theta,lower_bound,upper_bound);
    % This adds the data into the calibration structure
    cal_struct(n).x_cal=C_x;
    cal_struct(n).y_cal=C_y;
    % This saves the calibration type into the structure
    cal_struct(n).type='kriging';
    % This is the filename to save the calibration to
    filename_save=[pathname_read,'calibration_struct.mat'];
    % This save the calibration data to the drive
    save(filename_save,'cal_struct');
end;



function cal_struct=cubic_calibration(parameter_data);
% This function returns a calibration structure using a cubic calibration
% of the data in (X,Y) with a linear interpolation between the Z planes.

% This initializes the calibration structure
cal_struct=struct([]);
% This iterates through the cameras calculating the calibrations
for n=1:parameter_data(1).camera_num;
    % This is the calibration pathname to read into MATLAB
    pathname_read=parameter_data(n).calibration_dir;
    % This is the calibration filename to read into MATLAB
    filename_read=dir([pathname_read,'*_calibration.mat']);
    % This loads the nth camera's point correspondence vectors which are of
    % the form [x,y]=F([X,Y,Z])
    load([pathname_read,filename_read.name]);
    % This eliminates Z planes that are too close together
    [x,y,X,Y,Z]=eliminate_close_z_planes(x,y,X,Y,Z);
    % This groups the data into cell array based upon the Z plane, ie all cell
    % arrays of the same index are on the same Z plane
    [x_plane,y_plane,X_plane,Y_plane,Z_plane]=group_z_planes(x,y,X,Y,Z);
    % This calculates the calibration function and saves the function to the
    % calibration structure
    calibration_function=calculate_calibration_function(x_plane,y_plane,X_plane,Y_plane,Z_plane);
    % This saves the calibration function into the external calibration
    % structure
    cal_struct(n).calibration_function=calibration_function;
    % This saves the calibration type into the structure
    cal_struct(n).type='cubic';
    % This is the filename to save the calibration to
    filename_save=[pathname_read,'calibration_struct.mat'];
    % This save the calibration data to the drive
    save(filename_save,'cal_struct');
end;



function [x,y]=cubic_calibration_eval_01(calibration_function,X,Y,Z);
% This function evaluates the calibration function given by the structure
% using a cubic polynomial (in X and Y) and a linear interpolation (in Z).
% The function assigns each world coodinate (X,Y,Z) to one of five
% conditions given by the following:
%
%  1: D_Sum = 0             & size(D,2) = 1 (Extrapolation)
%  2: D_Sum = 0             & size(D,2) > 1 (Extrapolation)
%  3: 0 < D_Sum < size(D,2) & size(D,2) > 1 (Interpolation)
%  4: D_Sum = size(D,2)     & size(D,2) = 1 (Extrapolation)
%  5: D_Sum = size(D,2)     & size(D,2) > 1 (Extrapolation)
%
% For each of these five conditions the nearest calibration planes are
% found and linear weighting coefficients are calculated by the following:
%
%  1: w1 = 1 & w2 = 0       & plane_1 = 1           & plane_2 = 1
%  2: w1 = 1 & w2 = 0       & plane_1 = 1           & plane_2 = 1
%  3: w1 = L & w2 = 1 - L   & plane_1 = D_Sum       & plane_2 = D_Sum + 1
%  4: w1 = 0 & w2 = 1       & plane_1 = size(D,1)   & plane_2 = size(D,1))
%  5: w1 = 0 & w2 = 1       & plane_1 = size(D,1)   & plane_2 = size(D,1))
%
% The world coordinates (X,Y,Z) can have any dimensions.  The output image
% coordinates (x,y) have the same dimensions as the input world
% coordinates.

% These are the dimensions of the (X,Y,Z) world coordinates
[ii_dim,jj_dim,kk_dim]=size(X);

% This converts the world coordinates (X,Y,Z) to row vectors for processing
X=X(:)';
Y=Y(:)';
Z=Z(:)';

% This is a vector of the Z planes to determine which planes' calibration
% functions to use
Z_plane=zeros(size(calibration_function.z_world_value));
for ii=1:length(calibration_function.z_world_value);
    Z_plane_temp=calibration_function.z_world_value{ii};
    Z_plane(ii)=Z_plane_temp(1);
end;

% This creates matrices to calculate the distances between the evaluation
% points and each Z plane
[Z_Eval,Z_Planes]=meshgrid(Z,Z_plane);
% This is the distance matrix
D=Z_Eval-Z_Planes;
% These are the positive distances
D_Sum=sum(D>0,1);

% These are the subscripted indices (in the dimension 1 direction) of the
% first and second planes to interpolate between
plane_1_sub_index=D_Sum;
plane_2_sub_index=D_Sum+1;

% These are the interpolation points that are before the first calibration
% plane (and thus need to be extrapolated to)
plane_index_zero_length=(D_Sum==0);
% These are the interpolation points that are after the last calibration
% plane (and thus need to be extrapolated to)
plane_index_full_length=(D_Sum==size(D,1));

% Conditions 1 & 2: If the world coodinates are before the first plane,
% then the calibration planes are both set to the first plane
plane_1_sub_index(plane_index_zero_length)=1;
plane_2_sub_index(plane_index_zero_length)=1;

% Conditions 4 & 5: If the world coordinates are after the last plane, then
% the calibration planes are both set to the last plane
plane_1_sub_index(plane_index_full_length)=size(D,1);
plane_2_sub_index(plane_index_full_length)=size(D,1);

% Condition 3: If the world coordinates are between the first and last
% planes (inclusive) then the linear indices of the distances to these
% planes is calculated
plane_1_lin_index=sub2ind(size(D),plane_1_sub_index,1:size(D,2));
plane_2_lin_index=sub2ind(size(D),plane_2_sub_index,1:size(D,2));

% Condition 3: These are the weighting coefficients (some of which may be equal to
% infinity for extrapolated points)
w1=abs(D(plane_2_lin_index))./(D(plane_1_lin_index)-D(plane_2_lin_index));
w2=abs(D(plane_1_lin_index))./(D(plane_1_lin_index)-D(plane_2_lin_index));

% Condition 1 & 2: For points before the first calibration plane, this sets
% the weighting coefficients to use the calibration function of the first
% plane. This will remove some of the infinite weight coefficients assigned
% for condition 3.
w1(plane_index_zero_length)=1;
w2(plane_index_zero_length)=0;

% Condition 4 & 5: For points after the last calibration plane, this sets
% the weighting coefficients to use the calibration function of the last
% plane.  This will remove the rest of the infinite weight coefficients 
% assigned for condition 3.
w1(plane_index_full_length)=0;
w2(plane_index_full_length)=1;

% These are the linearly weighted x polynomial coefficients
Cx=bsxfun(@times,w1',calibration_function.x_coefficients(plane_1_sub_index,:))+bsxfun(@times,w2',calibration_function.x_coefficients(plane_2_sub_index,:));
% These are the linearly weighted y polynomial coefficients
Cy=bsxfun(@times,w1',calibration_function.y_coefficients(plane_1_sub_index,:))+bsxfun(@times,w2',calibration_function.y_coefficients(plane_2_sub_index,:));

% These are the current interpolated mean world coordinates
X_Mean=w1.*calibration_function.x_world_mean(plane_1_sub_index)+w2.*calibration_function.x_world_mean(plane_2_sub_index);
Y_Mean=w1.*calibration_function.y_world_mean(plane_1_sub_index)+w2.*calibration_function.y_world_mean(plane_2_sub_index);
% These are the current interpolated standard deviation of the world coordinates
X_Std=w1.*calibration_function.x_world_std(plane_1_sub_index)+w2.*calibration_function.x_world_std(plane_2_sub_index);
Y_Std=w1.*calibration_function.y_world_std(plane_1_sub_index)+w2.*calibration_function.y_world_std(plane_2_sub_index);
										
% This calculates the scaled X value
X_Scale=((X-X_Mean)./X_Std)';
% This calculates the scaled Y value
Y_Scale=((Y-Y_Mean)./Y_Std)';

% This is the polynomial independent vector
XY_Vect=[ones(size(X_Scale)),X_Scale,X_Scale.^2,X_Scale.^3,Y_Scale,Y_Scale.^2,Y_Scale.^3,X_Scale.*Y_Scale,X_Scale.^2.*Y_Scale,X_Scale.*Y_Scale.^2];																					

% This is the x image value
x=sum(Cx.*XY_Vect,2)';
% This is the y image value
y=sum(Cy.*XY_Vect,2)';

% This converts the dimensions of the image coordinates (x,y) to the same
% dimensions as the input world coordinates (X,Y,Z)
x=reshape(x,[ii_dim,jj_dim,kk_dim]);
y=reshape(y,[ii_dim,jj_dim,kk_dim]);



function calibration_function=calculate_calibration_function(x_plane,y_plane,X_plane,Y_plane,Z_plane);

% This initializes a data structure for saving the cubic calibration
% function into
calibration_function=struct;

% These are the minimization option parameters
min_options=optimset('MaxIter',1e4,'MaxFunEvals',1e4);

% This is the number of Z values
Z_Number=length(x_plane);
% This iterates through the Z values creating the calibration functions
for ii=1:Z_Number;
    
    % This saves the current Z value to the calibration structure
    calibration_function.z_world_value{ii}=Z_plane{ii};
    
    % This is the mean of the X and Y coordinates
    X_Mean=mean(X_plane{ii});
    Y_Mean=mean(Y_plane{ii});
    % These are the zero-mean X and Y coordinates
    X_Zero_Mean=X_plane{ii}-X_Mean;
    Y_Zero_Mean=Y_plane{ii}-Y_Mean;
    % This is the average X and Y distance from the origin
    X_Std=mean(abs(X_Zero_Mean));
    Y_Std=mean(abs(Y_Zero_Mean));
    % This scales the data to have an average distance from the origin of 1
    X_Scaled=X_Zero_Mean/X_Std;
    Y_Scaled=Y_Zero_Mean/Y_Std;

    % This save the mean values and standard deviation values of the world
    % data to the data structure
    calibration_function.x_world_mean(ii)=X_Mean;
    calibration_function.y_world_mean(ii)=Y_Mean;
    calibration_function.x_world_std(ii)=X_Std;
    calibration_function.y_world_std(ii)=Y_Std;  
    
    % These are the initial x coefficient guesses for the linear fit.  The
    % costant coefficient term is set to the mean data value and the
    % remaining terms are set to zero.
    Cx0=[mean(x_plane{ii});zeros(2,1)];
    % These are the initial y coefficient guesses for the linear fit.  The
    % costant coefficient term is set to the mean data value and the
    % remaining terms are set to zero.
    Cy0=[mean(y_plane{ii});zeros(2,1)];
    
    % This fits the linear polynomial to the x data
    [x_lin_coeff_vector,~]=fminsearch(@(C)linear_fit_residual(C,x_plane{ii},X_Scaled,Y_Scaled),Cx0,min_options);
   
    % This fits the linear polynomial to the y data
    [y_lin_coeff_vector,~]=fminsearch(@(C)linear_fit_residual(C,y_plane{ii},X_Scaled,Y_Scaled),Cy0,min_options);

    % These are the initial x coefficient guesses for the quadratic fit.  The
    % linear terms in the quadratic polynomial are taken from the
    % coefficients found in the linear fit.  The remaining coefficients are
    % set to zero.
    Cx0(1)=x_lin_coeff_vector(1);
    Cx0(2)=x_lin_coeff_vector(2);
    Cx0(3)=0;
    Cx0(4)=x_lin_coeff_vector(3);
    Cx0(5)=0;
    Cx0(6)=0;
    % These are the initial y coefficient guesses for the quadratic fit.  The
    % linear terms in the quadratic polynomial are taken from the
    % coefficients found in the linear fit.  The remaining coefficients are
    % set to zero.
    Cy0(1)=y_lin_coeff_vector(1);
    Cy0(2)=y_lin_coeff_vector(2);
    Cy0(3)=0;
    Cy0(4)=y_lin_coeff_vector(3);
    Cy0(5)=0;
    Cy0(6)=0;
    
    % This fits the quadratic polynomial to the x data
    [x_quad_coeff_vector,~]=fminsearch(@(C)quadratic_fit_residual(C,x_plane{ii},X_Scaled,Y_Scaled),Cx0,min_options);
 
    % This fits the quadratic polynomial to the y data
    [y_quad_coeff_vector,~]=fminsearch(@(C)quadratic_fit_residual(C,y_plane{ii},X_Scaled,Y_Scaled),Cy0,min_options);
    
    % These are the initial x coefficient guesses for the cubic fit.  The
    % quadratic terms in the cubic polynomial are taken from the
    % coefficients found in the quadratic fit.  The remaining coefficients are
    % set to zero.
    Cx0(1)=x_quad_coeff_vector(1);
    Cx0(2)=x_quad_coeff_vector(2);
    Cx0(3)=x_quad_coeff_vector(3);
    Cx0(4)=0;
    Cx0(5)=x_quad_coeff_vector(4);
    Cx0(6)=x_quad_coeff_vector(5);
    Cx0(7)=0;
    Cx0(8)=x_quad_coeff_vector(6);
    Cx0(9)=0;
    Cx0(10)=0;
    % These are the initial y coefficient guesses for the cubic fit.  The
    % quadratic terms in the cubic polynomial are taken from the
    % coefficients found in the quadratic fit.  The remaining coefficients are
    % set to zero.
    Cy0(1)=y_quad_coeff_vector(1);
    Cy0(2)=y_quad_coeff_vector(2);
    Cy0(3)=y_quad_coeff_vector(3);
    Cy0(4)=0;
    Cy0(5)=y_quad_coeff_vector(4);
    Cy0(6)=y_quad_coeff_vector(5);
    Cy0(7)=0;
    Cy0(8)=y_quad_coeff_vector(6);
    Cy0(9)=0;
    Cy0(10)=0;
  
    % This fits the cubic polynomial to the x data
    [x_coeff_vector,~]=fminsearch(@(C)cubic_fit_residual(C,x_plane{ii},X_Scaled,Y_Scaled),Cx0,min_options);
    % This saves the x coefficient vector to the data structure
    calibration_function.x_coefficients(ii,:)=x_coeff_vector';

    % This fits the cubic polynomial to the y data
    [y_coeff_vector,~]=fminsearch(@(C)cubic_fit_residual(C,y_plane{ii},X_Scaled,Y_Scaled),Cy0,min_options);
    % This saves the y coefficient vector to the data structure
    calibration_function.y_coefficients(ii,:)=y_coeff_vector';

end;



function r=cubic_fit_residual(C,x,X,Y);
% This function calculates the sum of the squares of the error of the cubic
% polynomial fit.

% This extracts the coefficents from the coefficient vector
CX0Y0=C(1);
CX1Y0=C(2);
CX2Y0=C(3);
CX3Y0=C(4);
CX0Y1=C(5);
CX0Y2=C(6);
CX0Y3=C(7);
CX1Y1=C(8);
CX2Y1=C(9);
CX1Y2=C(10);
% This is the polynomial evaluated at (X,Y)
xP = CX0Y0 + ...
    CX1Y0 * X + CX2Y0 * X.^2 + CX3Y0 * X.^3 + ...
    CX0Y1 * Y + CX0Y2 * Y.^2 + CX0Y3 * Y.^3 + ...
    CX1Y1 * X .* Y + CX2Y1 * X.^2 .* Y + CX1Y2 * X .* Y.^2;
% This is the sum of the squares of the error
r=sum((x-xP).^2);



function r=quadratic_fit_residual(C,x,X,Y);
% This function calculates the sum of the squares of the error of the
% quadratic polynomial fit.

% This extracts the coefficents from the coefficient vector
CX0Y0=C(1);
CX1Y0=C(2);
CX2Y0=C(3);
CX3Y0=0;
CX0Y1=C(4);
CX0Y2=C(5);
CX0Y3=0;
CX1Y1=C(6);
CX2Y1=0;
CX1Y2=0;
% This is the polynomial evaluated at (X,Y)
xP = CX0Y0 + ...
    CX1Y0 * X + CX2Y0 * X.^2 + CX3Y0 * X.^3 + ...
    CX0Y1 * Y + CX0Y2 * Y.^2 + CX0Y3 * Y.^3 + ...
    CX1Y1 * X .* Y + CX2Y1 * X.^2 .* Y + CX1Y2 * X .* Y.^2;
% This is the sum of the squares of the error
r=sum((x-xP).^2);



function r=linear_fit_residual(C,x,X,Y);
% This function calculates the sum of the squares of the error of the
% linear polynomial fit.

% This extracts the coefficents from the coefficient vector
CX0Y0=C(1);
CX1Y0=C(2);
CX2Y0=0;
CX3Y0=0;
CX0Y1=C(3);
CX0Y2=0;
CX0Y3=0;
CX1Y1=0;
CX2Y1=0;
CX1Y2=0;
% This is the polynomial evaluated at (X,Y)
xP = CX0Y0 + ...
    CX1Y0 * X + CX2Y0 * X.^2 + CX3Y0 * X.^3 + ...
    CX0Y1 * Y + CX0Y2 * Y.^2 + CX0Y3 * Y.^3 + ...
    CX1Y1 * X .* Y + CX2Y1 * X.^2 .* Y + CX1Y2 * X .* Y.^2;
% This is the sum of the squares of the error
r=sum((x-xP).^2);



function [x,y,X,Y,Z]=eliminate_close_z_planes(x,y,X,Y,Z);
% This function eliminates data of Z planes that are too close together;
% this could be due to round-off error in the computer when saving the
% calibration data.

% This is the tolerance for determining the unique Z planes
Z_Tol=1e-6;

% These are the unique values of Z to within machine tolerance
Z_Unique=unique(Z);
% This sorts the unique values of Z to determine the distance between
% adjacent values
Z_Unique_Sort=sort(Z_Unique,'ascend');
% This is the difference between adjacent Z values
Z_Diff=diff(Z_Unique_Sort);
% These are the indices of the Z value differences that are less then Z_Tol
Z_Nonunique_Index=find(Z_Diff<Z_Tol);

% This finds the points to remove if any Z values are too close together
if any(Z_Nonunique_Index);
    % This initializes the vector of coordinate point indices to remove
    remove_index=[];
    % These are the first non-unique value of Z
    Z_Values_One=Z_Unique_Sort(Z_Nonunique_Index);
    % These are the second non-unique value of Z
    Z_Values_Two=Z_Unique_Sort(Z_Nonunique_Index+1);
    % This iterates through the Z values that have differences less then Z_Tol
    for ii=1:length(Z_Values_One);
        % These are the indices of Z values equal to the first Z value in the
        % current difference
        Z_Value_One_Index=find(Z==Z_Values_One(ii));
        % These are the indices of Z values equal to the secod Z value in the
        % current difference
        Z_Value_Two_Index=find(Z==Z_Values_Two(ii));
        % This adds the indices of the smaller group of points to the
        % remove_index vector
        if numel(Z_Value_One_Index)>=numel(Z_Value_Two_Index);
            % This adds the indices of the second Z value to the vector of
            % points to remove
            remove_index=[remove_index;Z_Value_Two_Index];
        else;
            % This adds the indices of the first Z value to the vector of
            % points to remove
            remove_index=[remove_index;Z_Value_One_Index];
        end;
    end;
    % This removes the "non-unique" points from the calibration data
    x(remove_index)=[];
    y(remove_index)=[];
    X(remove_index)=[];
    Y(remove_index)=[];
    Z(remove_index)=[];
end;



function [x_plane,y_plane,X_plane,Y_plane,Z_plane]=group_z_planes(x,y,X,Y,Z);
% This function finds the unique Z values of the data and groups the data
% into cell arrays with all the same Z values.

% These are the unique Z values
Z_Unique=unique(Z);
% This is the number of unique Z values
Z_Number=length(Z_Unique);
% This initializes the coordinate cell arrays
x_plane=cell(Z_Number,1);
y_plane=cell(Z_Number,1);
X_plane=cell(Z_Number,1);
Y_plane=cell(Z_Number,1);
Z_plane=cell(Z_Number,1);
% This iterates through the different Z values saving the dat to the cell
% arrays
for ii=1:Z_Number;
    % These are the indices of the current Z values
    Current_Index=find(Z==Z_Unique(ii));
    % This saves these coordinates to the cell arrays
    x_plane{ii}=x(Current_Index);
    y_plane{ii}=y(Current_Index);
    X_plane{ii}=X(Current_Index);
    Y_plane{ii}=Y(Current_Index);
    Z_plane{ii}=Z(Current_Index);
end;



function cal_struct=linear_calibration(parameter_data);
% This function returns a calibration structure using a linear camera
% calibration model.  The form of the structure is
%
%  cal_struct.type='linear'
%  cal_struct.x_cal=P;
%  cal_struct.y_cal=P;
%
% The index of the structure corresponds the camera index.

% This is the maximum number of calibration points
point_max=1000;
% This initializes the calibration structure
cal_struct=struct([]);
% This iterates through the cameras calculating the calibrations
for n=1:parameter_data(1).camera_num;
    % This is the calibration pathname to read into MATLAB
    pathname_read=parameter_data(n).calibration_dir;
    % This is the calibration filename to read into MATLAB
    filename_read=dir([pathname_read,'*_calibration.mat']);
    % This loads the nth camera's point correspondence vectors which are of
    % the form [x,y]=F([X,Y,Z])
    load([pathname_read,filename_read.name]);
    % This checks the number of particles for calibration and downsamples
    % them if the particle number is to high
    if length(x)>point_max;
        % This generates a random index
        random_index=randperm(length(x));
        % This selects only the first point_max points
        random_index=random_index(1:point_max);
        % This extracts the random point_max points from the calibration
        % set
        x=x(random_index);
        y=y(random_index);
        X=X(random_index);
        Y=Y(random_index);
        Z=Z(random_index);
    end;
    % This creates a matrix w for inputing into the linear calibration
    % function
    w=ones(size(x));
    % This calculates the linear calibration matrix
    P=linear_calibration_04(x,y,w,X,Y,Z);
    % This adds the data into the calibration structure
    cal_struct(n).x_cal=P;
    cal_struct(n).y_cal=P;
    % This saves the calibration type into the structure
    cal_struct(n).type='linear';
    % This is the filename to save the calibration to
    filename_save=[pathname_read,'calibration_struct.mat'];
    % This save the calibration data to the drive
    save(filename_save,'cal_struct');
end;



function [x,y]=linear_calibration_eval(P,X,Y,Z);
% This function calculates the image coordinates of a set of vectors of
% world coordinates using the linear projective perspective matrix P.

% This initializes the output vectors
x=zeros(size(X));
y=zeros(size(X));
% This iterates through the list of points
for ii=1:length(X);
    cam_vect=P*[X(ii);Y(ii);Z(ii);1];
    x(ii)=cam_vect(1)/cam_vect(3);
    y(ii)=cam_vect(2)/cam_vect(3);
end;



function P=linear_calibration_04(x,y,w,X,Y,Z);
% This function is for generating the linear calibration matrix P given a
% set of homogeneous world coordinates [X,Y,Z,1] and homogeneous image
% coordinates [x,y,w].  A minimum of six world and image coordinate pairs
% is needed (actually 5.5 points as only 11 equations are needed to solve
% for P).  This algorithm is based on that described in "Multiple View
% Geometry in Computer Vision" by Richard Hartley and Andrew Zisserman in
% section 6.1
%
% This version is nearly identical to version 03 except that the null space
% of A is calculated using a singular value decomposition.  The command
% 'null' will not perform a least squares approximation while directly
% calculating the SVD allows for this ability.

% This extracts the length of the vectors (or the number of points)
n=length(x);
% This saves the inhomogeneous coordinates for some later fun
x_in=x;
y_in=y;
% This ensures that the image coordinates are homogenous
x=x./w;
y=y./w;
w=ones(size(w));
% This calculates the similarity transforms
[T,U]=similarity_transform(x,y,X,Y,Z);
% This initializes the transformed image coordinates
x_trans=zeros(size(x));
y_trans=zeros(size(y));
w_trans=zeros(size(w));
% This initializes the transformed world coordinates
X_Trans=zeros(size(X));
Y_Trans=zeros(size(Y));
Z_Trans=zeros(size(Z));
% This transforms the image and world coordinates
for ii=1:n;
    % This is the current image coordinate
    x_image=[x(ii);y(ii);1];
    % This is the current world coordinate
    X_World=[X(ii);Y(ii);Z(ii);1];
    % This transforms the image coordinates
    x_image=T*x_image;
    % This transforms the world coordinates
    X_World=U*X_World;
    % This saves the transformed image coordinates
    x_trans(ii)=x_image(1);
    y_trans(ii)=x_image(2);
    w_trans(ii)=x_image(3);
    % This saves the transformed world coordinates
    X_Trans(ii)=X_World(1);
    Y_Trans(ii)=X_World(2);
    Z_Trans(ii)=X_World(3);
end;
% This calculates the camera matrix in the transformed coordinates
P_Trans=camera_nullspace(x_trans,y_trans,w_trans,X_Trans,Y_Trans,Z_Trans);
% This converts the tranformed camera matrix back to the original
% coordinates sytem
P=inv(T)*P_Trans*U;
% This calculates the scaling factor of P
sigma=camera_scaling(x_in,y_in,X,Y,Z,P);
% This scales the camera matrix by the factor sigma
P=sigma*P;



function [T,U]=similarity_transform(x,y,X,Y,Z);
% This function yields similarity matricies that will scale the world and
% image coordinates so that the points centroid is at the origin and so
% that the image coordinates have an RMS distance of sqrt(2) from the
% origin and the world coordinates have an RMS distance of sqrt(3) from the
% origin.  This scaling is based on section 6.1 in "Multiple View Geometry 
% in Computer Vision" by Richard Hartley and Andrew Zisserman.

% This calculates the centroid of the image points
x_center=mean(x);
y_center=mean(y);
% This calculates the centroid of the world coordinates
X_Center=mean(X);
Y_Center=mean(Y);
Z_Center=mean(Z);
% This is the RMS distance of the image points from their centroid
d_vect=sqrt((x-x_center).^2+(y-y_center).^2);
d_rms=sqrt(mean(d_vect.^2));
% This is the RMS distance of the world points from their centroid
D_Vect=sqrt((X-X_Center).^2+(Y-Y_Center).^2+(Z-Z_Center).^2);
D_RMS=sqrt(mean(D_Vect.^2));
% These are the image translation and scale factors
s=sqrt(2)/d_rms;
tx=-x_center*s;
ty=-y_center*s;
% These are the world translation and scale factors
S=sqrt(3)/D_RMS;
TX=-X_Center*S;
TY=-Y_Center*S;
TZ=-Z_Center*S;
% This is the image similarity transform matrix
T=[s,0,tx;0,s,ty;0,0,1];
% This is the world similarity transform matrix
U=[S,0,0,TX;0,S,0,TY;0,0,S,TZ;0,0,0,1];



function P=camera_nullspace(x,y,w,X,Y,Z);
% This function calculates the camera calibration matrix from the
% homogeneous image and world coordinates.  This scaling is based on 
% section 6.1 in "Multiple View Geometry in Computer Vision" by Richard 
% Hartley and Andrew Zisserman.

% This extracts the length of the vectors (or the number of points)
n=length(x);
% This initializes the nullspace matrix A
A=zeros(2*n,12);
% This fills the matrix for each point
for ii=1:n;
    % This is the ii-th homogeneous world coordinate
    X_World=[X(ii),Y(ii),Z(ii),1];
    % This is the first equation for the current point
    A(2*(ii-1)+1,5:8)=-w(ii)*X_World;
    A(2*(ii-1)+1,9:12)=y(ii)*X_World;
    % This is the second equation for the current point
    A(2*ii,1:4)=w(ii)*X_World;
    A(2*ii,9:12)=-x(ii)*X_World;
end;
% This calculates the SVD of A
[U,S,V]=svd(A);
% The nullspace of A is given by any column vector in V whose corresponding
% singular value in S is zero; since the singular values in S are sorted in
% decreasing order, the best approximation to the null space of A is given
% by the last column vector of V
p=V(:,end);
% This initializes the matrix P
P=zeros(3,4);
% This is the first row of P
P(1,:)=p(1:4);
% This is the second row of P
P(2,:)=p(5:8);
% This is the third row of P
P(3,:)=p(9:12);



function sigma=camera_scaling(x_in,y_in,X,Y,Z,P);
% This function calculates the scaling factor of the camera matrix so the
% transformation from world to image coordinates is correctly scaled.  The
% function works by calculating the average ratio of the coordinates within
% a relative distance of mu fromt he median ratio.  This algorithm seams to
% work well, but it might go crazy if some of the entries of P are zero.


% This is the relative tolerance on estimating the camera matrix scale
% factor
mu=1e-2;
% This extracts the length of the vectors (or the number of points)
n=length(x_in);
% This initializes the calculated image coordinate vectors
x_calc=zeros(size(x_in));
y_calc=zeros(size(y_in));
% This calculates the image coordinates under the estimated camera matrix
% (in order to determine the scaling factor)
for ii=1:n;
    M_Calc=P*[X(ii);Y(ii);Z(ii);1];
    x_calc(ii)=M_Calc(1);
    y_calc(ii)=M_Calc(2);
end;
% This is the ratio of the scaled image coordinate system
qx=x_in./x_calc;
qy=y_in./y_calc;
% This concatenates the ratios into one matrix
q=[qx,qy];
% This ensures that q is a vector
q=q(:);
% This is the median of the ratios
m=median(q(:));
% This is a list of all components in q that are with a relative distance
% of mu*m away from the median m
median_list=logical((q>(m-mu*abs(m))).*(q<(m+mu*abs(m))));
% This is the average of all the components near the median
sigma=mean(q(median_list));



function [x_vect,y_vect,z_vect]=meshgrid_indexing_01(x,y,z,ii_vect);
% This function returns vectors that are equivalent to that given by the
% following commands
%
%   [X,Y,Z]=meshgrid(x,y,z);
%   X=X(:);
%   Y=Y(:);
%   Z=Z(:);
%   x_vect=X(ii_vect);
%   y_vect=Y(ii_vect);
%   z_vect=Z(ii_vect);
%
% except in a more memory efficient way so that for very large matricies
% out of memory errors do not occur.  The memory required to implement this
% function should be on the same order as that required to store the
% ii_vect vector.

% These are the lengths of the x and y vectors (the length of z is not
% required by the algorithm)
xn=length(x);
yn=length(y);
% These use some modular arithmetic to extract the indicies to which the
% x,y,z vectors are indexed into
ii_x=mod(floor((ii_vect-1)/yn),xn)+1;
ii_y=mod(ii_vect-1,yn)+1;
ii_z=ceil(ii_vect./(xn*yn));
% This extracts the indexed components of the vectors
x_vect=x(ii_x);
y_vect=y(ii_y);
z_vect=z(ii_z);



function M=interp_matrix_01(x,y,x_cam_res,y_cam_res);
% This function returns a sparse interpolation matrix M such that the value
% of the image at pixel coordinates [x(n),y(n)] is given by nth element of
% the product M*I(:) where I is the image matrix.

% This rounds the coordinates [x,y] to the nearest pixel
x_near=round(x);
y_near=round(y);
% This generates a list of interpolation points that are at the edge or
% beyond the edge of the image
x_near_lower=logical(x_near<=1);
x_near_upper=logical(x_cam_res<=x_near);
y_near_lower=logical(y_near<=1);
y_near_upper=logical(y_cam_res<=y_near);
% This is a list of any points at or beyond the image boundary
xy_boundary=or(or(x_near_lower,x_near_upper),or(y_near_lower,y_near_upper));
% Any interpolation point that is beyond the boundary of the image is set
% to the boundary
x_near(x_near_lower)=2;
x_near(x_near_upper)=x_cam_res-1;
y_near(y_near_lower)=2;
y_near(y_near_upper)=y_cam_res-1;
% This gives the fractional part of the coordinate [x,y]
r=x-x_near;
s=y-y_near;
% This computes the interpolation coefficients for each coordinate
h=nine_node_interp_01(r,s);
% This sets the interpolation coefficients of any points near or beyond the
% boundary of the image to zero
h(:,xy_boundary)=0;
% This initializes the linear index matrix
lin_indx_h=zeros(9,size(x,1));
% These are the linear indicies of each of the nearest pixels
lin_indx_h(1,:)=sub2ind([y_cam_res,x_cam_res],y_near+1,x_near+1)';
lin_indx_h(2,:)=sub2ind([y_cam_res,x_cam_res],y_near+1,x_near-1)';
lin_indx_h(3,:)=sub2ind([y_cam_res,x_cam_res],y_near-1,x_near-1)';
lin_indx_h(4,:)=sub2ind([y_cam_res,x_cam_res],y_near-1,x_near+1)';
lin_indx_h(5,:)=sub2ind([y_cam_res,x_cam_res],y_near+1,x_near+0)';
lin_indx_h(6,:)=sub2ind([y_cam_res,x_cam_res],y_near+0,x_near-1)';
lin_indx_h(7,:)=sub2ind([y_cam_res,x_cam_res],y_near-1,x_near+0)';
lin_indx_h(8,:)=sub2ind([y_cam_res,x_cam_res],y_near+0,x_near+1)';
lin_indx_h(9,:)=sub2ind([y_cam_res,x_cam_res],y_near+0,x_near+0)';
% This is the number of points to interpolate
interp_points=size(x,1);
% This is the vector of ii positions
ii_vect=repmat(1:interp_points,9,1);
ii_vect=ii_vect(:);
% This is the vector of jj positions
jj_vect=lin_indx_h(:);
% These are the values of the matrix at [ii,jj]
ss_vect=h(:);
% This creates the spares interpolation matrix
M=sparse(ii_vect,jj_vect,ss_vect,interp_points,x_cam_res*y_cam_res);



function h=nine_node_interp_01(r,s);
% This function creates interpolation coefficients that are functions of
% the nearest nine pixels.  In this way the interpolation can be expressed
% as a linear sum of pixel intensities.

% This initializes the interpolation coefficient matrix
h=zeros(9,size(r,1));
% This calculates the interpolation coefficients for each adjacent pixel
h(1,:) = ( (1/4)*(1+r).*(1+s).*r.*s)';
h(2,:) = (-(1/4)*(1-r).*(1+s).*r.*s)';
h(3,:) = ( (1/4)*(1-r).*(1-s).*r.*s)';
h(4,:) = (-(1/4)*(1+r).*(1-s).*r.*s)';
h(5,:) = ( (1/2)*(1-r).*(1+r).*s.*(1+s))';
h(6,:) = (-(1/2)*(1-s).*(1+s).*r.*(1-r))';
h(7,:) = (-(1/2)*(1-r).*(1+r).*s.*(1-s))';
h(8,:) = ( (1/2)*(1+s).*(1-s).*r.*(1+r))';
h(9,:) = ( (1+s).*(1-s).*(1+r).*(1-r))';