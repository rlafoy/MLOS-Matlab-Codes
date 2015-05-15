function [ghost_statistics_concise,ghost_statistics_verbose]=ghost_particle_statistics_02(testing_parameters,parameter_data);
% This function reads in both the true intensity field and the
% reconstructed intensity field and identifies all "particles" in both
% images using a regional maximum criterion.  A gaussian profile is fit to
% the particles to extract the subpixel positions, particle diameter, and
% particle intensity.  The true particles are differentiated from the ghost
% particles from the calibration data containing all the particle
% positions.  The particle statistics are then saved.
%
% This version is similar to version 01 except that it allows for the data
% locations to be passed in as a parameter structure.

% This is the directory containging the true 3D intensity fields
intensity_filename=dir([testing_parameters.intensity_directory,'*.mat']);
intensity_filename=[testing_parameters.intensity_directory,intensity_filename(1).name];
% This is the filename containing the tomographic reconstruction
reconstruction_filename=dir([parameter_data(1).full_reconstruction_dir,'*.mat']);
reconstruction_filename=[parameter_data(1).full_reconstruction_dir,reconstruction_filename(1).name];
% This is the filename containing the particle locations
particles_filename=dir([parameter_data(1).calibration_dir,'cam_*.mat']);
particles_filename=[parameter_data(1).calibration_dir,particles_filename(1).name];

% This loads the particle positions vectors
load(particles_filename);

% This displays the the true intensity field is being loaded
disp('Loading the true intensity field . . . ');
% This loads the true 3D intensity field
load(intensity_filename);
% This displays that the particles in the true intensity field are being
% identified
disp('Indentifying particles in the true intensity field . . . ');
% This extracts the gaussian parameters describing the particles in the
% true intensity field (for comparison to the reconstructed field)
[C_True,mu_true,sigma_true]=quadratic_particles(I);
% This deletes the intensity field from memory as it may be quite large
clear('I');

% These are the particle coordinates in the reconstructed volume
X_Mu_True=mu_true(:,2);
Y_Mu_True=mu_true(:,1);
Z_Mu_True=mu_true(:,3);
% This converts from voxel coordinates to world coordinates
[X_Mu,Y_Mu,Z_Mu]=voxel2world(X_Mu_True,Y_Mu_True,Z_Mu_True,testing_parameters);
% These are the indices corresponding to the true and ghost particles
[index_true_true,index_true_ghost,ii_true_match]=particle_matching(X,Y,Z,X_Mu,Y_Mu,Z_Mu,sigma_true);
% These are the number of actual and ghost particles in the true field
true_true_particle_num=length(index_true_true);
true_ghost_particle_num=length(index_true_ghost);
% These are the average intensities of the true and ghost particles in the
% true field
c_mean_true_true=mean(C_True(index_true_true,:),1);
c_mean_true_ghost=mean(C_True(index_true_ghost,:),1);
% These are the standard deviations of the intensities of the true and 
% ghost particles in the true field
c_std_true_true=std(C_True(index_true_true,:),0,1);
c_std_true_ghost=std(C_True(index_true_ghost,:),0,1);
% These are the average diameters of the true and ghost particles in the
% true field
sigma_mean_true_true=mean(sigma_true(index_true_true,:),1);
sigma_mean_true_ghost=mean(sigma_true(index_true_ghost,:),1);
% These are the standard deviations of the diameters of the true and 
% ghost particles in the true field
sigma_std_true_true=std(sigma_true(index_true_true,:),0,1);
sigma_std_true_ghost=std(sigma_true(index_true_ghost,:),0,1);
% These are the particle matches in the true intensity field
X_Mu_True_Match=X_Mu_True(ii_true_match);
Y_Mu_True_Match=Y_Mu_True(ii_true_match);
Z_Mu_True_Match=Z_Mu_True(ii_true_match);

% This displays the the reconstructed intensity field is being loaded
disp('Loading the reconstructed intensity field . . . ');
% This loads the reconstructed 3D intensity field
load(reconstruction_filename);
% This displays that the particles in the reconstructed intensity field 
% are being identified
disp('Indentifying particles in the reconstructed intensity field . . . ');
% This extracts the gaussian parameters describing the particles in the
% reconstructed intensity field
[C_Recon,mu_recon,sigma_recon]=quadratic_particles(I);

% These are the particle coordinates in the reconstructed volume
X_Mu_Recon=mu_recon(:,2);
Y_Mu_Recon=mu_recon(:,1);
Z_Mu_Recon=mu_recon(:,3);
% This converts from voxel coordinates to world coordinates
[X_Mu,Y_Mu,Z_Mu]=voxel2world(X_Mu_Recon,Y_Mu_Recon,Z_Mu_Recon,testing_parameters);
% These are the indices corresponding to the true and ghost particles
[index_recon_true,index_recon_ghost,ii_recon_match]=particle_matching(X,Y,Z,X_Mu,Y_Mu,Z_Mu,sigma_recon);
% These are the number of actual and ghost particles in the true field
recon_true_particle_num=length(index_recon_true);
recon_ghost_particle_num=length(index_recon_ghost);
% These are the average intensities of the true and ghost particles in the
% true field
c_mean_recon_true=mean(C_Recon(index_recon_true,:),1);
c_mean_recon_ghost=mean(C_Recon(index_recon_ghost,:),1);
% These are the standard deviations of the intensities of the true and 
% ghost particles in the true field
c_std_recon_true=std(C_Recon(index_recon_true,:),0,1);
c_std_recon_ghost=std(C_Recon(index_recon_ghost,:),0,1);
% These are the average diameters of the true and ghost particles in the
% true field
sigma_mean_recon_true=mean(sigma_recon(index_recon_true,:),1);
sigma_mean_recon_ghost=mean(sigma_recon(index_recon_ghost,:),1);
% These are the standard deviations of the diameters of the true and 
% ghost particles in the true field
sigma_std_recon_true=std(sigma_recon(index_recon_true,:),0,1);
sigma_std_recon_ghost=std(sigma_recon(index_recon_ghost,:),0,1);
% These are the particle matches in the true intensity field
X_Mu_Recon_Match=X_Mu_Recon(ii_recon_match);
Y_Mu_Recon_Match=Y_Mu_Recon(ii_recon_match);
Z_Mu_Recon_Match=Z_Mu_Recon(ii_recon_match);

% This is the probability of guessing correctly whether a particle is true
% or ghost based upon its intensity
%
% These are the intensities of the true and ghost particles
c_recon_true=mean(C_Recon(index_recon_true,:),2);
c_recon_ghost=mean(C_Recon(index_recon_ghost,:),2);
% These are the histograms of the intensities of the true and ghost
% particles
[n_true,x_true]=hist(c_recon_true,0.01:0.02:0.99);
[n_ghost,x_ghost]=hist(c_recon_ghost,0.01:0.02:0.99);
% This is the probability of guessing whether a particle is a true or ghost
% particle
P=sum(max(n_true,n_ghost))/(sum(n_true+n_ghost));

disp(['The estimation probability is ',num2str(P),'.']);

% This is the percentage overlap area compared to the total are of the true
% particles (this number should be as small as possible, if it is greater
% than one then the ghost particles and the true particles are
% indistinguishable)
Overlap_Percent=sum(min(n_true,n_ghost))/sum(n_true);

disp(['The percent overlap is ',num2str(Overlap_Percent)]);

% These are the particle position errors
X_Error=X_Mu_True_Match-X_Mu_Recon_Match;
Y_Error=Y_Mu_True_Match-Y_Mu_Recon_Match;
Z_Error=Z_Mu_True_Match-Z_Mu_Recon_Match;
% These are the average particle errors
X_Error_Mean=mean(X_Error);
Y_Error_Mean=mean(Y_Error);
Z_Error_Mean=mean(Z_Error);
% These are the error standard deviations
X_Error_STD=std(X_Error);
Y_Error_STD=std(Y_Error);
Z_Error_STD=std(Z_Error);

% This is the intensity of the true particle field
C_f=mean(c_mean_true_true);
% This is the standard deviation vector of the true particle field
Sigma_f=sigma_mean_true_true;

% This is the intensity of the true particles in the reconstructed
% particle field
C_g=mean(C_Recon(index_recon_true,:),2);
% This is the standard deviation vector of the true particles in the
% reconstructed particle field
Sigma_g=sigma_recon(index_recon_true,:);
% These are the correllations of the true particles
Corr_True=guassian_correlation(C_f,C_g,Sigma_f,Sigma_g);

% This is the intensity of the true particles in the reconstructed
% particle field
C_g=mean(C_Recon(index_recon_ghost,:),2);
% This is the standard deviation vector of the true particles in the
% reconstructed particle field
Sigma_g=sigma_recon(index_recon_ghost,:);
% These are the correllations of the true particles
Corr_Ghost=guassian_correlation(C_f,C_g,Sigma_f,Sigma_g);

% This initializes the structure for saving only some of the statistics
% from the ghost particle analysis
ghost_statistics_concise=struct;
% These are the number of actual and ghost particles identified in the true
% intensity field
ghost_statistics_concise.true_actual_particle_num=true_true_particle_num;
ghost_statistics_concise.true_ghost_particle_num=true_ghost_particle_num;
% These are the average intensities of the true and ghost particles in the 
% true intensity field
ghost_statistics_concise.true_actual_intensity_mean=c_mean_true_true;
ghost_statistics_concise.true_ghost_intensity_mean=c_mean_true_ghost;
% These are the standard deviations of the particle intensities in the true
% intensity field
ghost_statistics_concise.true_actual_intensity_std=c_std_true_true;
ghost_statistics_concise.true_ghost_intensity_std=c_std_true_ghost;
% These are the average particle diameters of the true intensity field (as
% defined by the standard deviation of the guassian fitting function)
ghost_statistics_concise.true_actual_diameter_mean=sigma_mean_true_true;
ghost_statistics_concise.true_ghost_diameter_mean=sigma_mean_true_ghost;
% These are the standard deviations of the particle diameters (as
% defined by the standard deviation of the guassian fitting function)
ghost_statistics_concise.true_actual_diameter_std=sigma_std_true_true;
ghost_statistics_concise.true_ghost_diameter_std=sigma_std_true_ghost;
% These are the number of actual and ghost particles identified in the
% reconstructed intensity field
ghost_statistics_concise.recon_actual_particle_num=recon_true_particle_num;
ghost_statistics_concise.recon_ghost_particle_num=recon_ghost_particle_num;
% These are the average intensities of the true and ghost particles in the 
% reconstructed intensity field
ghost_statistics_concise.recon_actual_intensity_mean=c_mean_recon_true;
ghost_statistics_concise.recon_ghost_intensity_mean=c_mean_recon_ghost;
% These are the standard deviations of the particle intensities in the
% reconstructed intensity field
ghost_statistics_concise.recon_actual_intensity_std=c_std_recon_true;
ghost_statistics_concise.recon_ghost_intensity_std=c_std_recon_ghost;
% These are the average particle diameters of the reconstructed intensity 
% field (as defined by the standard deviation of the guassian fitting 
% function)
ghost_statistics_concise.recon_actual_diameter_mean=sigma_mean_recon_true;
ghost_statistics_concise.recon_ghost_diameter_mean=sigma_mean_recon_ghost;
% These are the standard deviations of the particle diameters (as
% defined by the standard deviation of the guassian fitting function)
ghost_statistics_concise.recon_actual_diameter_std=sigma_std_recon_true;
ghost_statistics_concise.recon_ghost_diameter_std=sigma_std_recon_ghost;
% These are the mean correllation values of the actual and ghost particles
ghost_statistics_concise.recon_actual_corr_mean=mean(Corr_True);
ghost_statistics_concise.recon_ghost_corr_mean=mean(Corr_Ghost);
% This is the ratio of the sum of the correllations
ghost_statistics_concise.recon_corr_sum_ratio=sum(Corr_True)/sum(Corr_Ghost);
% These are the errors in the measured particle positions
ghost_statistics_concise.x_particle_position_error=X_Error;
ghost_statistics_concise.y_particle_position_error=Y_Error;
ghost_statistics_concise.z_particle_position_error=Z_Error;
% This is the probability of guessing whether the particle is true or a
% ghost based on intensity
ghost_statistics_concise.recon_particle_probability=P;
% This is the percentage overlap area of the ghost particles with the true
% particles divided by the total area of the true particles (this number
% should be as low as possible)
ghost_statistics_concise.recon_intensity_overlap_ratio=Overlap_Percent;


% This initializes the structure for saving all the relevant data
ghost_statistics_verbose=struct;
% These are the parameters for the true intensity field particles
ghost_statistics_verbose.true_intensity=C_True;
ghost_statistics_verbose.true_position=mu_true;
ghost_statistics_verbose.true_diameter=sigma_true;
% These are the parameters for the reconstructed intensity field particles
ghost_statistics_verbose.recon_intensity=C_Recon;
ghost_statistics_verbose.recon_position=mu_recon;
ghost_statistics_verbose.recon_diameter=sigma_recon;
% These are the indices of the actual and ghost particles in the
% true intensity image (with particles with NaN diameters ignored)
% These are the parameters for the true intensity field
ghost_statistics_verbose.true_actual_index=index_true_true;
ghost_statistics_verbose.true_ghost_index=index_true_ghost;
% These are the indices of the actual and ghost particles in the
% reconstructed image (with particles with NaN diameters ignored)
% These are the parameters for the true intensity field
ghost_statistics_verbose.recon_actual_index=index_recon_true;
ghost_statistics_verbose.recon_ghost_index=index_recon_ghost;
% These are the correlation magnitudes of the actual and ghost particles
% (indexed by recon_index_true and recon_index_ghost)
ghost_statistics_verbose.recon_corr_actual=Corr_True;
ghost_statistics_verbose.recon_corr_ghost=Corr_Ghost;
% These are the mean correllation values of the actual and ghost particles
ghost_statistics_verbose.recon_actual_corr_mean=mean(Corr_True);
ghost_statistics_verbose.recon_ghost_corr_mean=mean(Corr_Ghost);
% This is the ratio of the sum of the correllations
ghost_statistics_verbose.recon_corr_sum_ratio=sum(Corr_True)/sum(Corr_Ghost);
% These are the errors in the measured particle positions
ghost_statistics_verbose.x_particle_position_error=X_Error;
ghost_statistics_verbose.y_particle_position_error=Y_Error;
ghost_statistics_verbose.z_particle_position_error=Z_Error;
% This is the probability of guessing whether the particle is true or a
% ghost based on intensity
ghost_statistics_verbose.recon_particle_probability=P;
% This is the percentage overlap area of the ghost particles with the true
% particles divided by the total area of the true particles (this number
% should be as low as possible)
ghost_statistics_concise.recon_intensity_overlap_ratio=Overlap_Percent;

% % This displays the correlation statistics
% disp(['Mean True Correlation: ',num2str(mean(Corr_True))]);
% disp(['Mean Ghost Correlation: ',num2str(mean(Corr_Ghost))]);
% disp(['Sum of True Correlation / Sum of Ghost Correlation: ',num2str(sum(Corr_True)/sum(Corr_Ghost))]);
% 
% figure(1);
% [n_true,x_true]=hist(mean(C_Recon(index_recon_true,:),2),50);
% [n_ghost,x_ghost]=hist(mean(C_Recon(index_recon_ghost,:),2),50);
% bar(x_true,n_true);
% hold on;
% bar(x_ghost,n_ghost);
% hold off;
% h = findobj(gca,'Type','patch');
% set(h(1),'FaceColor','BLUE','EdgeColor','WHITE');
% set(h(2),'FaceColor','RED','EdgeColor','WHITE');
% axis([0,1,0,500]);
% 
% figure(4);
% [n_true,x_true]=hist(mean(sigma_recon(index_recon_ghost,:),2),1000);
% [n_ghost,x_ghost]=hist(mean(sigma_recon(index_recon_true,:),2),25);
% bar(x_true,log10(n_true));
% hold on;
% bar(x_ghost,log10(n_ghost));
% hold off;
% h = findobj(gca,'Type','patch');
% set(h(1),'FaceColor','RED','EdgeColor','WHITE');
% set(h(2),'FaceColor','BLUE','EdgeColor','WHITE');
% axis([0.5,5,1,5]);

% 
% figure(3);
% plot(X_Error,Y_Error,'.BLACK');
% axis equal;
% axis([-1.0,1.0,-1.0,1.0]);
% grid on;
% xlabel('X Error');
% ylabel('Y Error');
% 
% figure(4);
% plot(X_Error,Z_Error,'.BLACK');
% axis equal;
% axis([-1.2,1.2,-1.2,1.2]);
% grid on;
% xlabel('X Error');
% ylabel('Z Error');
% 
% keyboard;



function [index_true,index_ghost,ii_match]=particle_matching(X,Y,Z,X_Mu,Y_Mu,Z_Mu,sigma);
% This function returns the indices of the particles that match the
% positions of the true particles X,Y,Z the closest (in world coordinates).
% Particles with NaN of Inf values are excluded.

% This displays that the particles are being matched
disp('Performing particle matching process . . . ');
% This finds the nearest neighbor matches between the known particle
% locations (X,Y,Z) and the identified particle locations (X_Mu,Y_Mu,Z_Mu)
[~,ii_match]=nearest_neighbor_matching(X,Y,Z,X_Mu,Y_Mu,Z_Mu);
% These are indices of the true particles
index_true=ii_match;
% These are the indices of the ghost particles
index_ghost=setdiff(1:length(X_Mu),ii_match)';

% Since particles on the edge of the image do not have well defined
% diameters, the diameter of these particles is defined as NaN.  This finds
% any of these particles and removes them from the index list (as the NaN 
% will not provide useful statistics)
index_true(logical(sum(isnan(sigma(index_true,:))|isinf(sigma(index_true,:)),2)))=[];
index_ghost(logical(sum(isnan(sigma(index_ghost,:))|isinf(sigma(index_ghost,:)),2)))=[];



function Corr=guassian_correlation(C_f,C_g,Sigma_f,Sigma_g);
% This function determines the correllation of the two gaussian functions
% determined by
%
%  f(x,y,z) = C_f*exp(-(x^2)/(2*Sigma_f(1)^2)-(y^2)/(2*Sigma_f(2)^2)-(x^2)/(2*Sigma_f(3)^2))
%  g(x,y,z) = C_g*exp(-(x^2)/(2*Sigma_g(1)^2)-(y^2)/(2*Sigma_g(2)^2)-(x^2)/(2*Sigma_g(3)^2))
%
% where C_f and C_g are the amplitudes of the guassian functions and
% Sigma_f and Sigma_g and vectors containing the gaussian standard
% deviations.  C_f must be a scalar and Sigma_f must be 1 x 3.  C_g may be
% N x 1 and Sigma_g may be N x 3.

% This is the numerator of the correlation
Corr_Num=C_f*C_g*(2*pi)^(3/2);
% This is the denominator of the correlation
Corr_Den=sqrt((1/Sigma_f(1,1)^2+1./Sigma_g(:,1)).*(1/Sigma_f(1,2)^2+1./Sigma_g(:,2)).*(1/Sigma_f(1,3)^2+1./Sigma_g(:,3)));
% This is the total correllation
Corr=Corr_Num./Corr_Den;



function [d_match,ii_match]=nearest_neighbor_matching(X1,Y1,Z1,X2,Y2,Z2);
% This function finds the nearest neighbor match for all of X1,Y1,Z1 in the
% set of particles X2,Y2,Z2 and returns the corresponding index of X2,Y2,Z2
% and the distance to the nearest neighbor.

% This initializes the nearest particle index vector
ii_match=zeros(size(X1));
% This initializes the nearest particle distance vector
d_match=zeros(size(X1));
% This iterates through the known true particles (given by X1,Y1,Z1)
% finding the closest match in the identified particles (given by X2,Y2,Z2)
for n=1:length(X1);
    % This is the squared distance to all the identified particles
    d=(X2-X1(n)).^2+(Y2-Y1(n)).^2+(Z2-Z1(n)).^2;
    % This finds the smallest distance particle
    [d_min,min_index]=min(d);
    % This saves the minimum distance index
    ii_match(n)=min_index;
    % This saves the minimum distance magnitude
    d_match(n)=sqrt(d_min);
    % This is a short pause for breaking purposes
    pause(0.001);
end;



function [X_World,Y_World,Z_World]=voxel2world(X_Voxel,Y_Voxel,Z_Voxel,parameter_data);
% This function uses the data in the parameter_data structure to convert
% the voxel coordinates given by the vectors X_Voxel, Y_Voxel, Z_Voxel to
% world coordinates.

% This is the volume domain over which the reconstruction was performed
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
% This converts the particle voxel coordinates to world coordinates
X_World=(xmax-xmin)*(X_Voxel-1)/(xres_volume-1)+xmin;
Y_World=(ymax-ymin)*(Y_Voxel-1)/(yres_volume-1)+ymin;
Z_World=(zmax-zmin)*(Z_Voxel-1)/(zres_volume-1)+zmin;



function [X_Voxel,Y_Voxel,Z_Voxel]=world2voxel(X_World,Y_World,Z_World,parameter_data);
% This function converts from the world coordinates given by the vectors
% X_World, Y_World, Z_World to voxel coordinates given by the vectors
% X_Voxel, Y_Voxel, Z_Voxel using the parameters defined in parameter_data.

% This is the volume domain over which the reconstruction was performed
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
% This converts the world coorindates to voxel coordiantes
X_Voxel=(X_World-xmin)*(xres_volume-1)/(xmax-xmin)+1;
Y_Voxel=(Y_World-ymin)*(yres_volume-1)/(ymax-ymin)+1;
Z_Voxel=(Z_World-zmin)*(zres_volume-1)/(zmax-zmin)+1;



function [C,mu,sigma]=quadratic_particles(I);
% This function identifies all local maxima in the image I and fits a
% quadratic approximation to a gaussian function to the intensity values to
% find the guassian parameters.  The fit is performed in each of the i,j,k
% dimensions seperately.  At the edges of the image I the fit is not
% performed and only the maxima intensity is returned.  The output
% parameters C, mu, sigma correspond to the maxima intensity, coordinates,
% and standard deviation (of the gaussian function) respectively.  The
% output parameters are N x 3 in size where N is the number of regional
% maxima and the columns correspond to the i,j,k dimensions respectively.

% This extracts the regional maxima of the image
if any(isnan(I(:)));
    I(isnan(I(:)))=0;
    warning('3D image intensity field contains NaN values.  Setting all NaN values to zero and continuing.');
end;
region_maxima=imregionalmax(I);
% This initializes the output matrices
C=zeros(sum(region_maxima(:)),3);
mu=zeros(sum(region_maxima(:)),3);
sigma=zeros(sum(region_maxima(:)),3);
% These are the indices of the regional maxima
lin_center=find(region_maxima);
% This converts from linear indices to subscript indices
[ii_center,jj_center,kk_center]=ind2sub(size(I),lin_center);

% These are the indices of the maxima that lie along the edge of the volume
% (for which the quadratic fitting will not work)
edge_index=false(size(ii_center));
edge_index=edge_index|(ii_center==1)|(ii_center==size(I,1));
edge_index=edge_index|(jj_center==1)|(jj_center==size(I,2));
edge_index=edge_index|(kk_center==1)|(kk_center==size(I,3));
% This defines the fitting parameters for the maxima along the edge of the
% image (these should be replaced by something more advanced such as a
% forward or backward fitting scheme . . . )
C(edge_index,:)=repmat(I(lin_center(edge_index)),1,3);
mu(edge_index,:)=[ii_center(edge_index),jj_center(edge_index),kk_center(edge_index)];
sigma(edge_index,:)=NaN;
% This deletes the maxima that are on the edge from the vectors of maxima
% indices (since an out of range error will occur if the adjacent points
% are attempted to be accessed)
ii_center(edge_index)=[];
jj_center(edge_index)=[];
kk_center(edge_index)=[];

% These are the vectors of indices in the ii direction
ii_vector=[ii_center-1,ii_center,ii_center+1];
jj_vector=repmat(jj_center,1,3);
kk_vector=repmat(kk_center,1,3);
% These are the linear indices into the image matrix in the ii direction
lin_index=sub2ind(size(I),ii_vector,jj_vector,kk_vector);
% These are the intensity values in the ii direction
I_Vector=I(lin_index);
% This calculates the guassian fitting parameters in the ii direction
[C_ii,mu_ii,sigma_ii]=quadratic_fitting(I_Vector);

% These are the vectors of indices in the jj direction
ii_vector=repmat(ii_center,1,3);
jj_vector=[jj_center-1,jj_center,jj_center+1];
kk_vector=repmat(kk_center,1,3);
% These are the linear indices into the image matrix in the jj direction
lin_index=sub2ind(size(I),ii_vector,jj_vector,kk_vector);
% These are t?he intensity values in the jj direction
I_Vector=I(lin_index);
% This calculates the guassian fitting parameters in the jj direction
[C_jj,mu_jj,sigma_jj]=quadratic_fitting(I_Vector);

% These are the vectors of indices in the kk direction
ii_vector=repmat(ii_center,1,3);
jj_vector=repmat(jj_center,1,3);
kk_vector=[kk_center-1,kk_center,kk_center+1];
% These are the linear indices into the image matrix in the kk direction
lin_index=sub2ind(size(I),ii_vector,jj_vector,kk_vector);
% These are the intensity values in the kk direction
I_Vector=I(lin_index);
% This calculates the guassian fitting parameters in the kk direction
[C_kk,mu_kk,sigma_kk]=quadratic_fitting(I_Vector);

% This saves the fitting parameters to the output matrices
C(not(edge_index),:)=[C_ii,C_jj,C_kk];
mu(not(edge_index),:)=[mu_ii+ii_center,mu_jj+jj_center,mu_kk+kk_center];
sigma(not(edge_index),:)=[sigma_ii,sigma_jj,sigma_kk];



function [C,mu,sigma]=quadratic_fitting(I_Vect);
% This function fits a parabolic approximation of a gaussian function to
% the N x 3 vector I_Vect where each row of the vector corresponds to a
% fitting in one dimension.

% These are the coefficients of the quadratic polynomial for the fitting of
% in the ii direction
c0=I_Vect(:,2);
c1=0.5*(I_Vect(:,3)-I_Vect(:,1));
c2=0.5*(I_Vect(:,1)-2*I_Vect(:,2)+I_Vect(:,3));
% This initializes a constant for calculating the parameters
A=c1.^2-2*c0.*c2;
% These are the guassian parameters in the ii direction
C=c0.*exp((c1.^2)./(2*A));
mu=(c0.*c1)./A;
sigma=c0./sqrt(A);