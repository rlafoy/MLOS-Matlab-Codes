function particle_id_3d_01;
% This function is designed to identify and size the particles in a 3D
% intensity field.

% This is a test data set
filename_read='/mnt/current_storage/Projects2/Tomo_PIV/Simulated_Vortex_Ring/Tomo_Recon_Method_Tests/111031_Error_Analysis/PIV_Stokes_Vortex_Ring/Vortex_Ring_004662_Particles/MART_Reconstruction_Best/S00001.mat';

% This loads the 3D image I
load(filename_read);

% % These are the particle peaks in the image
% bw_peaks=imregionalmax(I,18);
% % This is a labeled array of the particle peaks
% bw_peaks_label_struct=bwconncomp(bw_peaks);
% bw_peak_label=bw_peaks_label_struct.PixelIdxList;
% bw_peak_num=bw_peaks_label_struct.NumObjects;
% % This calculates the mean particle intensity
% I_Peak_Mean=mean(I(bw_peaks(:)));

% This computes the watershed of the image
I_WS=watershed(I,18);

% This is a vector of the numerical linear (as opposed to logical) indices
nonlogical_linear_indices=1:numel(I);

I_Double=double(I);

% This iterates through the different watersheds fitting Quassian functions
% to the peaks
for watershed_number=2:max(I_WS(:));
    
    % These are the indices of the current watershed
    watershed_indices=(I_WS(:)==watershed_number);
    
    % These are the numerical linear (non-logical) indices of the current
    % watershed
    linear_vector=nonlogical_linear_indices(watershed_indices);
    
    % This is the maximum intensity of the current watershed
    [max_watershed,max_index]=max(I(linear_vector));
    
    % This converts the maximum index to subscript indices
    [max_ii,max_jj,max_kk]=ind2sub(size(I),linear_vector(max_index));
    
    % This is the initial estimate of the Gaussian peak parameters
    parameter_vector_zero=[max_ii,max_jj,max_kk,1,1,1];
    
    % This is the current list of indices corresponding to the current
    % watershed
    [ii_vector,jj_vector,kk_vector]=ind2sub(size(I),linear_vector);
    
    % This fits the current watershed to a Gaussian function
    parameter_vector=fminsearch(@(parameter_vector)gaussian_residual(parameter_vector,I_Double,ii_vector,jj_vector,kk_vector,linear_vector,max_watershed),parameter_vector_zero);
    
    keyboard;
    
end;
    
    

% % This is a list of random indicies to the peaks
% rand_indicies=randperm(bw_peak_num);
% 
% % This is a counting index
% count=0;
% % This calculates the peak diameter for a range of particles
% for n=1:particle_num;
%     
%     % These are the linear indicies of the current particle peak
%     peak_indicies=bw_peak_label{rand_indicies(n)};
%     % This takes the first index
%     peak_index=peak_indicies(1);
%     
%     % This is the location of the current particle in indicies
%     [ii_prt,jj_prt,kk_prt]=ind2sub(size(I),peak_index);
%     
%     try;
%         % This extracts a region around the peak
%         I_Win=I(ii_prt-ds:ii_prt+ds,jj_prt-ds:jj_prt+ds,kk_prt-ds:kk_prt+ds);
%         count=count+1;
%     catch;
%         continue;
%     end;
%     
%     % This is an approximation to the particle standard deviation size
%     %sigma(count)=sum(I_Win(:))/(I_Peak_Mean*sqrt(2*pi));
%     sigma(count)=(((3*sum(I_Win(:)>(I_Peak_Mean/2)))/(4*pi))^(1/3))/sqrt(2*log(2));
%     
% end;
% % This is the mean of the particle standard deviations
% I_Peak_Sigma=mean(sigma);



function r=gaussian_residual(parameter_vector,I,ii_vector,jj_vector,kk_vector,linear_vector,max_watershed);
% This returns the residual of the Gaussian fit to the actual particle
% image.

% % This is the peak intensity of the Gaussian function
% A=parameter_vector(1);

A=max_watershed;

% This is the ii location of the peak
ii_max=parameter_vector(1);
% This is the jj location of the peak
jj_max=parameter_vector(2);
% This is the kk_location of the peak
kk_max=parameter_vector(3);
% This is the ii standard deviation
ii_sigma=parameter_vector(4);
% This is the jj standard deviation
jj_sigma=parameter_vector(5);
% This is the kk standard deviation
kk_sigma=parameter_vector(6);

% This is the analytical Gaussian function
G=A*exp(-((ii_vector-ii_max).^2/(2*ii_sigma^2)+(jj_vector-jj_max).^2/(2*jj_sigma^2)+(kk_vector-kk_max).^2/(2*kk_sigma^2)));

% This is the sum of the squares of the errors to the actual image
r=sum((G-I(linear_vector)).^2);


