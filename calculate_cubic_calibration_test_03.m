function calculate_cubic_calibration_test_02;
% This function is for the development and testing of a cubic camera
% calibration function.

% This is a calibration data set for testing purposes
filename_read='/mnt/current_storage/Projects/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Calibration_01/cam_4/cam_4_calibration.mat';

% This is the filename to save the data
filename_save='/mnt/current_storage/Projects/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Calibration_01/cam_1/cam_1_calibration_function.mat';


% This loads the data file
load(filename_read);

% This eliminates Z planes that are too close together
[x,y,X,Y,Z]=eliminate_close_z_planes(x,y,X,Y,Z);

% This groups the data into cell array based upon the Z plane, ie all cell
% arrays of the same index are on the same Z plane
[x_plane,y_plane,X_plane,Y_plane,Z_plane]=group_z_planes(x,y,X,Y,Z);

index_vect=1:2:length(x_plane);

for ii=1:length(index_vect);
    x_plane_new{ii}=x_plane{index_vect(ii)};
    y_plane_new{ii}=y_plane{index_vect(ii)};
    X_plane_new{ii}=X_plane{index_vect(ii)};
    Y_plane_new{ii}=Y_plane{index_vect(ii)};
    Z_plane_new{ii}=Z_plane{index_vect(ii)};
end;

x_plane=x_plane_new;
y_plane=y_plane_new;
X_plane=X_plane_new;
Y_plane=Y_plane_new;
Z_plane=Z_plane_new;

x=[];
y=[];
X=[];
Y=[];
Z=[];
for ii=1:length(x_plane);
    x=[x;x_plane{ii}];
    y=[y;y_plane{ii}];
    X=[X;X_plane{ii}];
    Y=[Y;Y_plane{ii}];
    Z=[Z;Z_plane{ii}];
end;

% This is the filename of the new calibration data to write that excludes
% the in between grid points which cause all sorts of nasty problems (at
% least for the moment)
filename_write='/mnt/current_storage/Projects/Tomo_PIV/101107_Vortex_Ring_Processing/Camera_Position_01/Calibration_01/cam_4/cam_4_better_calibration.mat';

save(filename_write,'x','y','X','Y','Z');
return;

% This eliminates Z planes that are too close together
[x,y,X,Y,Z]=eliminate_close_z_planes(x,y,X,Y,Z);

% This groups the data into cell array based upon the Z plane, ie all cell
% arrays of the same index are on the same Z plane
[x_plane,y_plane,X_plane,Y_plane,Z_plane]=group_z_planes(x,y,X,Y,Z);

% This calculates the calibration function and saves the function to the
% calibration structure
calibration_function=calculate_calibration_function(x_plane,y_plane,X_plane,Y_plane,Z_plane);
% 
% save(filename_save,'calibration_function');

% load(filename_save);

n=1000;
X_Eval=zeros(n,1);
Y_Eval=zeros(n,1);
Z_Eval=linspace(unique(Z_plane{1}),unique(Z_plane{end}),n);


% % X=10*randn(5,1);
% % Y=10*randn(5,1);
% % Z=10*randn(5,1);
% n=2000;
% % X_Eval=X(n);
% % Y_Eval=Y(n);
% % Z_Eval=Z(n);
% X_Eval=10*randn(1,1e2);
% Y_Eval=10*randn(1,1e2);
% Z_Eval=10*randn(1,1e2);
% X=X_Eval;
% Y=Y_Eval;
% Z=Z_Eval;
% tic;
% [xfit_0,yfit_0]=cubic_calibration_eval(calibration_function,X_Eval,Y_Eval,Z_Eval);
% disp(toc);
% tic;
[xfit_1,yfit_1]=cubic_calibration_eval_01(calibration_function,X_Eval,Y_Eval,Z_Eval);

plot(xfit_1,yfit_1,'.BLACK');
axis equal;
axis([1,1024,1,1024]);
grid on;

keyboard;



% fprintf('\n\n\n\n');
% disp('[ x(n),    xfit_0,   x_fit_1,  y(n),     yfit_0,   yfit_1    ]');
% disp([xfit_0',xfit_1',yfit_0',yfit_1']);

%keyboard;





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



function test_calibration_function(x,y,X,Y,Cx,Cy);

CX0Y0=Cx(1);
CX1Y0=Cx(2);
CX2Y0=Cx(3);
CX3Y0=Cx(4);
CX0Y1=Cx(5);
CX0Y2=Cx(6);
CX0Y3=Cx(7);
CX1Y1=Cx(8);
CX2Y1=Cx(9);
CX1Y2=Cx(10);

for jj=1:length(X);
   
    x_fit(jj) = CX0Y0 + ...
        CX1Y0 * X(jj) + CX2Y0 * X(jj).^2 + CX3Y0 * X(jj).^3 + ...
        CX0Y1 * Y(jj) + CX0Y2 * Y(jj).^2 + CX0Y3 * Y(jj).^3 + ...
        CX1Y1 * X(jj) .* Y(jj) + CX2Y1 * X(jj).^2 .* Y(jj) + CX1Y2 * X(jj) .* Y(jj).^2;

end;

CX0Y0=Cy(1);
CX1Y0=Cy(2);
CX2Y0=Cy(3);
CX3Y0=Cy(4);
CX0Y1=Cy(5);
CX0Y2=Cy(6);
CX0Y3=Cy(7);
CX1Y1=Cy(8);
CX2Y1=Cy(9);
CX1Y2=Cy(10);

for jj=1:length(X);

    y_fit(jj) = CX0Y0 + ...
        CX1Y0 * X(jj) + CX2Y0 * X(jj).^2 + CX3Y0 * X(jj).^3 + ...
        CX0Y1 * Y(jj) + CX0Y2 * Y(jj).^2 + CX0Y3 * Y(jj).^3 + ...
        CX1Y1 * X(jj) .* Y(jj) + CX2Y1 * X(jj).^2 .* Y(jj) + CX1Y2 * X(jj) .* Y(jj).^2;
    
end;

figure(1);
plot(x,y,'oBLUE');
hold on;
plot(x_fit,y_fit,'oRED');
hold off;
axis equal;
grid on;
drawnow;



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





function [x,y]=cubic_calibration_eval(calibration_function,X,Y,Z);
% This function evaluates the cubic calibration functions at the world
% point (X,Y,Z) to yield the image point (x,y).

% This is a vector of the Z planes to determine which planes' calibration
% functions to use
Z_plane=zeros(size(calibration_function.z_world_value));
for ii=1:length(calibration_function.z_world_value);
    Z_plane_temp=calibration_function.z_world_value{ii};
    Z_plane(ii)=Z_plane_temp(1);
end;

% This initializes the output vectors
x=zeros(size(X));
y=zeros(size(X));
% This iterates through the points
for ii=1:length(X);
    % These are the current values
    X_Current=X(ii);
    Y_Current=Y(ii);
    Z_Current=Z(ii);
    % This determines the location of the current Z point with respect to
    % the Z calibration planes
    min_plane_index=sum(Z_plane<Z_Current);
    max_plane_index=min_plane_index+1;
    % These are the weighting coefficients
    if min_plane_index==0;
        min_plane_index=1;
        max_plane_index=1;
        w1=1;
        w2=0;
    elseif min_plane_index==length(Z_plane);
        max_plane_index=min_plane_index;
        w1=1;
        w2=0;
    else;
        % This is the first weighting coefficient
        w1=(Z_plane(max_plane_index)-Z_Current)/(Z_plane(max_plane_index)-Z_plane(min_plane_index));
        % This is the second weighting coefficient
        w2=(Z_Current-Z_plane(min_plane_index))/(Z_plane(max_plane_index)-Z_plane(min_plane_index));
    end;
    
    % disp(['w1 = ',num2str(w1),' and w2 = ',num2str(w2)]);

    
    % These are the linearly weighted x polynomial coefficients
    Cx=w1*calibration_function.x_coefficients(min_plane_index,:)+w2*calibration_function.x_coefficients(max_plane_index,:);
    % These are the linearly weighted y polynomial coefficients
    Cy=w1*calibration_function.y_coefficients(min_plane_index,:)+w2*calibration_function.y_coefficients(max_plane_index,:);
    
    % These are the current interpolated mean world coordinates
    X_Mean=w1*calibration_function.x_world_mean(min_plane_index)+w2*calibration_function.x_world_mean(max_plane_index);
    Y_Mean=w1*calibration_function.y_world_mean(min_plane_index)+w2*calibration_function.y_world_mean(max_plane_index);
    % These are the current interpolated standard deviation of the world coordinates
    X_Std=w1*calibration_function.x_world_std(min_plane_index)+w2*calibration_function.x_world_std(max_plane_index);
    Y_Std=w1*calibration_function.y_world_std(min_plane_index)+w2*calibration_function.y_world_std(max_plane_index);
    
    % This calculates the scaled X value
    X_Scale=(X_Current-X_Mean)/X_Std;
    % This calculates the scaled Y value
    Y_Scale=(Y_Current-Y_Mean)/Y_Std;
     
    % This is the polynomial independent vector
    XY_Vect=[1,X_Scale,X_Scale.^2,X_Scale.^3,Y_Scale,Y_Scale.^2,Y_Scale.^3,X_Scale.*Y_Scale,X_Scale.^2.*Y_Scale,X_Scale.*Y_Scale.^2];    

    
    % This is the x image value
    x(ii)=XY_Vect*Cx';
    % This is the y image value
    y(ii)=XY_Vect*Cy';
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

    
    
