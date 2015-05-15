function calculate_linear_calibration_01;
% This function takes a series of world and image point correspondences in
% the form of
%
%  x = f(X,Y,Z)
%  z = g(X,Y,Z)
%
% and generates the linear camera calibration matrix in a decomposed form
% given by
%
%  P = K * [ R | t ]
%
% such that both the intrinsic and extrinsic camera parameters may be
% determined.  The is accomplished by first calculating the linear matrix P
% from the list of point correspondences and then decomposing P into its
% constituent intrinsic and extrinsic components.

% This is the volume domain over which to perform the reconstruction
wrld_xmin=-0.28;
wrld_xmax=0.34;
wrld_ymin=-0.34;
wrld_ymax=0.28;
wrld_zmin=0.14;
wrld_zmax=1.78;
% This creates the volume domain vector
v_domain(1)=wrld_xmin;
v_domain(2)=wrld_xmax;
v_domain(3)=wrld_ymin;
v_domain(4)=wrld_ymax;
v_domain(5)=wrld_zmin;
v_domain(6)=wrld_zmax;
% This is the camera resolution
xres_cam=1024;
yres_cam=1024;
% This is the camera domain vector (this is somewhat meaningless, but is
% included for backwards compatibility)
cam_domain=[-1,1,-1,1];
% This states that the images should be saved
camera_images_save=true;
% This states whether to graph the camera calibration
display_calibration_gui=false;
% This is the calibration directory
calibration_directory='/mnt/current_storage/Projects/Argonne/Argonne_Tests_072910/Calibration_02/Calibration_Data_01/Simulated_Calibration_Images/';

% These are the calibration directories
calibration_files{1}='/mnt/current_storage/Projects/Argonne/Argonne_Tests_072910/Calibration_02/Calibration_Data_01/Calibration_Data/cam_1/cam_1_calibration.mat';
calibration_files{2}='/mnt/current_storage/Projects/Argonne/Argonne_Tests_072910/Calibration_02/Calibration_Data_01/Calibration_Data/cam_2/cam_2_calibration.mat';

% This creates the camera data structure
camera_data=struct([]);

% This iterates through the the calibration files (or the cameras)
for ii=1:length(calibration_files);
    % This loads the calibration data
    load(calibration_files{ii});
    % This creates a matrix w for inputing into the linear calibration
    % function
    w=ones(size(x));
    % This calculates the linear calibration matrix P
    P=linear_calibration_04(x,y,w,X,Y,Z);
    % This calculates the intrisic camera matrix, the camera rotation matrix,
    % and the translation vector
    [K,R,t]=linear_decomposition(P);
    % This enters the data into the camera data structure
    camera_data(ii).intrinsic=K;
    camera_data(ii).rotation=R;
    camera_data(ii).translation=t;
    camera_data(ii).perspective=P;
    camera_data(ii).cam_domain=cam_domain;
    camera_data(ii).resolution=[xres_cam,yres_cam];
end;

% This initializes the structure to save the data in
image_parameters=struct;
% This saves the volumne domain to the parameter structure
image_parameters.vol_domain=v_domain;
% This saves the boolean value determining whether or not to save the 2D
% simulated camera images
image_parameters.images_save=camera_images_save;
% This saves the image resolution to the structure
image_parameters.image_res=[xres_cam,yres_cam];
% testing_parameters.camera_angles=[x_camera_theta,y_camera_theta];
% This saves the camera calibration gui display boolean value
image_parameters.display_gui=display_calibration_gui;
% This saves the directory containging the calibration data to the
% structure
image_parameters.calibration_directory=calibration_directory;
% This saves the camera data to the image_parameters structure
image_parameters.camera_data=camera_data;

% This generates the calibration images
calibration_grid_simulation_01(image_parameters);




function graph_calibration(v_domain,camera_data);
% This function graphs the interogation domain as well as the camera
% positions, directions, et cetera.  Basically this function gives a
% pictorial representation of the camera calibration with respect to the
% interogation volume.

% This extracts the domain limits from the input vector
wrld_xmin=v_domain(1);
wrld_xmax=v_domain(2);
wrld_ymin=v_domain(3);
wrld_ymax=v_domain(4);
wrld_zmin=v_domain(5);
wrld_zmax=v_domain(6);
% This is a vector of the origin coordinates of the interogation volume
origin_vect=[(wrld_xmin+wrld_xmax)/2;(wrld_ymin+wrld_ymax)/2;(wrld_zmin+wrld_zmax)/2];
% This is the vertex matrix for drawing the interogation volume using the
% patch function
[x_vert,y_vert,z_vert]=meshgrid([wrld_xmin,wrld_xmax],[wrld_ymin,wrld_ymax],[wrld_zmin,wrld_zmax]);
vertex_matrix=[x_vert(:),y_vert(:),z_vert(:)];
% This is the face matrix for drawing the interogation volume using the
% patch funciton
face_matrix=[   1,3,7,5;
                3,4,8,7;
                4,2,6,8;
                2,1,5,6;
                1,3,4,2;
                5,7,8,6     ];
% This is the current figure
figure(1);
% This clears the current figure
clf;
% This creates the patch object of the cube
h=patch('Vertices',vertex_matrix,'Faces',face_matrix,'FaceColor','BLACK','EdgeColor','BLACK');
% This sets the transparency of the cube
alpha(h,0.1);
% These are the colors to iterate through for the different cameras
color_vect='cmyrgb';
% This is a vector used to project the cameras back into space
lambda=0:0.1:4;
% This iterates through the cameras plot thier positions and orientations
% in space
for ii=1:length(camera_data);
    % This is the camera center
    C=camera_data(ii).translation;
    % This is the camera rotation matrix
    R=camera_data(ii).rotation;
    % This is the intrinsic matrix of the camera
    K=camera_data(ii).intrinsic;
    % These are the minimum and maximum coordinates in each camera plane
    cam_domain=camera_data(ii).cam_domain;
    cam_xmin=cam_domain(1);
    cam_xmax=cam_domain(2);
    cam_ymin=cam_domain(3);
    cam_ymax=cam_domain(4);
    % This is the resolution of the camera
    cam_res=camera_data(ii).resolution;
    xres=cam_res(1);
    yres=cam_res(2);
    % These are the scaling factors for calculating the intrinsic matrix
    alpha_x=(xres-1)/(cam_xmax-cam_xmin);
    alpha_y=(yres-1)/(cam_ymax-cam_ymin);
    
    % These are the instrinsic parameters of K in pixel coordinates for
    % converting to world coordinates
    fk_u=K(1,1)/alpha_x;
    fk_v=K(2,2)/alpha_y;
    gamma=-K(1,2)/alpha_x;
    u_0=(alpha_x*cam_xmax+1-K(1,3))/alpha_x;
    v_0=(alpha_y*cam_ymax+1-K(2,3))/alpha_y;
    % This converts the intrinsic matrix from pixel coordinates to world
    % coordinates
    K=[-fk_u,gamma,u_0;0,-fk_v,v_0;0,0,1];
       
    % This generates vectors corresponding to the limits (the corners) of
    % the camera coordinate plane
    [x_vect,y_vect]=meshgrid([cam_xmin,cam_xmax],[cam_ymin,cam_ymax]);
    x_vect=x_vect(:);
    y_vect=y_vect(:);
    % This permutes the vectors so that when plotted they form a square
    % rather than two triangles point to point
    x_vect=x_vect([1;2;4;3],:);
    y_vect=y_vect([1;2;4;3],:);
    % This iterates through the corners of each camera
    for jj=1:4;
        % This is the color to plot the current camera
        plot_color=color_vect(mod(ii,length(color_vect))+1);
        % This generates the first ray of the current camera
        %
        % This is the point corresponding to the corner of the camera
        % coordinate system
        p=[x_vect(jj);y_vect(jj);1];
        % This calculates the projection of the camera corner
        ray_projection=repmat(C,1,size(lambda,2))+repmat(lambda,3,1).*repmat((K*R)\p,1,size(lambda,2));
        x_wrld_ray_1=ray_projection(1,:);
        y_wrld_ray_1=ray_projection(2,:);
        z_wrld_ray_1=ray_projection(3,:);
        % This generates the second ray of the current camera
        %
        % If the index equals four, then it is matched with the first index
        % (to create a closed square)
        if jj+1>4;
            jj=0;
        end;
        % This is the point corresponding to the corner of the camera
        % coordinate system
        p=[x_vect(jj+1);y_vect(jj+1);1];
        % This calculates the projection of the camera corner
        ray_projection=repmat(C,1,size(lambda,2))+repmat(lambda,3,1).*repmat((K*R)\p,1,size(lambda,2));
        x_wrld_ray_2=ray_projection(1,:);
        y_wrld_ray_2=ray_projection(2,:);
        z_wrld_ray_2=ray_projection(3,:);
        % This specifies the vertex matrix of the edge of the camera
        % pyramid
        vertex_matrix=[ x_wrld_ray_1(1),    y_wrld_ray_1(1),    z_wrld_ray_1(1);
                        x_wrld_ray_1(end),  y_wrld_ray_1(end),  z_wrld_ray_1(end);
                        x_wrld_ray_2(end),  y_wrld_ray_2(end),  z_wrld_ray_2(end)   ];
        % This specifies the face matrix of the camera pyramid
        face_matrix=[1,2,3];
        % This plots the current surface of the camera pyramid
        h=patch('Vertices',vertex_matrix,'Faces',face_matrix,'FaceColor',plot_color,'EdgeColor','BLACK');
        % This sets the transparency of the pyramid
        alpha(h,0.5);
        % This plots the camera center point
        hold on;
        plot3(C(1),C(2),C(3),['o',plot_color]);
        hold off;
    end;
    % This plots the principal ray
    %
    % This sets the independent variable generating the ray to a length
    % equal to the distance between the camera center and the origin of the
    % interogation volume
    lambda_principal=[0,norm(C-origin_vect)];
    % This is the center of the camera image plane coordinate
    p=[0;0;1];
    % This calculates the principal ray projection
    ray_projection=repmat(C,1,size(lambda_principal,2))+repmat(lambda_principal,3,1).*repmat(inv(R)*inv(K)*p,1,size(lambda_principal,2));
    x_wrld_ray=ray_projection(1,:);
    y_wrld_ray=ray_projection(2,:);
    z_wrld_ray=ray_projection(3,:);   
    % This plots the principal ray projection
    hold on;
    plot3(x_wrld_ray,y_wrld_ray,z_wrld_ray,'--oBLACK');
    hold off;
end;
% This sets the axis to equal proportions, turns the grid on, and sets the
% default view to roughly isometric
axis equal;
grid on;
view(3);
% This changes the axis coordinates such that a cube centered around the
% cameras and interogation volume is plotted
%
% This is the current axis coordinates
axis_vect=axis;
% This is the range of the axis coordinates
axis_dx=diff(axis_vect(1:2));
axis_dy=diff(axis_vect(3:4));
axis_dz=diff(axis_vect(5:6));
% This sets the remaining two axis to the same range as the largest range
% axis (this may lead to axes that are too large due to the dumb way that
% matlab's 'axis equal' command works)
if all([axis_dx==axis_dy,axis_dx==axis_dz,axis_dy==axis_dz]);
    axis(axis_vect);
elseif axis_dx==max([axis_dx,axis_dy,axis_dz]);
    axis_vect(3)=axis_vect(3)-(axis_dx-axis_dy)/2;
    axis_vect(4)=axis_vect(4)+(axis_dx-axis_dy)/2;
    axis_vect(5)=axis_vect(5)-(axis_dx-axis_dz)/2;
    axis_vect(6)=axis_vect(6)+(axis_dx-axis_dz)/2;
    axis(axis_vect);
elseif axis_dy==max([axis_dx,axis_dy,axis_dz]);
    axis_vect(1)=axis_vect(1)-(axis_dy-axis_dx)/2;
    axis_vect(2)=axis_vect(2)+(axis_dy-axis_dx)/2;
    axis_vect(5)=axis_vect(5)-(axis_dy-axis_dz)/2;
    axis_vect(6)=axis_vect(6)+(axis_dy-axis_dz)/2;
    axis(axis_vect);
elseif axis_dz==max([axis_dx,axis_dy,axis_dz]);
    axis_vect(1)=axis_vect(1)-(axis_dz-axis_dx)/2;
    axis_vect(2)=axis_vect(2)+(axis_dz-axis_dx)/2;
    axis_vect(3)=axis_vect(3)-(axis_dz-axis_dy)/2;
    axis_vect(4)=axis_vect(4)+(axis_dz-axis_dy)/2;
    axis(axis_vect);
end;
% This labels the coordinate axis
xlabel('X');
ylabel('Y');
zlabel('Z');



function [K,R,t]=linear_decomposition(P);
% This function decomposes the linear camera calibration P into the
% intrinsic and extrinsic matricies.  This is taken from section 6.2.4 of 
% "Multiple View Geometry in Computer Vision" by Richard Hartley and Andrew
% Zisserman.

% These are the linear matrix column vectors
p1=P(:,1);
p2=P(:,2);
p3=P(:,3);
p4=P(:,4);
% These are the coordinates of the camera center (calculated from the fact
% that P*C = 0
C_X=det([p2,p3,p4]);
C_Y=-det([p1,p3,p4]);
C_Z=det([p1,p2,p4]);
C_T=-det([p1,p2,p3]);
% This is the inhomogeneous camera center
C=[C_X;C_Y;C_Z;C_T];
% This the camera center
C_Tilde=C(1:3)/C(4);
% This is the left 3 x 3 submatrix of P
M=P(:,1:3);
% This performs an RQ matrix decomposition of the left 3 x 3 matrix
[K,R]=rq_decomposition(M);
% This is the translation matrix
t=-R*C_Tilde;



function [R,Q]=rq_decomposition(A);
% This function performs the RQ decomposition of the matrix A such that
%
%  A = R * Q
%
% where R is a right upper-triangular matrix and Q is an orthogonal matrix.
% The matrix A must be a n x n matrix.

% This is the size of A
n=size(A,1);
% This creates a permutation matrix
p=fliplr(eye(n));
% This performs the QR decomposition
[Q,R]=qr((p*A)');
% This permutes the R matrix
R=p*R'*p;
% This permutes the Q matrix
Q=p*Q';



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
