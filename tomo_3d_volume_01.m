function tomo_3d_volume_01;
% This function uses the MART algorithm along with Elsinga's formulation to
% reconstruct a 3d intensity field from piv images.

% This is the domain over which to reconstruct the volume
xmin=-1;
xmax=1;
ymin=-1;
ymax=1;
zmin=11;
zmax=13;
% zmin=12-0.025;
% zmax=12+0.025;
% This is the width of the voxels; in the actual implementation this will
% be set approximately equal to the pixel magnification in the images
dx=0.05;
dy=0.05;
dz=0.05;
% This then resets [dx,dy,dz] to values that give integral divisions of the
% volume so that the minimum/maximum voxel coordinates correspond to the
% volume domain minimum/maximum
dx_voxel=(xmax-xmin)/round((xmax-xmin)/dx);
dy_voxel=(ymax-ymin)/round((ymax-ymin)/dy);
dz_voxel=(zmax-zmin)/round((zmax-zmin)/dz);
% This is the radius of the sphere centered on the voxels
r=(3*dx_voxel*dy_voxel*dz_voxel/(4*pi))^(1/3);
% This is the volume of the sphere centered on the voxels
V_Sphere=dx_voxel*dy_voxel*dz_voxel;
% This creates a coordinate vectors for the volume
X=xmin:dx_voxel:xmax;
Y=ymin:dy_voxel:ymax;
Z=zmin:dz_voxel:zmax;
% This size of the volume matrix
L=length(X);
M=length(Y);
N=length(Z);
% This initializes the volume intensity matrix (in linear form)
% E=zeros(L*M*N,1);

% imagesc(w(2000:3000,30000:31000));
% colormap('gray');
% axis equal;

[X_Full,Y_Full,Z_Full]=meshgrid(X,Y,Z);

load('/home/rlafoy/Documents/AEThER/TomoPIV/Image_Simulation_01/camera_data_01.mat');

% P=[-14,0,0,0;0,-14,0,-7;0,0,1,0];
% x0_vect=-0.5:0.01:0.5;
% y0_vect=-0.5:0.01:0.5;
% r=0.1;
% tic;
% for ii=1:length(x0_vect);
%     for jj=1:length(y0_vect);
%         I=pixel_distance_test(x0_vect(ii),y0_vect(jj),x_cam,y_cam,r);
%     end;
% end;
% t=toc/(length(x0_vect)*length(y0_vect));
% disp(t);
% I=pixel_distance_test(x0,y0,x_cam,y_cam,r);
% x_plot=X_Full(I(:));
% y_plot=Y_Full(I(:));
% z_plot=Z_Full(I(:));
% plot3(x_plot,y_plot,z_plot,'oBLUE');
% grid on;
% axis([xmin,xmax,ymin,ymax,zmin,zmax]);



% return;


% camera_data(1).resolution=[1024,1024];
% camera_data(2).resolution=[1024,1024];
% 
% camera_data(1).perspective=[-14,0,0,0;0,-14,0,-3.5;0,0,1,0];
% camera_data(2).perspective=[-14,0,0,0;0,-14,0,-7;0,0,1,0];
% 
% camera_data(1).cam_domain=[-2,2,-2,2];
% camera_data(2).cam_domain=[-2,2,-2,2];



% This initializes a variable that is the total number of pixels
pixel_num=0;
% This iterates through the cameras to calculate the size of w
for n=1:length(camera_data);
    % This extracts the resolution of the nth camera
    cam_res=camera_data(n).resolution;
    % This adds this resolution to the total pixel number
    pixel_num=pixel_num+prod(cam_res);
end;
% This creates the weighting matrix
w=sparse(pixel_num,L*M*N);
% This intializes the pixel indexing
pixel_index=0;
% This iterates through the cameras calculating the weight matrix
% components

tic;

for n=1:length(camera_data);
    % This loads the perspective projection matrix
    P=camera_data(n).perspective;
    % This is the pseudo-inverse of P
    P_Inv=(P')*inv(P*(P'));
    % This generates the camera center point
    C=null(P);
    % This extracts the resolution of the nth camera
    cam_res=camera_data(n).resolution;
    x_cam_res=cam_res(1);
    y_cam_res=cam_res(2);
    % This is the spatial extent of the image in global coordinates
    cam_extent=camera_data(n).cam_domain;
    xmin_cam=cam_extent(1);
    xmax_cam=cam_extent(2);
    ymin_cam=cam_extent(3);
    ymax_cam=cam_extent(4);
    % These are the widths of the individual pixels
    dx_pixel=(xmax_cam-xmin_cam)/x_cam_res;
    dy_pixel=(ymax_cam-ymin_cam)/y_cam_res;
    % This is the radius of the cylinder in terms of the pixel sizes
    % (specifically the area of the cylinder and the pixel are set equal)
    R=sqrt(dx_pixel*dy_pixel/pi);
    
    % This calculates the voxel coordinates in terms of the current cameras
    % coordinate system
    [x_cam,y_cam]=world_to_cam_convert(X_Full(:),Y_Full(:),Z_Full(:),P);
    
    % These vectors specify the position of the pixels in camera
    % coordiantes
%     x_pos=xmax_cam-(xmax_cam-xmin_cam)*(1-(0:x_cam_res-1))/(x_cam_res-1);
%     y_pos=ymax_cam-(ymax_cam-ymin_cam)*(1-(0:y_cam_res-1))/(y_cam_res-1);
    
    x_pos=((xmax_cam-xmin_cam)*(1:x_cam_res)+xmin_cam*x_cam_res-xmax_cam)/(x_cam_res-1);
    y_pos=((ymax_cam-ymin_cam)*(1:y_cam_res)+ymin_cam*y_cam_res-ymax_cam)/(y_cam_res-1);
    
    % This iterates through the pixels of the camera to calculate the
    % components of the weighting matrix
    for ii=1:y_cam_res;
        
        disp(' ');
        disp([num2str(round(100*ii/y_cam_res)),'% of the calculation of camera number ',num2str(n),' interpolation matrix complete . . .']);
        disp(['Elapsed time equals ',num2str(toc),' seconds . . . ']);
        
        for jj=1:x_cam_res;
            % This increments the pixel indexing
            pixel_index=pixel_index+1;
            % This is the current pixel coordinate
            m=[x_pos(jj);y_pos(ii);1]; %%%%% I am not sure if this is the correct form of m
            % This generates two points along the line of sight of the
            % pixel for generation of the infinitely extended cylinder
%             P1=P_Inv*m;
%             P2=P1+C;

            P1=C;
            P2=P_Inv*m+C;
            
            P1=P1(1:3)/P1(4);
            P2=P2(1:3)/P2(4);
            
            % This is the first point specifying the cylinder
            xc1=P1(1);
            yc1=P1(2);
            zc1=P1(3);
            % This is the second point specifying the cylinder
            xc2=P2(1);
            yc2=P2(2);
            zc2=P2(3);
            % This returns a logical index to all reconstruction voxels
            % that are within 1 radius of the line of sight of the current
            % pixel
            
%             fprintf('\n\n\n\n');
%             disp([min(x_cam(:)),x_pos(jj),max(x_cam(:))]);
%             disp([min(y_cam(:)),y_pos(ii),max(y_cam(:))]);
%             pause(0.01);
            
            voxel_index=pixel_distance_test(x_pos(jj),y_pos(ii),x_cam,y_cam,3*R);
            
            
            
            % These are the voxel coordinates
            xs=X_Full(voxel_index(:));
            ys=Y_Full(voxel_index(:));
            zs=Z_Full(voxel_index(:));
            
%             if sum(voxel_index)>0;
%                 t_plot_min=min([(zmin-zc2)/(zc1-zc2),(zmax-zc2)/(zc1-zc2)]);
%                 t_plot_max=max([(zmin-zc2)/(zc1-zc2),(zmax-zc2)/(zc1-zc2)]);
%                 t_plot_vect=[t_plot_min,t_plot_max];
%                 x_plot_vect=(xc1-xc2)*t_plot_vect+xc2;
%                 y_plot_vect=(yc1-yc2)*t_plot_vect+yc2;
%                 z_plot_vect=(zc1-zc2)*t_plot_vect+zc2;
%                 plot3(xs,ys,zs,'oBLUE');
%                 hold on;
%                 plot3(x_plot_vect,y_plot_vect,z_plot_vect,'-BLACK');
%                 hold off;
%                 grid on;
%                 axis([xmin,xmax,ymin,ymax,zmin,zmax]);
%                 pause(2);
%             end;
            
            % This computes the impact diameters b (the distance between
            % the spheres centered on the voxel positions and the center of
            % the cylinder)
            dxc=xc2-xc1;
            dyc=yc2-yc1;
            dzc=zc2-zc1;
            t_min=(dxc*(xs-xc1)+dyc*(ys-yc1)+dzc*(zs-zc1))/(dxc^2+dyc^2+dzc^2);
            b=sqrt((xs-(xc1+dxc*t_min)).^2+(ys-(yc1+dyc*t_min)).^2+(zs-(zc1+dzc*t_min)).^2);
            % This finds the volumes of interesction of the cylinders and
            % spheres
            V_Intersect=cylinder_sphere_intersect_05(r*ones(size(b)),R*ones(size(b)),b);
            % The values of the weighting functions for these voxels are
            w_temp=V_Intersect/V_Sphere;
            % This saves the interpolation values given by the volume
            % intersections to huge interpolation matrix
            w(pixel_index,voxel_index)=w_temp;
            
            % I now have to save these values to their approriate place in
            % w and then I am done and everything will work wonderfully and
            % it will only take seconds and there will be rainbows and
            % sunshine and unicorns and . . .
            
        end;
    end;
         
end;

% imagesc(w(2000:3000,30000:31000));
% colormap('gray');
% axis equal;

whos('w');
disp(max(w(:)));
disp(sum(w(:)));




function I=pixel_distance_test(x0,y0,x_cam,y_cam,r);
% This function returns the indicies of all points in the variables x_cam
% and y_cam that are less then or equal to a distnace r away from x0 and
% y0.  The variables x0, y0, and r must be scalars.  The variables x_cam
% and y_cam may be scalars, vectors, or matricies.  The returned indicies
% will be of the same dimension of x_cam and y_cam.
%
% These are the logical indicies
I=(x_cam-x0).^2+(y_cam-y0).^2<=r^2;



function [x_cam,y_cam]=world_to_cam_convert(x_world,y_world,z_world,P);
% This function maps a set of points from world coordinates to camera
% coordinates given the perspective projection matrix P.  The inputs
% x_world, y_world, z_world may be scalars, vectors, or matricies.  The
% outputs x_cam and y_cam will be of the same size as x_world, y_world, and
% z_world.
%
% This is the size of the world coordinate vectors
[L,M,N]=size(x_world);
% If the input world coordinates are in the form of matricies
if (L>1)&&(M>1);
    % This is the size of the full matrix (or the length of it's vector
    % form)
    LMN=L*M*N;
    % This converts the matricies to vectors
    x_world_vect=x_world(:)';
    y_world_vect=y_world(:)';
    z_world_vect=z_world(:)';
    % This is the world coordinate matrix
    w=[x_world_vect;y_world_vect;z_world_vect;ones(size(x_world_vect))];
    % This is the temporary camera coordinate matric
    m_temp=P*w;
    % Then since this matrix may be quite large, we first try a fast way
    % and if that doesn't work, try the slow and stupid way which will work
    try;
        % This attempts to divide m by the z coordinate using faster matrix
        % multiplication
        %
        % This is a diagonal matrix for dividing P*w by z_world
        D=sparse(1:LMN,1:LMN,1./z_world_vect,LMN,LMN,LMN);
        
%         D=sparse(1:LMN,1:LMN,1./m_temp(3,:),LMN,LMN,LMN);
        
        % This is the camer coordinate matrix
        m=m_temp*D;
    catch;
        % If the matrix multiplication fails, a slow for loop is used which
        % required less memory
        m=zeros(3,LMN);
        for n=1:LMN;
            m(:,n)=m_temp(:,n)/z_world_vect(n);
            
%             m(:,n)=m_temp(:,n)/m_temp(3,n);
            
        end;
    end;
    % These are the camera coordinates
    x_cam_vect=m(1,:);
    y_cam_vect=m(2,:);
    % This reshapes the output to the same size as the input
    x_cam=reshape(x_cam_vect,L,M,N);
    y_cam=reshape(y_cam_vect,L,M,N);
elseif ((L==1)&&(N==1))||((M==1)&&(N==1));
    % This is the length of the vectors
    LMN=length(x_world);
    % If the coordinates are column vectors, they are converted to row
    % vectors to aid in the maths
    if M==1;
        x_world_vect=x_world';
        y_world_vect=y_world';
        z_world_vect=z_world';
    end;
    % This is the world coordinate matrix
    w=[x_world_vect;y_world_vect;z_world_vect;ones(size(x_world_vect))];
    % This is the temporary camera coordinate matric
    m_temp=P*w;
    % Then since this matrix may be quite large, we first try a fast way
    % and if that doesn't work, try the slow and stupid way which will work
    try;
        % This attempts to divide m by the z coordinate using faster matrix
        % multiplication
        %
        % This is a diagonal matrix for dividing P*w by z_world
        D=sparse(1:LMN,1:LMN,1./z_world_vect,LMN,LMN,LMN);
        
%         D=sparse(1:LMN,1:LMN,1./m_temp(3,:),LMN,LMN,LMN);
        
        % This is the camer coordinate matrix
        m=m_temp*D;
    catch;
        % If the matrix multiplication fails, a slow for loop is used which
        % required less memory
        m=zeros(3,LMN);
        for n=1:LMN;
            m(:,n)=m_temp(:,n)/z_world_vect(n);
            
%             m(:,n)=m_temp(:,n)/m_temp(3,n);
            
        end;
    end;
    % These are the camera coordinates
    x_cam=m(1,:);
    y_cam=m(2,:);
    % If the coordinates were column vectors, they are converted back to
    % column vectors
    if M==1;
        x_cam=x_cam';
        y_cam=y_cam';
    end;
end;



function V=cylinder_sphere_intersect_05(r,R,b);
% This function calculates the volume of intersection of an arbitrarily
% defined sphere in three-space and an infinitely extended cylinder.  This
% version is the same as version 04 except that it attempts to be a bit
% more optimized.  The vectors r,R,b must be column vectors.
%
%   r - Sphere Radius
%   R - Cylinder Radius
%   b - Impact Diameter
%
% This is the elliptic integral tolerance
tol=1e-4;
% These are the constants used for the volume calculation
A=max([r.^2,(b+R).^2],[],2);
B=min([r.^2,(b+R).^2],[],2);
C=(b-R).^2;
k_sqr=(B-C)./(A-C);
k=sqrt(k_sqr);
alpha_sqr=-(B-C)./C;
s=(b+R).*(b-R);
% This computes the Heaviside step function of 'R-b'
H=heaviside_step(R-b);
% This computes the first complete elliptic integral
K=first_elliptic(k,tol);
% This computes the second complete elliptic integral
E=second_elliptic(k,tol);
% This computes the third complete elliptic integral
PI=third_elliptic(k,alpha_sqr,tol);

% This calculates some useful values for the volume calculation (in order
% to reduce redundant calculations)
V_Sphere=(4*pi*r.^3)/3;
AC_Sqrt=sqrt(A-C);
A_Sqrt=sqrt(A);

% This initializes the V matrix
V_Temp=zeros(size(r,1),8);
% If b==0 and r>=R
V_Temp(:,1)=V_Sphere-(4*pi/3)*(r.^2-R.^2).^(3/2);
% If (b+r)<R
V_Temp(:,2)=V_Sphere;
% If b>(r+R)
V_Temp(:,3)=zeros(size(r));
% If r<(b+R) and b~=R
V_Temp(:,4)=V_Sphere.*H+((4/3)./AC_Sqrt).*((PI.*s.*B.^2)./C+K.*(s.*(A-2*B)+(A-B).*(3*B-C-2*A)/3)+E.*(A-C).*((2*A+2*C-4*B)/3-s));
% If r<(b+R) and b==R
V_Temp(:,5)=V_Sphere/2+((4/3)./A_Sqrt).*(K.*(A-B).*(3*B-2*A)/3+E.*A.*(2*A-4*B)/3);
% If r==b+R
V_Temp(:,6)=V_Sphere.*H+(V_Sphere/pi).*atan((2*sqrt(b.*R))./(b-R))-(4/3)*AC_Sqrt.*(s+(2/3)*(A-C));
% If r>(b+R) and b~=R
V_Temp(:,7)=V_Sphere.*H+((4/3)./AC_Sqrt).*((PI.*s.*A.^2)./C-K.*(A.*s-(A-B).*(A-C)/3)-E.*(A-C).*((4*A-2*B-2*C)/3+s));
% If r>(b+R) and b==R
V_Temp(:,8)=V_Sphere/2+((4/3)./A_Sqrt).*(K.*A.*(A-B)/3-E.*A.*(4*A-2*B)/3);

% This initializes the conditional matrix
Cond=zeros(size(r,1),8);
% This is the conditional matrix
Cond(:,1)=((b==0).*(r>=R));
Cond(:,2)=((b+r)<R);
Cond(:,3)=(b>(r+R));
Cond(:,4)=((r<(b+R)).*(b~=R).*(b~=0).*(R<=(b+r)).*(b<=(r+R)));
Cond(:,5)=((r<(b+R)).*(b==R));
Cond(:,6)=((r==(b+R)).*(b~=0));
Cond(:,7)=((r>(b+R)).*(b~=R).*(b~=0));
Cond(:,8)=((r>(b+R)).*(b==R).*(b~=0));

% This removes indeterminate values from V_Temp
V_Temp(logical(isnan(V_Temp)))=0;
% This is the output volume
V=sum(V_Temp.*Cond,2);



function H=heaviside_step(x);
% This function computes the Heaviside step function of x
H=0*(x<=0)+1*(x>0);



function K=first_elliptic(k,tol);
% This is the first complete elliptic integral.  The vector k must be a
% column vector.
%
% This checks whether any values of k are infinite or indeterminate
k_indeterminate=logical(isinf(k)+isnan(k)+(k==1));
k(k_indeterminate)=0;
% This evaluates the elliptic integral
K=symmetric_first(zeros(size(k)),1-k.^2,ones(size(k)),tol);
% This replaces the indeterminate values of K
K(k_indeterminate)=NaN;



function E=second_elliptic(k,tol);
% This is the second complete elliptic integral.  The vector k must be a
% column vector.
%
% This checks whether any values of k are infinite or indeterminate
k_indeterminate=logical(isinf(k)+isnan(k)+(k==1));
k(k_indeterminate)=0;
% This evaluates the elliptic integral
E=(symmetric_third_degen(zeros(size(k)),1-k.^2,ones(size(k)),tol)+...
    symmetric_third_degen(zeros(size(k)),ones(size(k)),1-k.^2,tol)).*(1-k.^2)/3;
% This replaces the indeterminate values of E
E(k_indeterminate)=NaN;


function Pi=third_elliptic(k,n,tol);
% This is the third imcomplete elliptic integral.  The vectors k and n must
% be column vectors.
%
% This checks whether any values of k or n are infinite or indeterminate
kn_indeterminate=logical(isinf(k)+isnan(k)+(k==1)+isinf(n)+isnan(n));
k(kn_indeterminate)=0;
n(kn_indeterminate)=0;
% This evaluates the elliptic integrals
RF=symmetric_first(zeros(size(k)),1-k.^2,ones(size(k)),tol);
RJ=symmetric_third(zeros(size(k)),1-k.^2,ones(size(k)),1-n,tol);
Pi=RF+RJ.*n/3;
% This replaces the indeterminate values of Pi
Pi(kn_indeterminate)=NaN;




function RF=symmetric_first(x,y,z,tol);
% This function calculates a symmetric integral of the first kind to a
% relative error less then tol.  The vectors x,y,z must be column vectors.

% This is the maximum allowed number of iterations
nmax=20;
% This saves the variables initial values
x0=x;
y0=y;
z0=z;
% This is the first A variable
A=(x0+y0+z0)/3;
% This saves the initial A value
A0=A;
% This is the convergence test variable
Q=max(abs([A0-x0,A0-y0,A0-z0]),[],2)*(3*tol)^(-1/6);
% This initializes a counting variable
n=0;
% This iterates through the algorithm
while true;
    % This saves some variables for faster computation
    x_sqrt=sqrt(x);
    y_sqrt=sqrt(y);
    z_sqrt=sqrt(z);
    % This is the algorithm
    n=n+1;
    lambda=x_sqrt.*y_sqrt+x_sqrt.*z_sqrt+y_sqrt.*z_sqrt;
    A=(A+lambda)/4;
    if all(Q*4^(-n)<abs(A));
        break;
    end;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    z=(z+lambda)/4;
    % This checks to see if too many iterations are being carried out
    if n>nmax;
        error('The solution has not converged in the maximum number of iterations . . . ');
    end;
end;
% These are the coefficient arguments
X=(A0-x0)./(A*4^n);
Y=(A0-y0)./(A*4^n);
Z=-X-Y;
% These are the polynomial coefficients
E2=X.*Y-Z.^2;
E3=X.*Y.*Z;
% This is the elliptic integral approximation
RF=(1-E2/10+E3/14+(E2.^2)/24-3*E2.*E3/44)./sqrt(A);



function RC=symmetric_first_degen(x,y,tol);
% This function evaluates the degenerate case of symmetric integral of the
% first kind where RC(x,y) = RF(x,y,y).  The variables x,y can be vectors.

% This is the maximum allowed number of iterations
nmax=20;
% This saves the variables initial values
x0=x;
y0=y;
% This is the first A variable
A=(x0+2*y0)/3;
% This saves the initial A value
A0=A;
% This is the convergence test variable
Q=abs(A0-x0)*(3*tol)^(-1/8);
% This initializes a counting variable
n=0;
% This iterates through the algorithm
while true;
    n=n+1;
    lambda=2*sqrt(x).*sqrt(y)+y;
    A=(A+lambda)/4;
    if all(Q*4^(-n)<abs(A));
        break;
    end;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    % This checks to see if too many iterations are being carried out
    if n>nmax;
        error('The solution has not converged in the maximum number of iterations . . . ');
    end;
end;
% This is the independent variable of the polynomial
s=(y0-A0)./(A*4^n);
% These are the polynomial coefficients
P=[9/8,159/208,9/22,3/8,1/7,3/10,0,1];
% This is the elliptic integral approximation
RC=polyval(P,s)./sqrt(A);


function RJ=symmetric_third(x,y,z,p,tol);
% This function calculates a symmetric integral of the third kind to a
% relative error less then tol.  The vectors x,y,z,p must be column
% vectors.

% This is the maximum allowed number of iterations
nmax=20;
% This saves the variables initial values
x0=x;
y0=y;
z0=z;
p0=p;
% This is the first A variable
A=(x0+y0+z0+2*p0)/5;
% This saves the initial A value
A0=A;
% This is the delta variable
delta=(p0-x0).*(p0-y0).*(p0-z0);
% This is the convergence test variable
Q=max(abs([A0-x0,A0-y0,A0-z0,A0-p0]),[],2)*(tol/4)^(-1/6);
% This initializes a counting variable
n=0;
% This initializes the sum variable
sum=zeros(size(x));
% This iterates through the algorithm
while true;
    % This saves some variables for faster computation
    x_sqrt=sqrt(x);
    y_sqrt=sqrt(y);
    z_sqrt=sqrt(z);
    p_sqrt=sqrt(p);
    % This is the algorithm
    d=(p_sqrt+x_sqrt).*(p_sqrt+y_sqrt).*(p_sqrt+z_sqrt);
    e=(delta*4^(-3*n))./(d.^2);
    RC=symmetric_first_degen(1,1+e,tol);
    sum=sum+RC*(4^(-n))./d;
    n=n+1;
    lambda=x_sqrt.*y_sqrt+x_sqrt.*z_sqrt+y_sqrt.*z_sqrt;
    A=(A+lambda)/4;
    if all(Q*4^(-n)<abs(A));
        break;
    end;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    z=(z+lambda)/4;
    p=(p+lambda)/4;
    % This checks to see if too many iterations are being carried out
    if n>nmax;
        error('The solution has not converged in the maximum number of iterations . . . ');
    end;
end;
% These are the coefficient arguments
X=(A0-x0)./(A*4^n);
Y=(A0-y0)./(A*4^n);
Z=(A0-z0)./(A*4^n);
P=(-X-Y-Z)/2;
% These are the polynomial coefficients
E2=X.*Y+X.*Z+Y.*Z-3*P.^2;
E3=X.*Y.*Z+2*E2.*P+4*P.^3;
E4=(2*X.*Y.*Z+E2.*P+3*P.^3).*P;
E5=X.*Y.*Z.*P.^2;
% This is the elliptic integral approximation
RJ=(4^(-n))*(A.^(-3/2)).*(1-3*E2/14+E3/6+9*(E2.^2)/88-3*E4/22-9*E2.*E3/52+3*E5/26)+6*sum;



function RD=symmetric_third_degen(x,y,z,tol);
% This function calculates a symmetric integral of the third kind to a
% relative error less then tol.  The variables x,y,z must be column
% vectors.

% This is the maximum allowed number of iterations
nmax=20;
% This saves the variables initial values
x0=x;
y0=y;
z0=z;
% This is the first A variable
A=(x0+y0+3*z0)/5;
% This saves the initial A value
A0=A;
% This is the convergence test variable
Q=max(abs([A0-x0,A0-y0,A0-z0]),[],2)*(tol/4)^(-1/6);
% This initializes a counting variable
n=0;
% This initializes the sum variable
sum=zeros(size(x));
% This iterates through the algorithm
while true;
    % This saves some variables for faster computation
    x_sqrt=sqrt(x);
    y_sqrt=sqrt(y);
    z_sqrt=sqrt(z);
    % This is the algorithm
    lambda=x_sqrt.*y_sqrt+x_sqrt.*z_sqrt+y_sqrt.*z_sqrt;
    sum=sum+(4^(-n))./((z+lambda).*z_sqrt);
    n=n+1;
    A=(A+lambda)/4;
    if all(Q*4^(-n)<abs(A));
        break;
    end;
    x=(x+lambda)/4;
    y=(y+lambda)/4;
    z=(z+lambda)/4;
    % This checks to see if too many iterations are being carried out
    if n>nmax;
        error('The solution has not converged in the maximum number of iterations . . . ');
    end;
end;
% These are the coefficient arguments
X=(A0-x0)./(A*4^n);
Y=(A0-y0)./(A*4^n);
Z=-(X+Y)/3;
% These are the polynomial coefficients
E2=X.*Y-6*Z.^2;
E3=(3*X.*Y-8*Z.^2).*Z;
E4=3*(X.*Y-Z.^2).*Z.^2;
E5=X.*Y.*Z.^3;
% This is the elliptic integral approximation
RD=(4^(-n))*(A.^(-3/2)).*(1-3*E2/14+E3/6+9*(E2.^2)/88-3*E4/22-9*E2.*E3/52+3*E5/26)+3*sum;




