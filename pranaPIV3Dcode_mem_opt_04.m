function pranaPIV3Dcode_mem_opt_04(Data);
% This function runs the 3D version of Prana using the data in the
% structure Data as the parameters.  This function and the structure Data
% are created by the function 'pranaPIV3D_parallel_code_runner'.
%
% This version has better memory management.  Additionally any functions 
% that are not called have been deleted from the version.  Any switches 
% that are not programmed in yet in 3D will return errors if called.
%
% PIVcubed has been optimized significantly.  In this version as much data
% is precomputed as possible and the fft calls are optimized with the fftw
% library.  There is a potential bug in PIVcubed when the indices of the
% windows are calculated; this potential bug is documented in PIVcubed.
% Additionally there is a bug that incorrectly saves multiple correlation
% peaks if requested by the user.  The "correct" peak is saved to the
% output array 'C' but the other peaks are saved as NaN values currently.
%
% Authors: Brady Drew, Sam Raben, Rod La Foy, et al . . .
% Last Modified: 17 July 2012 by Rod La Foy

% This initializes the processing
fprintf('\n-------------- Processing Dataset ------------------\n');
pranaprocessing(Data);
fprintf('---------------- Job Completed ---------------------\n');



function pranaprocessing(Data,I1,I2,maskname)
%% --- Read Formatted Parameters ---
%input/output directory
if ispc
    imbase=[Data.imdirec '\' Data.imbase];
    maskbase=[Data.maskdirec '\' Data.maskbase];
    pltdirec=[Data.outdirec '\'];
else
    imbase=[Data.imdirec '/' Data.imbase];
    maskbase=[Data.maskdirec '/' Data.maskbase];
    pltdirec=[Data.outdirec '/'];
end

if nargin<3
    I1 = str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend);
    I2 = I1+str2double(Data.imcstep);
end

%processing mask
if strcmp(Data.masktype,'none')
    itemp = load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(1))],'I');
    mask = 1+0*itemp.I;
    maskname=[];
    clear itemp
elseif strcmp(Data.masktype,'static')
    mask = double(imread(Data.staticmaskname));
    mask = flipud(mask);
    maskname=[];
elseif strcmp(Data.masktype,'dynamic')
    if nargin<4
        maskfend=str2double(Data.maskfstart)+str2double(Data.maskfstep)*length(str2double(Data.imfstart):str2double(Data.imfstep):str2double(Data.imfend))-1;
        maskname=str2double(Data.maskfstart):str2double(Data.maskfstep):maskfend;
    end
end

%method and passes
P=str2double(Data.passes);
Method={'Multipass','Multigrid','Deform','Ensemble','Multiframe'};
M=Method(str2double(Data.method));

%algorithm options
Velinterp=str2double(Data.velinterp);
Iminterp=str2double(Data.iminterp);
Nmax=str2double(Data.framestep);
ds=str2double(Data.PIVerror);

%physical parameters
Mag = str2double(Data.wrmag);
dt = str2double(Data.wrsep);
Freq = str2double(Data.wrsamp);

%initialization
Wres=zeros(P,3);
Wsize=zeros(P,3);
Gres=zeros(P,3);
Gbuf=zeros(P,3);
Corr=zeros(P,1);
D=zeros(P,1);
Zeromean=zeros(P,1);
Peaklocator=zeros(P,1);
Velsmoothswitch=zeros(P,1);
Velsmoothfilt=zeros(P,1);
Valswitch=zeros(P,1);
UODswitch=zeros(P,1);
Bootswitch=zeros(P,1);
Threshswitch=zeros(P,1);
Writeswitch=zeros(P,1);
Peakswitch=zeros(P,1);
UODwinsize=zeros(P,3,1);
UODthresh=zeros(P,1);
Bootper=zeros(P,1);
Bootiter=zeros(P,1);
Bootkmax=zeros(P,1);

Uthresh=zeros(P,2);
Vthresh=zeros(P,2);
Wthresh=zeros(P,2);

extrapeaks=zeros(P,1);
PeakNum=zeros(P,1);
PeakMag=zeros(P,1);
PeakVel=zeros(P,1);
wbase=cell(0);

%read data info for each pass
for e=1:P
    
    %create structure for pass "e"
    eval(['A=Data.PIV' num2str(e) ';'])
    
    %store bulk window offset info
    if e==1
        strfindloc = strfind(A.BWO,',');
        BWO=[str2double(A.BWO(1:strfindloc(1)-1)) str2double(A.BWO(strfindloc(1)+1:strfindloc(2)-1)) str2double(A.BWO(strfindloc(2)+1:end))];
    end
    
    %window and grid resolution
    strfindloc = strfind(A.winres,',');
    Wres(e,:)=[str2double(A.winres(1:strfindloc(1)-1)) str2double(A.winres(strfindloc(1)+1:strfindloc(2)-1)) str2double(A.winres(strfindloc(2)+1:end))];
    strfindloc = strfind(A.winsize,',');
    Wsize(e,:)=[str2double(A.winsize(1:strfindloc(1)-1)) str2double(A.winsize(strfindloc(1)+1:strfindloc(2)-1)) str2double(A.winsize(strfindloc(2)+1:end))];
    strfindloc = strfind(A.gridres,',');
    Gres(e,:)=[str2double(A.gridres(1:strfindloc(1)-1)) str2double(A.gridres(strfindloc(1)+1:strfindloc(2)-1)) str2double(A.gridres(strfindloc(2)+1:end))];
    strfindloc = strfind(A.gridbuf,',');
    Gbuf(e,:)=[str2double(A.gridbuf(1:strfindloc(1)-1)) str2double(A.gridbuf(strfindloc(1)+1:strfindloc(2)-1)) str2double(A.gridbuf(strfindloc(2)+1:end))];
    Corr(e)=str2double(A.corr)-1;
    D(e)=str2double(A.RPCd);
    Zeromean(e)=str2double(A.zeromean);
    Peaklocator(e)=str2double(A.peaklocator);
    Velsmoothswitch(e)=str2double(A.velsmooth);
    Velsmoothfilt(e)=str2double(A.velsmoothfilt);
    
    %validation and thresholding
    Valswitch(e)=str2double(A.val);
    UODswitch(e)=str2double(A.uod);
    Bootswitch(e)=str2double(A.bootstrap);
    Threshswitch(e)=str2double(A.thresh);
    Writeswitch(e)=str2double(A.write);
    
    vpass=[0 strfind(A.uod_window,';') length(A.uod_window)+1];
    for q=1:(length(vpass)-1)
        B=A.uod_window((vpass(q)+1):(vpass(q+1)-1));
        strfindloc = strfind(B,',');
        UODwinsize(e,:,q)=[str2double(B(1:strfindloc(1)-1)) str2double(B(strfindloc(1)+1:strfindloc(2)-1))  str2double(B(strfindloc(2)+1:end))];
        UODthresh(e,q)=str2double(A.uod_thresh(1+2*(q-1)));
    end
    
    Bootper(e)=str2double(A.bootstrap_percentsampled);
    Bootiter(e)=str2double(A.bootstrap_iterations);
    Bootkmax(e)=str2double(A.bootstrap_passes);
    
    if str2double(A.thresh)==1
        Uthresh(e,:)=[str2double(A.valuthresh(1:(strfind(A.valuthresh,',')-1))) str2double(A.valuthresh((strfind(A.valuthresh,',')+1):end))];
        Vthresh(e,:)=[str2double(A.valvthresh(1:(strfind(A.valvthresh,',')-1))) str2double(A.valvthresh((strfind(A.valvthresh,',')+1):end))];
        Wthresh(e,:)=[str2double(A.valwthresh(1:(strfind(A.valwthresh,',')-1))) str2double(A.valwthresh((strfind(A.valwthresh,',')+1):end))];
    else
        Uthresh(e,:)=[-inf,inf];
        Vthresh(e,:)=[-inf,inf];
        Wthresh(e,:)=[-inf,inf];
    end
    
    extrapeaks(e)=str2double(A.valextrapeaks);
    
    %peak information
    Peakswitch(e)=str2double(A.savepeakinfo);
    PeakNum(e)=str2double(A.corrpeaknum);
    PeakMag(e)=str2double(A.savepeakmag);
    PeakVel(e)=str2double(A.savepeakvel);
    
    %output directory
    wbase(e,:)={A.outbase};
    
end


clear('mask');

%% --- Evaluate Image Sequence ---
switch char(M)
    
    case {'Multipass','Multigrid','Deform'}
        
        for q=1:length(I1)
            tf=tic;
            frametitle=['Frame' sprintf(['%0.' Data.imzeros 'i'],I1(q)) ' and Frame' sprintf(['%0.' Data.imzeros 'i'],I2(q))];
            
            %load image pair and flip coordinates
            switch lower(Data.imext)
                case 'tif'
                    im1=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]));
                    im2=double(imread([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]));
                    im1=flipud(im1);
                    im2=flipud(im2);
                case 'mat'
                    it  = load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I1(q))]);
                    im1 = double(it.I(end:-1:1,:,:));
                    it  = load([imbase sprintf(['%0.' Data.imzeros 'i.' Data.imext],I2(q))]);
                    im2 = double(it.I(end:-1:1,:,:));
                    clear('it');
            end
            
            % This clears 'it' from memory
            clear('it');
            
            L=size(im1);
            
            %load dynamic mask and flip coordinates
            if strcmp(Data.masktype,'dynamic')
                mask = double(imread([maskbase sprintf(['%0.' Data.maskzeros 'i.' Data.maskext],maskname(q))]));
                mask = flipud(mask);
            end
            
            for e=1:P
                t1=tic;
                [X,Y,Z]=IMgrid(L,Gres(e,:),Gbuf(e,:));
            
                % This dimensions of the array are accessed individually in
                % case any of the dimensions are singular
                S(1)=size(X,1);
                S(2)=size(X,2);
                S(3)=size(X,3);
                
                X=X(:);Y=Y(:);Z=Z(:);
                
                if strcmp(M,'Multipass')
                    
                    Ub=BWO(1)*ones(prod(L),1);
                    Vb=BWO(2)*ones(prod(L),1);
                    Wb=BWO(3)*ones(prod(L),1);
                    
                elseif e==1
                    
                    UI = BWO(1)*ones(L);
                    Ub = reshape(UI(Y(1):Gres(e,2):Y(end),X(1):Gres(e,1):X(end),Z(1):Gres(e,3):Z(end)),length(X),1);
                    clear('UI');
                    
                    VI = BWO(2)*ones(L);
                    Vb = reshape(VI(Y(1):Gres(e,2):Y(end),X(1):Gres(e,1):X(end),Z(1):Gres(e,3):Z(end)),length(X),1);
                    clear('VI');
                    
                    WI = BWO(3)*ones(L);
                    Wb = reshape(WI(Y(1):Gres(e,2):Y(end),X(1):Gres(e,1):X(end),Z(1):Gres(e,3):Z(end)),length(X),1);
                    clear('WI');
                    
                else
                    
                    Ub=BWO(1)*ones(prod(L),1);
                    Vb=BWO(2)*ones(prod(L),1);
                    Wb=BWO(3)*ones(prod(L),1);
                    
                end
                
                mask=ones(L);
                Eval=reshape(mask(Y(1):Gres(e,2):Y(end),X(1):Gres(e,1):X(end),Z(1):Gres(e,3):Z(end)),length(X),1);
                Eval(Eval==0)=-1;
                Eval(Eval>0)=0;
                clear('mask');
                
                
                %correlate image pair
                if (e~=1) && strcmp(M,'Deform')         %then don't offset windows, images already deformed
                    fprintf('******************************************\nImage Deformation not yet programed in 3D\n******************************************\n');
                    error(' ');
                    
                else                                    %either first pass, or not deform
                    
                    if Corr(e)<2
                        disp(' ');
                        disp(['Processing correlations for pass ',num2str(e),'/',num2str(P),'.']);
                        [~,~,~,Uc,Vc,Wc,Cc]=PIVcubed(im1,im2,Corr(e),Wsize(e,:),Wres(e,:),0,D(e),Zeromean(e),Peaklocator(e),Peakswitch(e) || (Valswitch(e) && extrapeaks(e)),X(Eval>=0),Y(Eval>=0),Z(Eval>=0),Ub(Eval>=0),Vb(Eval>=0),Wb(Eval>=0));
                        disp(['Completed pass ',num2str(e),'/',num2str(P),'.']);
                    else
                        fprintf('******************************************\nPhase Correlation not yet programed in 3D\n******************************************\n')
                        error(' ');
                    end
                end
                
                if Corr(e)<2
                    if Peakswitch(e) || (Valswitch(e) && extrapeaks(e))
                        U=zeros(size(X,1),3);
                        V=zeros(size(X,1),3);
                        W=zeros(size(X,1),3);
                        U(repmat(Eval>=0,[1 3]))=Uc;
                        V(repmat(Eval>=0,[1 3]))=Vc;
                        W(repmat(Eval>=0,[1 3]))=Wc;
                        C=zeros(size(X,1),3); C(repmat(Eval>=0,[1 3]))=Cc;
                    else
                        U=zeros(size(X));
                        V=zeros(size(X));
                        W=zeros(size(X));
                        C=[];
                        U(Eval>=0)=Uc;
                        V(Eval>=0)=Vc;
                        W(Eval>=0)=Wc;
                    end
                else
                    U=zeros(size(X));
                    V=zeros(size(X));
                    W=zeros(size(X));
                    U(Eval>=0)=Uc;
                    V(Eval>=0)=Vc;
                    W(Eval>=0)=Wc;
                    if Peakswitch(e)
                        C=zeros(size(X,1),3);
                        C(repmat(Eval>=0,[1 3]))=Cc;
                    else
                        C=[];
                    end
                end
                
                corrtime(e)=toc(t1);
                
                %validation
                if Valswitch(e)
                    t1=tic;
                    
                    [Uval,Vval,Wval,Evalval,Cval]=VAL(X,Y,Z,U,V,W,Eval,C,Threshswitch(e),UODswitch(e),Bootswitch(e),extrapeaks(e),...
                        Uthresh(e,:),Vthresh(e,:),Wthresh(e,:),UODwinsize(e,:,:),UODthresh(e,UODthresh(e,:)~=0)',Bootper(e),Bootiter(e),Bootkmax(e));
                    
                    valtime(e)=toc(t1);
                else
                    Uval=U(:,1);Vval=V(:,1);Wval=W(:,1);Evalval=Eval(:,1);
                    if ~isempty(C)
                        Cval=C(:,1);
                    else
                        Cval=[];
                    end
                end
                
                %write output
                if Writeswitch(e)
                    t1=tic;
                    
                    if Peakswitch(e)
                        if PeakVel(e) && Corr(e)<2
                            U=[Uval,U(:,1:PeakNum(e))];
                            V=[Vval,V(:,1:PeakNum(e))];
                            W=[Wval,W(:,1:PeakNum(e))];
                        else
                            U=Uval; V=Vval; W=Wval;
                        end
                        if PeakMag(e)
                            C=[Cval,C(:,1:PeakNum(e))];
                        else
                            C=Cval;
                        end
                    else
                        U=Uval; V=Vval; W=Wval; C=Cval;
                    end
                    Eval=Evalval;
                    
                    %convert to physical units
                    Xval=X;Yval=Y;Zval=Z;
                    X=X*Mag;Y=Y*Mag;Z=Z*Mag;
                    U=U*Mag/dt;V=V*Mag/dt;W=W*Mag/dt;
                    
                    %convert to matrix if necessary
                    if size(X,2)==1
                        [X,Y,Z,U,V,W,Eval,C]=matrixform(X,Y,Z,U,V,W,Eval,C);
                    end
                    
                    %remove nans from data, replace with zeros
                    U(Eval<0)=0;V(Eval<0)=0;W(Eval<0)=0;
                    
                    if str2double(Data.datout)
                        time=I1(q)/Freq;
                        write_dat_val_C([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.dat' ],I1(q))],X,Y,Z,U,V,W,Eval,C,e,time,frametitle);
                    end
                    
                    if str2double(Data.multiplematout)
                        save([pltdirec char(wbase(e,:)) sprintf(['%0.' Data.imzeros 'i.mat' ],I1(q))],'X','Y','Z','U','V','W','Eval','C')
                    end
                    X=Xval;Y=Yval;Z=Zval;
                    
                    savetime(e)=toc(t1);
                end
                U=Uval; V=Vval; W=Wval;
                
                if e~=P
                    %reshape from list of grid points to matrix
                    X=reshape(X,[S(1),S(2),S(3)]);
                    Y=reshape(Y,[S(1),S(2),S(3)]);
                    Z=reshape(Z,[S(1),S(2),S(3)]);
                    U=reshape(U(:,1),[S(1),S(2),S(3)]);
                    V=reshape(V(:,1),[S(1),S(2),S(3)]);
                    W=reshape(W(:,1),[S(1),S(2),S(3)]);
                    
                    if strcmp(M,'Multigrid') || strcmp(M,'Deform')
                        t1=tic;
                        
                        %velocity smoothing
                        if Velsmoothswitch(e)==1
                            [U,V,W]=VELfilt(U,V,W,Velsmoothfilt(e));
                        end
                        
                        %velocity interpolation
                        [Xvf,Yvf,Zvf]=IMgrid2(L,Gres(e+1,:),Gbuf(e+1,:)); %
                        UI = VFinterp(X,Y,Z,U,Xvf,Yvf,Zvf,Velinterp);
                        UI = UI(2:end-1,2:end-1,2:end-1);
                        VI = VFinterp(X,Y,Z,V,Xvf,Yvf,Zvf,Velinterp);
                        VI = VI(2:end-1,2:end-1,2:end-1);
                        WI = VFinterp(X,Y,Z,W,Xvf,Yvf,Zvf,Velinterp);
                        WI = WI(2:end-1,2:end-1,2:end-1);
                        interptime(e)=toc(t1);
                        
                        if strcmp(M,'Deform')
                            fprintf('******************************************\nImage Deformation not yet programed in 3D\nSwitching to Mutligrid\n******************************************\n')
                            error(' ');
                        end
                    else
                        UI=U;VI=V;WI=W;
                    end
                end
            end
            
            eltime=toc(tf);
            %output text
            fprintf('\n----------------------------------------------------\n')
            fprintf(['Job: ',Data.batchname,'\n'])
            fprintf([frametitle ' Completed (' num2str(q) '/' num2str(length(I1)) ')\n'])
            fprintf('----------------------------------------------------\n')
            for e=1:P
                fprintf('correlation...                   %0.2i:%0.2i.%0.0f\n',floor(corrtime(e)/60),floor(rem(corrtime(e),60)),rem(corrtime(e),60)-floor(rem(corrtime(e),60)))
                if Valswitch(e)
                    fprintf('validation...                    %0.2i:%0.2i.%0.0f\n',floor(valtime(e)/60),floor(rem(valtime(e),60)),rem(valtime(e),60)-floor(rem(valtime(e),60)))
                end
                if Writeswitch(e)
                    fprintf('save time...                     %0.2i:%0.2i.%0.0f\n',floor(savetime(e)/60),floor(rem(savetime(e),60)),rem(savetime(e),60)-floor(rem(savetime(e),60)))
                end
                if strcmp(M,'Multigrid') || strcmp(M,'Deform')
                    if e~=P
                        fprintf('velocity interpolation...        %0.2i:%0.2i.%0.0f\n',floor(interptime(e)/60),floor(rem(interptime(e),60)),rem(interptime(e),60)-floor(rem(interptime(e),60)))
                        if strcmp(M,'Deform')
                            fprintf('image deformation...             %0.2i:%0.2i.%0.0f\n',floor(deformtime(e)/60),floor(rem(deformtime(e),60)),rem(deformtime(e),60)-floor(rem(deformtime(e),60)))
                        end
                    end
                end
            end
            fprintf('total frame time...              %0.2i:%0.2i.%0.0f\n',floor(eltime/60),floor(rem(eltime,60)),rem(eltime,60)-floor(rem(eltime,60)))
            frametime(q)=eltime;
            comptime=mean(frametime)*(length(I1)-q);
            fprintf('estimated job completion time... %0.2i:%0.2i:%0.2i\n\n',floor(comptime/3600),floor(rem(comptime,3600)/60),floor(rem(comptime,60)))
        end
        
    case 'Ensemble'
        fprintf('******************************************\nEnsemble Correlation not yet programed in 3D\n******************************************\n')
        error(' ');
        
    case 'Multiframe'
        fprintf('******************************************\nMultiframe Correlation not yet programed in 3D\n******************************************\n')
        error(' ');
end

%signal job complete
beep,pause(0.2),beep



function [X,Y,Z,U,V,W,C]=PIVcubed(im1,im2,corr,window,res,zpad,D,Zeromean,Peaklocator,Peakswitch,X,Y,Z,Uin,Vin,Win)
% --- DPIV Correlation ---

% This is the size of the image
L=size(im1);

% This normalizes the images (in case they are very large or very small
% values in the images due to prior scaling - this avoids numerical errors
% in the correlation)
im_min=min([min(im1(:)),min(im2(:))]);
im_max=max([max(im1(:)),max(im2(:))]);
im1=(im1-im_min)/(im_max-im_min);
im2=(im2-im_min)/(im_max-im_min);

%convert to gridpoint list
X=X(:);
Y=Y(:);
Z=Z(:);

%correlation and window mask types
ctype    = {'SCC','RPC'};
tcorr = char(ctype(corr+1));

%preallocate velocity fields and grid format
Nx = window(1);
Ny = window(2);
Nz = window(3);
if nargin <= 14
    Uin = zeros(length(X),1);
    Vin = zeros(length(X),1);
    Win = zeros(length(X),1);
end

if Peakswitch
    Uin=repmat(Uin(:,1),[1 3]);
    Vin=repmat(Vin(:,1),[1 3]);
    Win=repmat(Win(:,1),[1 3]);
    U = zeros(length(X),3);
    V = zeros(length(X),3);
    W = zeros(length(X),3);
    C = zeros(length(X),3);
else
    U = zeros(length(X),1);
    V = zeros(length(X),1);
    W = zeros(length(X),1);
    C = [];
end

%sets up extended domain size
if zpad~=0
    Sx=2*Nx;
    Sy=2*Ny;
    Sz=2*Nz;
else
    Sx=Nx;
    Sy=Ny;
    Sz=Nz;
end

%fftshift indicies
fftindy = [Sy/2+1:Sy 1:Sy/2];
fftindx = [Sx/2+1:Sx 1:Sx/2];
fftindz = [Sz/2+1:Sz 1:Sz/2];

% These are the centers of the windows (including the discrete window
% offset) . . .
%
% If Data.PIVn.savepeakinfo is set to '0' then the size of Uin, Vin, Win
% will be [length(X),1].  If Data.PIVn.savepeakinfo is set to '1' then the
% size of Uin, Vin, Win will be [length(X),3].  I do not know what the
% difference is, but in the old versions of the code
% pranaPIV3Dcode_mem_opt_01 and pranaPIV3Dcode_mem_opt_02 (which came from 
% pranaPIV3Dcode) only the value from the first column is used.  This 
% convention is retained here, although it may not be correct . . .
x1 = X - floor(round(Uin(:,1))/2);
x2 = X +  ceil(round(Uin(:,1))/2);
y1 = Y - floor(round(Vin(:,1))/2);
y2 = Y +  ceil(round(Vin(:,1))/2);
z1 = Z - floor(round(Win(:,1))/2);
z2 = Z +  ceil(round(Win(:,1))/2);

% These are the ranges of the windows (which may extend beyond the edge of
% the image)
xmin1 = x1-Nx/2+1;
xmax1 = x1+Nx/2;
xmin2 = x2-Nx/2+1;
xmax2 = x2+Nx/2;
ymin1 = y1-Ny/2+1;
ymax1 = y1+Ny/2;
ymin2 = y2-Ny/2+1;
ymax2 = y2+Ny/2;
zmin1 = z1-Nz/2+1;
zmax1 = z1+Nz/2;
zmin2 = z2-Nz/2+1;
zmax2 = z2+Nz/2;

% These takes only the section of the windows that are inside the range of
% the first image
ii_min_1=bsxfun(@max,1,ymin1);
ii_max_1=bsxfun(@min,L(1),ymax1);
jj_min_1=bsxfun(@max,1,xmin1);
jj_max_1=bsxfun(@min,L(2),xmax1);
kk_min_1=bsxfun(@max,1,zmin1);
kk_max_1=bsxfun(@min,L(3),zmax1);

% These takes only the section of the windows that are inside the range of
% the second image
ii_min_2=bsxfun(@max,1,ymin2);
ii_max_2=bsxfun(@min,L(1),ymax2);
jj_min_2=bsxfun(@max,1,xmin2);
jj_max_2=bsxfun(@min,L(2),xmax2);
kk_min_2=bsxfun(@max,1,zmin2);
kk_max_2=bsxfun(@min,L(3),zmax2);

% This sets the fft library for optimal speed calculation
fftw('planner','measure');

%%
switch upper(tcorr)

    %Standard Cross Correlation
    case 'SCC'
        
        for n=1:length(X);
            
            % This displays the current vector being processed
            if mod(n,100)==0;
                disp(['Calculating vector ',num2str(n),'/',num2str(length(X)),'.']);
            end;
            
            % This is the n-th window to process from the first image
            zone_1=im1(ii_min_1(n):ii_max_1(n),jj_min_1(n):jj_max_1(n),kk_min_1(n):kk_max_1(n));
            % This does some zero-padding magic
            if size(zone_1,1)~=Ny || size(zone_1,2)~=Nx || size(zone_1,3)~=Nz;
                w1 = zeros(Ny,Nx,Nz);
                w1( 1+max([0 1-ymin1(n)]):Ny-max([0 ymax1(n)-L(1)]),1+max([0 1-xmin1(n)]):Nx-max([0 xmax1(n)-L(2)]),...
                    1+max([0 1-zmin1(n)]):Nz-max([0 zmax1(n)-L(3)]) ) = zone_1;
                zone_1 = w1;
            end;
            % This clears the array 'w1' from memory
            w1=[];
            
            % This is the n-th window to process from the second image
            zone_2=im2(ii_min_2(n):ii_max_2(n),jj_min_2(n):jj_max_2(n),kk_min_2(n):kk_max_2(n));
            % This does some zero-padding magic
            if size(zone_2,1)~=Ny || size(zone_2,2)~=Nx || size(zone_2,3)~=Nz;
                w2 = zeros(Ny,Nx,Nz);
                w2( 1+max([0 1-ymin2(n)]):Ny-max([0 ymax2(n)-L(1)]),1+max([0 1-xmin2(n)]):Nx-max([0 xmax2(n)-L(2)]),...
                    1+max([0 1-zmin2(n)]):Nz-max([0 zmax2(n)-L(3)]) ) = zone_2;
                zone_2 = w2;
            end;
            % This clears the array 'w2' from memory
            w2=[];
            
            % If called for, this subtracts the mean from the images
            if Zeromean==1;
                zone_1=zone_1-mean(zone_1(:));
                zone_2=zone_2-mean(zone_2(:));
            end;

            %apply the image spatial filter
            region1 = (zone_1).*windowmask([Nx Ny Nz],[res(1) res(2) res(3)]);
            % This clears 'zone1' from memory
            zone_1=[];
            % This applies the image spatial filter
            region2 = (zone_2).*windowmask([Nx Ny Nz],[res(1) res(2) res(3)]);
            % This clears 'zone2' from memory
            zone_2=[];
            
            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx Sz]);
            % This clears 'region1' from memory
            region1=[];
            
            % This applies the fft to image 2
            f2   = fftn(region2,[Sy Sx Sz]);
            % This clears 'region2' from memory
            region2=[];
            
            % This is the cross correlation in Fourier space
            P21  = f2.*conj(f1);
            
            % This clears the Fourier transforms from memory
            f1=[];
            f2=[];

            %Standard Fourier Based Cross-Correlation
            G = ifftn(P21,'symmetric');
            % This clears the cross correlation from memory
            P21=[];
            % This extracts a window of G
            G = G(fftindy,fftindx,fftindz);
            % This is the absolute value of the matrix (or the actual
            % cross-correlation function)
            G = abs(G);
            
            %subpixel estimation
            [U(n,:),V(n,:),W(n,:),Ctemp]=subpixel(G,Nx,Ny,Nz,ones(Ny,Nx,Nz),Peaklocator,Peakswitch);
            
            % This clears the cross-correlation from memory
            G=[];
            
            % This saves the peak values
            if Peakswitch;
                C(n,:)=Ctemp;
            end;
        end;

    %Robust Phase Correlation
    case 'RPC'
        
        % This creates the specral RPC filter
        H_Spectral=fftshift(double(energyfilt_fast(Sx,Sy,Sz,D,0)));

        for n=1:length(X);
            
            % This displays the current vector being processed
            if mod(n,100)==0;
                disp(['Calculating vector ',num2str(n),'/',num2str(length(X)),'.']);
            end;
            
            % This is the n-th window to process from the first image
            zone_1=im1(ii_min_1(n):ii_max_1(n),jj_min_1(n):jj_max_1(n),kk_min_1(n):kk_max_1(n));
            % This does some zero-padding magic
            if size(zone_1,1)~=Ny || size(zone_1,2)~=Nx || size(zone_1,3)~=Nz;
                w1 = zeros(Ny,Nx,Nz);
                w1( 1+max([0 1-ymin1(n)]):Ny-max([0 ymax1(n)-L(1)]),1+max([0 1-xmin1(n)]):Nx-max([0 xmax1(n)-L(2)]),...
                    1+max([0 1-zmin1(n)]):Nz-max([0 zmax1(n)-L(3)]) ) = zone_1;
                zone_1 = w1;
            end;
            % This clears the array 'w1' from memory
            w1=[];
            
            % This is the n-th window to process from the second image
            zone_2=im2(ii_min_2(n):ii_max_2(n),jj_min_2(n):jj_max_2(n),kk_min_2(n):kk_max_2(n));
            % This does some zero-padding magic
            if size(zone_2,1)~=Ny || size(zone_2,2)~=Nx || size(zone_2,3)~=Nz;
                w2 = zeros(Ny,Nx,Nz);
                w2( 1+max([0 1-ymin2(n)]):Ny-max([0 ymax2(n)-L(1)]),1+max([0 1-xmin2(n)]):Nx-max([0 xmax2(n)-L(2)]),...
                    1+max([0 1-zmin2(n)]):Nz-max([0 zmax2(n)-L(3)]) ) = zone_2;
                zone_2 = w2;
            end;
            % This clears the array 'w2' from memory
            w2=[];
            
            % If called for, this subtracts the mean from the images
            if Zeromean==1;
                zone_1=zone_1-mean(zone_1(:));
                zone_2=zone_2-mean(zone_2(:));
            end;

            %apply the image spatial filter
            region1 = (zone_1).*windowmask([Nx Ny Nz],[res(1) res(2) res(3)]);
            % This clears 'zone1' from memory
            zone_1=[];
            % This applies the spatial filter
            region2 = (zone_2).*windowmask([Nx Ny Nz],[res(1) res(2) res(3)]);
            % This clears 'zone2' from memory
            zone_2=[];
            
            %FFTs and Cross-Correlation
            f1   = fftn(region1,[Sy Sx Sz]);
            % This clears 'region1' from memory
            region1=[];
            % This is the FFT of the second image
            f2   = fftn(region2,[Sy Sx Sz]);
            % This clears 'region2' from memory
            region2=[];
            % This is the Fourier space cross-correlaiton
            P21  = f2.*conj(f1);
            % This clears the Fourier transforms from memory
            f1=[];
            f2=[];
           
            %Phase Correlation
            Weden = sqrt(P21.*conj(P21));
            R = P21;
            % These are the non-zero indices of P21
            P21_Nonzero=(P21~=0);
            % This rescales R
            R(P21_Nonzero)=P21(P21_Nonzero)./Weden(P21_Nonzero);
            
            % This clears 'P21', 'P21_Nonzero', and 'Weden' from memory
            P21=[];
            P21_Nonzero=[];
            Weden=[];
            
            % This applies the spectral energy filter
            R_Filt=R.*H_Spectral;
            % This clears 'R' from memory
            R=[];
            % This converts out of Fourier space to give the correlation
            % volume
            G = ifftn(R_Filt,'symmetric');
            % This clears 'R' from memory
            R_Filt=[];
            % This extracts a window of the correlation function
            G = G(fftindy,fftindx,fftindz);
            % This is the absolute value of the correlation function
            G = abs(G);

            %subpixel estimation
            [U(n,:),V(n,:),W(n,:),Ctemp]=subpixel(G,Nx,Ny,Nz,ones(Ny,Nx,Nz),Peaklocator,Peakswitch);

            % This clears 'G' from memory
            G=[];

            % This saves the peak values
            if Peakswitch;
                C(n,:)=Ctemp;
            end;
            
        end;
end;

%add DWO to estimation
U = round(Uin)+U;
V = round(Vin)+V;
W = round(Win)+W;



function [X,Y,Z]=IMgrid(L,S,G)
% --- Grid Generation Subfunction ---

%grid buffer
if nargin<3
    G=[0 0 0 0 0 0];
end

S=[S(2) S(1) S(3)];
G=[G(2) G(1) G(3) L(1)-G(2)+1 L(2)-G(1)+1 L(3)-G(3)+1];

%form grid
if max(S)==0
    %pixel grid
    y=(1:L(1))';
    x=1:L(2);
    z(1,1,:)=1:L(3);
else
    if G(1)==0
        %buffers 1/2 grid spacing
        y=(ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)))';
    else
        %predefined grid buffer
        y=(G(1):S(1):G(4))';
    end
    if G(2)==0
        %buffers 1/2 grid spacing
        x=ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2));
    else
        %predefined grid buffer
        x=(G(2):S(2):G(5));
    end
    if G(3)==0
        %buffers 1/2 grid spacing
        z(1,1,:)=ceil((L(3)-(floor(L(3)/S(3))-2)*S(3))/2):S(3):(L(3)-S(3));
    else
        %predefined grid buffer
        z(1,1,:)=(G(3):S(3):G(6));
    end
end
%vector2matrix conversion
X=x(ones(length(y),1),:,ones(length(z),1));
Y=y(:,ones(1,length(x)),ones(length(z),1));
Z=z(ones(1,length(y)),ones(length(x),1),:);



function [X,Y,Z]=IMgrid2(L,S,G)
% --- Grid Generation Subfunction ---

%grid buffer
if nargin<3
    G=[0 0 0 0 0 0];
end

S=[S(2) S(1) S(3)];
G=[G(2) G(1) G(3) L(1)-G(2)+1 L(2)-G(1)+1 L(3)-G(3)+1];

%form grid
if max(S)==0
    %pixel grid
    y=(1:L(1))';
    x=1:L(2);
    z(1,1,:)=1:L(3);
else
    if G(1)==0
        %buffers 1/2 grid spacing
        y=([1 ceil((L(1)-(floor(L(1)/S(1))-2)*S(1))/2):S(1):(L(1)-S(1)) L(1)])';
    else
        %predefined grid buffer
        y=([1 G(1):S(1):G(4) L(1)])';
    end
    if G(2)==0
        %buffers 1/2 grid spacing
        x=[1 ceil((L(2)-(floor(L(2)/S(2))-2)*S(2))/2):S(2):(L(2)-S(2)) L(2)];
    else
        %predefined grid buffer
        x=([1 G(2):S(2):G(5) L(2)]);
    end
    if G(3)==0
        %buffers 1/2 grid spacing
        z(1,1,:)=[1 ceil((L(3)-(floor(L(3)/S(3))-2)*S(3))/2):S(3):(L(3)-S(3)) L(3)];
    else
        %predefined grid buffer
        z(1,1,:)=([1 G(3):S(3):G(6) L(3)]);
    end
end
%vector2matrix conversion
X=x(ones(length(y),1),:,ones(length(z),1));
Y=y(:,ones(1,length(x)),ones(length(z),1));
Z=z(ones(1,length(y)),ones(length(x),1),:);



function [UI]=VFinterp(X,Y,Z,U,XI,YI,ZI,M)
% --- Velocity Interpolation Subfunction

%find grid sizes
Method={'nearest','linear','cubic'};
L=[max(YI(:)) max(XI(:)) max(ZI(:))];

% This dimensions of the array are accessed individually in case any of the
% dimensions are singular
S(1)=size(X,1);
S(2)=size(X,2);
S(3)=size(X,3);

%buffer matrix with nearest neighbor approximation for image boundaries
Xf = zeros(S+2);
Xf(2:end-1,2:end-1,2:end-1) = X;
Xf(2:end-1,2:end-1,1)   = X(:,:,1);
Xf(2:end-1,2:end-1,end) = X(:,:,1);
Xf(1,2:end-1,:)   = permute(repmat(X(1,:,1),[S(3)+2 1]),[3 2 1]);
Xf(end,2:end-1,:) = permute(repmat(X(1,:,1),[S(3)+2 1]),[3 2 1]);
Xf(:,1,:)   = permute(ones([S(1)+2 S(3)+2]),[1 3 2]);
Xf(:,end,:) = L(2).*permute(ones([S(1)+2 S(3)+2]),[1 3 2]);

Yf = zeros(S+2);
Yf(2:end-1,2:end-1,2:end-1) = Y;
Yf(2:end-1,2:end-1,1)   = Y(:,:,1);
Yf(2:end-1,2:end-1,end) = Y(:,:,1);
Yf(2:end-1,1,:)   = permute(repmat(Y(:,1,1),[1 S(3)+2]),[3 1 2]);
Yf(2:end-1,end,:) = permute(repmat(Y(:,1,1),[1 S(3)+2]),[3 1 2]);
Yf(1,:,:)   = permute(ones([S(2)+2 S(3)+2]),[3 1 2]);
Yf(end,:,:) = L(1).*permute(ones([S(2)+2 S(3)+2]),[3 1 2]);

Zf = zeros(S+2);
Zf(2:end-1,2:end-1,2:end-1) = Z;
Zf(:,:,1)   = ones([S(1)+2 S(2)+2]);
Zf(:,:,end) = L(3).*ones([S(1)+2 S(2)+2]);
Zf(1,2:end-1,2:end-1)   = Z(1,:,:);
Zf(end,2:end-1,2:end-1) = Z(end,:,:);
Zf(2:end-1,1,2:end-1) = Z(:,1,:);
Zf(2:end-1,end,2:end-1) = Z(:,end,:);
Zf([1 end],[1 end],2:end-1) = repmat(Z(1,1,:),[2 2 1]);

Uf = zeros(S+2);
Uf(2:end-1,2:end-1,2:end-1) =  U;
Uf(1,2:end-1,2:end-1)       = (U(1,:,:)   - U(2,:,:))./(Y(1,:,:)-Y(2,:,:)).*(1-Y(2,:,:))+U(1,:,:);
Uf(end,2:end-1,2:end-1)     = (U(end,:,:) - U(end-1,:,:))./(Y(end,:,:)-Y(end-1,:,:)).*(L(1)-Y(end-1,:,:))+U(end,:,:);
Uf(2:end-1,1,2:end-1)       = (U(:,1,:)   - U(:,2,:))./(X(:,1,:)-X(:,2,:)).*(1-X(:,2,:))+U(:,1,:);
Uf(2:end-1,end,2:end-1)     = (U(:,end,:) - U(:,end-1,:))./(X(:,end,:)-X(:,end-1,:)).*(L(2)-X(:,end-1,:))+U(:,end,:);
Uf(2:end-1,2:end-1,1)       = mean(U(:,:,1:2),3);
Uf(1,2:end-1,2:end-1)       = mean(U(1:2,:,:),1);
Uf(2:end-1,1,2:end-1)       = mean(U(:,1:2,:),2);
Uf(2:end-1,2:end-1,end)     = mean(U(:,:,end-1:end),3);
Uf(end,2:end-1,2:end-1)     = mean(U(end-1:end,:,:),1);
Uf(2:end-1,end,2:end-1)     = mean(U(:,end-1:end,:),2);
Uf(1,2:end-1,1)     = mean([U(1:2,:,1);         U(2,:,1)  ],1);
Uf(end,2:end-1,1)   = mean([U(end-1:end,:,1);   U(end,:,1)],1);
Uf(1,2:end-1,end)   = mean([U(1:2,:,end);       U(2,:,end)  ],1);
Uf(end,2:end-1,end) = mean([U(end-1:end,:,end); U(end,:,end)],1);
Uf(2:end-1,1,1)     = mean([U(:,1:2,1),         U(:,2,1)  ],2);
Uf(2:end-1,end,1)   = mean([U(:,end-1:end,1),   U(:,end,1)],2);
Uf(2:end-1,1,end)   = mean([U(:,1:2,end),       U(:,2,end)  ],2);
Uf(2:end-1,end,end) = mean([U(:,end-1:end,end), U(:,end,end)],2);
Uf(1,1,:)     = mean(permute([Uf(2:3,2,:);           Uf(2,2,:)        ],[1 3 2]),1);
Uf(1,end,:)   = mean(permute([Uf(2:3,end-1,:);       Uf(2,end-1,:)    ],[1 3 2]),1);
Uf(end,1,:)   = mean(permute([Uf(end-2:end-1,2,:);   Uf(end-1,2,:)    ],[1 3 2]),1);
Uf(end,end,:) = mean(permute([Uf(end-2:end-1,end,:); Uf(end-1,end-1,:)],[1 3 2]),1);

%velocity interpolation
try;
    % This tries the velocity interpolation with the current method
    UI=interp3(Xf,Yf,Zf,Uf,XI,YI,ZI,char(Method(M)));
catch;
    % If the user specified interpolation method fails, a lower order
    % method is attempted
    if M>=1;
        % This displays a warning telling the user that a lower order
        % method is being attempted
        warning(['Velocity interpolation failed with method ''',char(Method(M)),'''. Switching to interpolation method ''',char(Method(M-1)),''' and re-attempting the interpolation.']);
        try;
            % This tries calculating the interpolation again
            UI=interp3(Xf,Yf,Zf,Uf,XI,YI,ZI,char(Method(M-1)));
        catch;
            % If the lower order interpolation method fails, an even lower order
            % method is attempted
            if M>=2;
                % This displays a warning telling the user that a lower order
                % method is being attempted
                warning(['Velocity interpolation failed with method ''',char(Method(M-1)),'''. Switching to interpolation method ''',char(Method(M-2)),''' and re-attempting the interpolation.']);
                try;
                    % This tries calculating the interpolation again
                    UI=interp3(Xf,Yf,Zf,Uf,XI,YI,ZI,char(Method(M-2)));
                catch;
                    % This displays an error stating that the interpolation
                    % cannot be completed
                    error(['Velocity interpolation failed with method ''',char(Method(M-2)),'''.']);
                end;
            else;
                % This displays an error stating that the interpolation
                % cannot be completed
                error(['Velocity interpolation cannot be completed with a lower order method.']);
            end;
        end;
    else;
        % This displays an error stating that the interpolation
        % cannot be completed
        error(['Velocity interpolation cannot be completed with a lower order method.']);
    end;
end;



function [Uf,Vf,Wf]=VELfilt(U,V,W,C)
% --- Velocity Smoothing Subfunction ---
% This function smooths the velocity field input data by filtering the
% velocity field data with a Guassian type convolution kernal.  If the
% input velocity field contains NaN values, these values will be expanded
% to neighboring indices during the filtering operation.  This issue is
% resolved by interpolating over the NaN values before the filtering
% operation.  This issue may also be resolved on some machines by disabling
% the Intel IPP library with the command:
%
%  iptsetpref('UseIPPL',false);
%
% however the due to variation in machine architecture this command may not
% reliably work.
%
% Authors:              Sam Raben, Rod La Foy
% Last Modified On:     16 July 2012
% Last Modified By:     Rod La Foy

% This creates the 3D Gaussian convolution kernal for smoothing the image
std = C;
siz = ([7 7 7]-1)/2;
[xx,yy,zz] = meshgrid(-siz(2):siz(2),-siz(1):siz(1),-siz(3):siz(3));
arg   = -(xx.*xx + yy.*yy + zz.*zz)/(2*std*std);
A     = exp(arg);
A(A<eps*max(A(:))) = 0;
sumA = sum(A(:));
if sumA ~= 0
    A  = A/sumA;
end

% This is a list of the indices of NaN values in the vector field
nan_indices=isnan(U(:))|isnan(V(:))|isnan(W(:));

% If NaN values exist in the velocity field, this replaces the NaN values
% with interpolated values
if any(nan_indices);
    
    % These are the coordinates of the velocity field (or the indices)
    [X,Y,Z]=meshgrid(1:size(U,2),1:size(U,1),1:size(U,3));
    % These are the coordinates of the NaN values
    X_NaN=X(nan_indices);
    Y_NaN=Y(nan_indices);
    Z_NaN=Z(nan_indices);
    % This removes the values of the coordinates where the NaN values exist
    X(nan_indices)=[];
    Y(nan_indices)=[];
    Z(nan_indices)=[];
    % This creates temporary copies of the velocity field (which can have the
    % NaN values removed)
    U_Temp=U;
    V_Temp=V;
    W_Temp=W;
    % This removes the values of the velocity field where the NaN values exist
    % in the temporary variables
    U_Temp(nan_indices)=[];
    V_Temp(nan_indices)=[];
    W_Temp(nan_indices)=[];
    
    % This creates an interpolation structure for the U velocity
    F_U_Interp=TriScatteredInterp(X',Y',Z',U_Temp','nearest');
    % This interpolates the U velocity over the NaN values
    UI=F_U_Interp(X_NaN,Y_NaN,Z_NaN);
    % This saves the interpolated values to the original U field
    U(nan_indices)=UI;
    
    % This creates an interpolation structure for the V velocity
    F_V_Interp=TriScatteredInterp(X',Y',Z',V_Temp','nearest');
    % This interpolates the V velocity over the NaN values
    VI=F_V_Interp(X_NaN,Y_NaN,Z_NaN);
    % This saves the interpolated values to the original U field
    V(nan_indices)=VI;
    
    % This creates an interpolation structure for the W velocity
    F_W_Interp=TriScatteredInterp(X',Y',Z',W_Temp','nearest');
    % This interpolates the W velocity over the NaN values
    WI=F_W_Interp(X_NaN,Y_NaN,Z_NaN);
    % This saves the interpolated values to the original W field
    W(nan_indices)=WI;
    
end;

% This performs the smoothing operation with the Guassian smoothing kernal
Uf=imfilter(U,A,'replicate');
Vf=imfilter(V,A,'replicate');
Wf=imfilter(W,A,'replicate');



function [W]=windowmask(N,R)
% --- Gaussian Window Mask Subfunction ---

% %generic indices
x  = -1:2/(N(1)-1):1;
y  = (-1:2/(N(2)-1):1)';
z  = -1:2/(N(3)-1):1;
%
% %gaussian window sizes
% px = (1.224*N(1)/R(1))^1.0172;
% py = (1.224*N(2)/R(2))^1.0172;
[px]=findwidth(R(1)/N(1));
[py]=findwidth(R(2)/N(2));
[pz]=findwidth(R(3)/N(3));
%
% %generate 2D window
wx=exp(-px^2.*x.^2/2);
wy=exp(-py^2.*y.^2/2);
wz=exp(-pz^2.*z.^2/2);

Wxy  = repmat(wy*wx,[1 1 length(z)]);
W = zeros(size(Wxy));
for i = 1:length(wz)
    W(:,:,i)=Wxy(:,:,i)*wz(i);
end



function [W]=energyfilt_fast(Nx,Ny,Nz,d,q)
% --- RPC Spectral Filter Subfunction ---
%
% This version is designed to run faster by pre-computing values as well as
% performing fewer operations.
%
% Authors: Adric Eckstein, Sam Raben, Rod La Foy
% Last Modified: 18 June 2012 by Rod La Foy

% This is the bit-depth of the image
I0=2^8-1;

%assume no aliasing
if nargin<4
    q = 0;
end

% These are a some pre-calculated constants for later referencing
%
% This is d^2
d_sqr=d^2;
% This is the constant in front of each exponential
C1=(pi*I0*d_sqr/8)^2;
% This is 2*pi
pi_2=2*pi;
% This is d^2/16
d_sqr_16=d_sqr/16;
% This is pi^2
pi_sqr=pi^2;

%initialize indices
[k1,k2,k3]=ndgrid(-pi:pi_2/Ny:pi-pi_2/Ny,-pi:pi_2/Nx:pi-pi_2/Nx,-pi:pi_2/Nz:pi-pi_2/Nz);

% These are some more pre-calculated constants
%
% This is k1^2+k2^2+k3^2
k_sqr_sum=k1.^2+k2.^2+k3.^2;
% This is 2*pi*(n_i+n_j+n_k)
pi_2_n_sum=pi_2*[6,2,2,-2,4,0,4,0,2,-2,-2,-6,0,-4,0,-4,4,0,0,-4,2,-2,2,-2];
% This is pi^2*(n_i^2+n_j^2+n_k^2)
pi_sqr_n_sqr_sum=pi_sqr*[12,12,12,12,8,8,8,8,12,12,12,12,8,8,8,8,8,8,8,8,4,4,4,4];

%particle-image spectrum
Ep = C1*exp(-d_sqr_16*k_sqr_sum);

%aliased particle-image spectrum
Ea = C1*(exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(01)+pi_sqr_n_sqr_sum(01)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(02)+pi_sqr_n_sqr_sum(02)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(03)+pi_sqr_n_sqr_sum(03)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(04)+pi_sqr_n_sqr_sum(04)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(05)+pi_sqr_n_sqr_sum(05)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(06)+pi_sqr_n_sqr_sum(06)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(07)+pi_sqr_n_sqr_sum(07)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(08)+pi_sqr_n_sqr_sum(08)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(09)+pi_sqr_n_sqr_sum(09)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(10)+pi_sqr_n_sqr_sum(10)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(11)+pi_sqr_n_sqr_sum(11)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(12)+pi_sqr_n_sqr_sum(12)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(13)+pi_sqr_n_sqr_sum(13)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(14)+pi_sqr_n_sqr_sum(14)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(15)+pi_sqr_n_sqr_sum(15)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(16)+pi_sqr_n_sqr_sum(16)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(17)+pi_sqr_n_sqr_sum(17)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(18)+pi_sqr_n_sqr_sum(18)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(19)+pi_sqr_n_sqr_sum(19)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(20)+pi_sqr_n_sqr_sum(20)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(21)+pi_sqr_n_sqr_sum(21)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(22)+pi_sqr_n_sqr_sum(22)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(23)+pi_sqr_n_sqr_sum(23)))+...
         exp(-d_sqr_16*(k_sqr_sum+pi_2_n_sum(24)+pi_sqr_n_sqr_sum(24))));

% This clears the index matrices 'k1', 'k2', and 'k3'
clear('k1','k2','k3');

%noise spectrum
En = pi/4*Nx*Ny;

%DPIV SNR spectral filter
W  = Ep./((1-q)*En+(q)*Ea);

% This clears 'Ep' and 'Ea' from memory
clear('Ep','Ea');

% Adric has a transpose but for symmetric particles this is not nessiary
% I don't know how to do a transpose in 3D
% W  = W'/max(max(W));keyboard
W  = W/max(W(:));



function [u,v,w,M]=subpixel(G,ccsizex,ccsizey,ccsizez,W,Method,Peakswitch)
%intialize indices
cc_x = -ccsizex/2:ccsizex/2-1;
cc_y = -ccsizey/2:ccsizey/2-1;
cc_z = -ccsizez/2:ccsizez/2-1;

%find maximum correlation value
[M,I] = max(G(:));

%if correlation empty
if M==0
    if Peakswitch
        u=zeros(1,3);
        v=zeros(1,3);
        w=zeros(1,3);
        M=zeros(1,3);
    else
        u=0; v=0; w=0;
    end
else
    if Peakswitch
        %Locate peaks using imregionalmax
        if any(isnan(G(:)));
            warning('NaN values found while trying to perform a regional maxima calculation.  Returning zero values and continuing processing.');
            u=0;
            v=0;
            w=0;
            return;
        end;
        A=imregionalmax(G);
        peakmat=G.*A;
        for i=2:3
            peakmat(peakmat==M(i-1))=0;
            [M(i),I(i)]=max(peakmat(:));
        end
        j=length(M);
    else
        j=1;
    end
    
    for i=1:j
        % I hav only set up 3point Gaussian right now.  Will employ the
        % other methods later
        method=Method;
        
        %find x, y and z indices
        %         shift_locy = 1+mod(I(i)-1,ccsizey);
        %         shift_locx = ceil(I(i)/ccsizey);
        [shift_locy, shift_locx, shift_locz] = ind2sub(size(G),I(i));
        
        shift_errx=[];
        shift_erry=[];
        shift_errz=[];
        
        %find subpixel displacement in x
        if shift_locx == 1
            %boundary condition 1
            shift_errx =  G( shift_locy , shift_locx+1 , shift_locz )/M(i); method=1;
        elseif shift_locx == ccsizex
            %boundary condition 2
            shift_errx = -G( shift_locy , shift_locx-1 , shift_locz )/M(i); method=1;
        elseif G( shift_locy , shift_locx+1 ) == 0
            %endpoint discontinuity 1
            shift_errx = -G( shift_locy , shift_locx-1 , shift_locz )/M(i); method=1;
        elseif G( shift_locy , shift_locx-1 ) == 0
            %endpoint discontinuity 2
            shift_errx =  G( shift_locy , shift_locx+1 , shift_locz )/M(i); method=1;
        end
        if shift_locy == 1
            %boundary condition 1
            shift_erry = -G( shift_locy+1 , shift_locx , shift_locz )/M(i); method=1;
        elseif shift_locy == ccsizey
            %boundary condition 2
            shift_erry =  G( shift_locy-1 , shift_locx , shift_locz )/M(i); method=1;
        elseif G( shift_locy+1 , shift_locx , shift_locz ) == 0
            %endpoint discontinuity 1
            shift_erry =  G( shift_locy-1 , shift_locx , shift_locz )/M(i); method=1;
        elseif G( shift_locy-1 , shift_locx , shift_locz ) == 0
            %endpoint discontinuity 2
            shift_erry = -G( shift_locy+1 , shift_locx , shift_locz )/M(i); method=1;
        end
        if shift_locz == 1
            %boundary condition 1
            shift_errz = -G( shift_locy , shift_locx , shift_locz+1 )/M(i); method=1;
        elseif shift_locz == ccsizez
            %boundary condition 2
            shift_errz =  G( shift_locy , shift_locx , shift_locz-1 )/M(i); method=1;
        elseif G( shift_locy , shift_locx , shift_locz+1 ) == 0
            %endpoint discontinuity 1
            shift_errz =  G( shift_locy , shift_locx , shift_locz-1 )/M(i); method=1;
        elseif G( shift_locy , shift_locx , shift_locz-1 ) == 0
            %endpoint discontinuity 2
            shift_errz = -G( shift_locy , shift_locx , shift_locz+1 )/M(i); method=1;
        end
        
        if method==2
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            % Gaussian Least Squares %
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %Find a suitable window around the peak (5x5 preferred)
            x_min=shift_locx-2; x_max=shift_locx+2;
            y_min=shift_locy-2; y_max=shift_locy+2;
            z_min=shift_locz-2; z_max=shift_locz+2;
            if x_min<1
                x_min=1;
            end
            if x_max>ccsizex
                x_max=ccsizex;
            end
            if y_min<1
                y_min=1;
            end
            if y_max>ccsizey
                y_max=ccsizey;
            end
            if z_min<1
                z_min=1;
            end
            if z_max>ccsizez
                z_max=ccsizez;
            end
            
            points=G(y_min:y_max,x_min:x_max,z_min:z_max).*W(y_min:y_max,x_min:x_max,z_min:z_max);
            
            %Options for the lsqnonlin solver
            options=optimset('MaxIter',1200,'MaxFunEvals',5000,'TolX',5e-6,'TolFun',5e-6,...
                'LargeScale','off','Display','off','DiffMinChange',1e-7,'DiffMaxChange',1,...
                'Algorithm','levenberg-marquardt');
            
            %Initial values for the solver
            x0=[M(i) 1 shift_locx shift_locy shift_locz];
            
            [xloc yloc,zloc]=meshgrid(x_min:x_max,y_min:y_max,z_min:z_max);
            
            %Run solver; default to 3-point gauss if it fails
            try
                xvars=lsqnonlin(@leastsquares3D,x0,[],[],options,points(:),[yloc(:),xloc(:),zloc(:)]);
                shift_errx=xvars(3)-shift_locx;
                shift_erry=xvars(4)-shift_locy;
                shift_errz=xvars(5)-shift_locz;
            catch EE
                method=1;
                fprintf('Leastsquares Method Failed\nPerforming 3 Point Gaussian Subpixel Interpolations\n%s',EE.message)
            end
            
        end
        
        if method==1
            
            %%%%%%%%%%%%%%%%%%%%
            % 3-Point Gaussian %
            %%%%%%%%%%%%%%%%%%%%
            
            if isempty(shift_errx)
                %gaussian fit
                lCm1 = log(G( shift_locy , shift_locx-1 , shift_locz )*W( shift_locy , shift_locx-1 , shift_locz ));
                lC00 = log(G( shift_locy , shift_locx   , shift_locz )*W( shift_locy , shift_locx   , shift_locz ));
                lCp1 = log(G( shift_locy , shift_locx+1 , shift_locz )*W( shift_locy , shift_locx+1 , shift_locz ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_errx = 0;
                else
                    shift_errx = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                end
            end
            if isempty(shift_erry)
                lCm1 = log(G( shift_locy-1 , shift_locx , shift_locz )*W( shift_locy-1 , shift_locx , shift_locz ));
                lC00 = log(G( shift_locy   , shift_locx , shift_locz )*W( shift_locy   , shift_locx , shift_locz ));
                lCp1 = log(G( shift_locy+1 , shift_locx , shift_locz )*W( shift_locy+1 , shift_locx , shift_locz ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_erry = 0;
                else
                    shift_erry = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                end
            end
            if isempty(shift_errz)
                lCm1 = log(G( shift_locy , shift_locx , shift_locz-1 )*W( shift_locy , shift_locx , shift_locz-1 ));
                lC00 = log(G( shift_locy , shift_locx , shift_locz   )*W( shift_locy , shift_locx , shift_locz ));
                lCp1 = log(G( shift_locy , shift_locx , shift_locz+1 )*W( shift_locy , shift_locx , shift_locz+1 ));
                if (2*(lCm1+lCp1-2*lC00)) == 0
                    shift_errz = 0;
                else
                    shift_errz = (lCm1-lCp1)/(2*(lCm1+lCp1-2*lC00));
                end
            end
            
        end
        
        u(i)=cc_x(shift_locx)+shift_errx;
        v(i)=cc_y(shift_locy)+shift_erry;
        w(i)=cc_z(shift_locz)+shift_errz;
        
        if isinf(u(i)) || isinf(v(i)) || isinf(w(i))
            u(i)=0; v(i)=0; w(i)=0;
        end
    end
end



function F = leastsquares3D(x,mapint_i,locxyz_i)
%This function is called by lsqnonlin if the least squares or continuous
%least squares method has been chosen. x contains initial guesses[I0, betas, x_c,
%y_c]. mapint_i is a matrix containing pixel intensity values, and locxy_i
%is a 1x2 vector containing the row/column coordinates of the top left
%pixel in mapint_i
%
%F is the variable being minimized - the difference between the gaussian
%curve and the actual intensity values.
%
%Adapted from M. Brady's 'leastsquaresgaussfit' and 'mapintensity'
%B.Drew - 7.18.2008

I0=x(1);
betas=x(2);
x_centroid=x(3);
y_centroid=x(4);
z_centroid=x(5);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Just like in the continuous four-point method, lsqnonlin tries negative
%values for x(2), which will return errors unless the abs() function is
%used in front of all the x(2)'s.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

num1=(I0*pi)/4;
num2=sqrt(abs(betas));

gauss_int = zeros(size(mapint_i));
xp = zeros(size(mapint_i));
yp = zeros(size(mapint_i));
zp = zeros(size(mapint_i));
for ii = 1:length(mapint_i)
    xp(ii) = locxyz_i(ii,2);
    yp(ii) = locxyz_i(ii,1);
    zp(ii) = locxyz_i(ii,3);
end

% map an intensity profile of a gaussian function:
for rr = 1:size(xp,1)
    gauss_int(rr)=I0*exp(-abs(betas)*(((xp(rr))-x_centroid)^2 + ...
        ((yp(rr))-y_centroid)^2 + ((zp(rr))-z_centroid)^2));
end

% compare the Gaussian curve to the actual pixel intensities
F=mapint_i-gauss_int;



function [Uval,Vval,Wval,Evalval,Cval]=VAL(X,Y,Z,U,V,W,Eval,C,Threshswitch,UODswitch,Bootswitch,extrapeaks,Uthresh,Vthresh,Wthresh,UODwinsize,UODthresh,Bootper,Bootiter,Bootkmax)
% --- Validation Subfunction ---
if extrapeaks
    j=3;
else
    j=1;
end

[X,Y,Z,U,V,W,Eval,C]=matrixform(X,Y,Z,U,V,W,Eval,C);
Uval=U(:,:,:,1);Vval=V(:,:,:,1);Wval=W(:,:,:,1);Evalval=Eval(:,:,:,1);
if ~isempty(C)
    Cval=C(:,:,:,1);
else
    Cval=[];
end

% The size of the array is accessed by each index individually in case any
% of the dimensions are singular
S(1)=size(X,1);
S(2)=size(X,2);
S(3)=size(X,3);

if Threshswitch || UODswitch
    for i=1:j
        %Thresholding
        if Threshswitch
            [Uval,Vval,Wval,Evalval] = Thresh(Uval,Vval,Wval,Uthresh,Vthresh,Wthresh,Evalval);
        end
        
        %Univeral Outlier Detection
        if UODswitch
            t=permute(UODwinsize,[2 3 1]);
            t=t(:,t(1,:)~=0);
            [Uval,Vval,Wval,Evalval] = UOD(Uval,Vval,Wval,t',UODthresh,Evalval);
        end
        %         disp([num2str(sum(sum(Evalval>0))),' bad vectors'])
        %Try additional peaks where validation failed
        if i<j
            Utemp=U(:,:,:,i+1);Vtemp=V(:,:,:,i+1);Wtemp=W(:,:,:,i+1);Evaltemp=Eval(:,:,:,i+1);Ctemp=C(:,:,:,i+1);
            Uval(Evalval>0)=Utemp(Evalval>0);
            Vval(Evalval>0)=Vtemp(Evalval>0);
            Wval(Evalval>0)=Wtemp(Evalval>0);
            Evalval(Evalval>0)=Evaltemp(Evalval>0);
            Cval(Evalval>0)=Ctemp(Evalval>0);
        end
    end
end

% .........................................................................
% The commented version is an attempt to break an infinite loop where the
% window size continually grows, but it fails if there are not a
% sufficiently large number of points to interpolate from.
% .........................................................................

% %replacement
% for i=1:S(1)
%     for j=1:S(2)
%         for k=1:S(3)
%             if Evalval(i,j,k)>0 && Evalval(i,j,k)<200
%                 %initialize replacement search size
%                 q=0;
%                 s=0;
%                 
%                 % This initializes the extraction domain of U
%                 U_Domain=zeros(1,6);
%                 
%                 %get replacement block with at least 8 valid points
%                 while s==0
%                     q=q+1;
%                     Imin = max([i-q 1   ]);
%                     Imax = min([i+q S(1)]);
%                     Jmin = max([j-q 1   ]);
%                     Jmax = min([j+q S(2)]);
%                     Kmin = max([k-q 1   ]);
%                     Kmax = min([k+q S(3)]);
%                     
%                     % This is a vector of the current domain
%                     U_Domain_Temp=[Imin,Imax,Jmin,Jmax,Kmin,Kmax];
%                     
%                     % This is the difference between the domain of the
%                     % current iteration and the last iteration
%                     r_domain=sum((U_Domain-U_Domain_Temp).^2);
%                     
%                     % If the domain is equal to the domain from the
%                     % last iteration, then this exits the loop and uses
%                     % the current domain for the UOD calculation
%                     if r_domain<eps(1);
%                         % This extracts the same domains from V and W
%                         Xblock = X(Iind,Jind,Kind)-X(i,j,k);
%                         Yblock = Y(Iind,Jind,Kind)-Y(i,j,k);
%                         Zblock = Z(Iind,Jind,Kind)-Z(i,j,k);
%                         Vblock = Vval(Iind,Jind,Kind);
%                         Wblock = Wval(Iind,Jind,Kind);
%                         % This sets s equal to 1 to exit the loop
%                         s=1;
%                     end;
%                     
%                     
%                     % This sets the temporary current domain vector to
%                     % the non-temporary current domain vector
%                     U_Domain=U_Domain_Temp;
%                     
%                     % These are the indices to extract from U
%                     Iind = Imin:Imax;
%                     Jind = Jmin:Jmax;
%                     Kind = Kmin:Kmax;
%                     % This extracts the current domain
%                     Ublock = Uval(Iind,Jind,Kind);
%                     
%                     if length(Ublock(~isnan(Ublock)))>=8
%                         Xblock = X(Iind,Jind,Kind)-X(i,j,k);
%                         Yblock = Y(Iind,Jind,Kind)-Y(i,j,k);
%                         Zblock = Z(Iind,Jind,Kind)-Z(i,j,k);
%                         Vblock = Vval(Iind,Jind,Kind);
%                         Wblock = Wval(Iind,Jind,Kind);
%                         s=1;
%                     end
%                 end
%                 
%                 %distance from erroneous vector
%                 Dblock = (Xblock.^2+Yblock.^2+Zblock.^2).^0.5;
%                 Dblock(isnan(Ublock))=nan;
%                 
%                 %validated vector
%                 Uval(i,j,k) = nansum(nansum(nansum(Dblock.*Ublock)))/nansum(nansum(nansum(Dblock)));
%                 Vval(i,j,k) = nansum(nansum(nansum(Dblock.*Vblock)))/nansum(nansum(nansum(Dblock)));
%                 Wval(i,j,k) = nansum(nansum(nansum(Dblock.*Wblock)))/nansum(nansum(nansum(Dblock)));
%             end
%         end
%     end
% end

%replacement
for i=1:S(1)
    for j=1:S(2)
        for k=1:S(3)
            if Evalval(i,j,k)>0 && Evalval(i,j,k)<200
                %initialize replacement search size
                q=0;
                s=0;
                
                %get replacement block with at least 8 valid points
                while s==0
                    q=q+1;
                    Imin = max([i-q 1   ]);
                    Imax = min([i+q S(1)]);
                    Jmin = max([j-q 1   ]);
                    Jmax = min([j+q S(2)]);
                    Kmin = max([k-q 1   ]);
                    Kmax = min([k+q S(3)]);
                    Iind = Imin:Imax;
                    Jind = Jmin:Jmax;
                    Kind = Kmin:Kmax;
                    Ublock = Uval(Iind,Jind,Kind);
                    if length(Ublock(~isnan(Ublock)))>=8
                        Xblock = X(Iind,Jind,Kind)-X(i,j,k);
                        Yblock = Y(Iind,Jind,Kind)-Y(i,j,k);
                        Zblock = Z(Iind,Jind,Kind)-Z(i,j,k);
                        Vblock = Vval(Iind,Jind,Kind);
                        Wblock = Wval(Iind,Jind,Kind);
                        s=1;
                    end
                end
                
                %distance from erroneous vector
                Dblock = (Xblock.^2+Yblock.^2+Zblock.^2).^0.5;
                Dblock(isnan(Ublock))=nan;
                
                %validated vector
                Uval(i,j,k) = nansum(nansum(nansum(Dblock.*Ublock)))/nansum(nansum(nansum(Dblock)));
                Vval(i,j,k) = nansum(nansum(nansum(Dblock.*Vblock)))/nansum(nansum(nansum(Dblock)));
                Wval(i,j,k) = nansum(nansum(nansum(Dblock.*Wblock)))/nansum(nansum(nansum(Dblock)));
            end
        end
    end
end

%Bootstrapping
if Bootswitch
    [Uval,Vval,Wval,Evalval] = bootstrapping(X,Y,Z,Uval,Vval,Wval,Bootper/100,Bootiter,Bootkmax,Evalval);
end

%convert back to vector
[Uval,Vval,Wval,Evalval,Cval]=vectorform(X,Y,Z,Uval,Vval,Wval,Evalval,Cval);



function [U,V,W,Eval] = Thresh(U,V,W,uthreshold,vthreshold,wthreshold,Eval)
% --- Thresholding Validation Subfunction ---

%neglect u and v threshold
if nargin<=4
    uthreshold = [-inf inf];
    vthreshold = [-inf inf];
    wthreshold = [-inf inf];
end

S(1)=size(U,1);
S(2)=size(U,2);
S(3)=size(U,3);

%thresholding
for i=1:S(1)
    for j=1:S(2)
        for k = 1:S(3)
            if Eval(i,j,k)==0
                %velocity threshold condition
                if U(i,j,k)<uthreshold(1) || U(i,j,k)>uthreshold(2) || V(i,j,k)<vthreshold(1) || V(i,j,k)>vthreshold(2) || W(i,j,k)<wthreshold(1) || W(i,j,k)>wthreshold(2)
                    U(i,j,k)=nan;
                    V(i,j,k)=nan;
                    W(i,j,k)=nan;
                    Eval(i,j,k)=100;
                end
            elseif Eval(i,j,k)==-1
                %boundary condition
                U(i,j,k)=nan;
                V(i,j,k)=nan;
                W(i,j,k)=nan;
            end
        end
    end
end



function [U,V,W,Eval] = UOD(U,V,W,t,tol,Eval)
% --- Universal Outlier Detection Validation Subfunction ---
%
% This function performs the Universal Outlier Detection algorithm on the
% 3D velocity fields denoted by U, V, W.  The output matrix Eval contains
% nonzero values at the points where an outlier was detected.
%
% In this version a break was inserted to exit infinite loops in the while
% structure.  It is uncertain how this will affect the performance of the
% UOD function, however it appears that it is "rare" that an infinite loop
% is encountered.
%
% Authors: Rod La Foy and others . . .
% Created: . . . 
% Last Modified: 16 April 2013

%number of validation passes
pass = length(tol);

% This is the size of the array to process (the indices are accessed
% individually in case any of the dimensions are singular)
S(1)=size(U,1);
S(2)=size(U,2);
S(3)=size(U,3);

%outlier searching
for m=1:pass
    
    q = (t(m,:)-1)/2;
    
    for i=1:S(1)
        for j=1:S(2)
            for k=1:S(3)
                if Eval(i,j,k)==0
                    %get evaluation block with at least 8 valid points
                    s=0;
                    % This initializes the extraction domain of U
                    U_Domain=zeros(1,6);
                    while s==0
                        Imin = max([i-q(2) 1   ]);
                        Imax = min([i+q(2) S(1)]);
                        Jmin = max([j-q(1) 1   ]);
                        Jmax = min([j+q(1) S(2)]);
                        Kmin = max([k-q(3) 1   ]);
                        Kmax = min([k+q(3) S(3)]);
                        
                        % This is a vector of the current domain
                        U_Domain_Temp=[Imin,Imax,Jmin,Jmax,Kmin,Kmax];
                        
                        % This is the difference between the domain of the
                        % current iteration and the last iteration
                        r_domain=sum((U_Domain-U_Domain_Temp).^2);
                        
                        % If the domain is equal to the domain from the
                        % last iteration, then this exits the loop and uses
                        % the current domain for the UOD calculation
                        if r_domain<eps(1);
                            % This extracts the same domains from V and W
                            Vblock = V(Iind,Jind,Kind);
                            Wblock = W(Iind,Jind,Kind);
                            % This sets s equal to 1 to exit the loop
                            s=1;
                        end;
                        
                        % This sets the temporary current domain vector to
                        % the non-temporary current domain vector
                        U_Domain=U_Domain_Temp;
                        
                        % These are the indices to extract from U
                        Iind = Imin:Imax;
                        Jind = Jmin:Jmax;
                        Kind = Kmin:Kmax;
                        % This extracts the current domain
                        Ublock = U(Iind,Jind,Kind);
                        
                        if length(Ublock(~isnan(Ublock(:))))>=26
                            %                         Xblock = X(Iind,Jind)-X(i,j);
                            %                         Yblock = Y(Iind,Jind)-Y(i,j);
                            Vblock = V(Iind,Jind,Kind);
                            Wblock = W(Iind,Jind,Kind);
                            s=1;
                        else
                            q=q+1;
                        end
                    end
                    
                    %                 %distance from vector location
                    %                 Dblock = (Xblock.^2+Yblock.^2).^0.5;
                    %                 Dblock(isnan(Ublock))=nan;
                    
                    %universal outlier detection
                    Ipos = find(Iind==i);
                    Jpos = find(Jind==j);
                    Kpos = find(Kind==k);
                    [Ru]=UOD_sub(Ublock,Ipos,Jpos,Kpos);
                    [Rv]=UOD_sub(Vblock,Ipos,Jpos,Kpos);
                    [Rw]=UOD_sub(Wblock,Ipos,Jpos,Kpos);
                    
                    if Ru > tol(m) || Rv > tol(m) || Rw > tol(m)
                        %UOD threshold condition
                        U(i,j,k)=nan;
                        V(i,j,k)=nan;
                        W(i,j,k)=nan;
                        Eval(i,j,k)=m;
                    end
                    
                end
                
            end
        end
    end
end



function [R]=UOD_sub(W,p,q,r)
% --- Universal Outlier Detection Algorithm ---

%minimum variance assumption
e=0.1;

%remove value from query point
x=W(p,q,r);
W(p,q,r)=nan;

%remove any erroneous points
P=W(:);
Ps = sort(P);
Psfull = Ps(~isnan(Ps));
N=length(Psfull);

if N<=floor(length(W)/3)
    %return negative threshold value if no valid vectors in block
    R = inf;
else
    %return the median deviation normalized to the MAD
    if mod(N,2)==0
        M = (Psfull(N/2)+Psfull(N/2+1))/2;
        MADfull = sort(abs(Psfull-M));
        Q = (MADfull(N/2)+MADfull(N/2+1))/2;
        R = abs(x-M)/(Q+e);
    else
        M = Psfull((N+1)/2);
        MADfull = sort(abs(Psfull-M));
        Q = MADfull((N+1)/2);
        R = abs(x-M)/(Q+e);
    end
end



function [U,V,W,Eval] = bootstrapping(x,y,z,u,v,w,per,iter,kmax,Eval)
% Bootstrapping Validation Subfunction
%
% function [U,V,W,Eval] = bootstrapping(x,y,z,u,v,w,per,iter,kmax,Eval)
%
% per  = percent removed for each interpolation (0-1)
% iter = number of interpolations per frame (for histogram)
% kmax = number of passes

n = size(x);

M = zeros(n(1),n(2),n(3),iter);

tol = 0.3;
ktol = 1;

while tol > 0 && ktol <= kmax+1
    U = zeros(n(1),n(2),n(3),iter);
    V = zeros(n(1),n(2),n(3),iter);
    W = zeros(n(1),n(2),n(3),iter);
    
    for i = 1:iter
        clear S m Up Vp Wp Ui Vi Wi
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Data Removal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [m]= bootstrapping_dataremove(size(u),per,sum(Eval,4));
        
        S(:,1) = x(m==1);
        S(:,2) = y(m==1);
        S(:,3) = z(m==1);
        
        Up = u(m==1);
        Vp = v(m==1);
        Wp = w(m==1);
        M(:,:,:,i) = m;
        
        %         Ui = gridfit(S(:,1),S(:,2),S(:,3),Up,x(1,:),y(:,1),z(:,1));
        %         Vi = gridfit(S(:,1),S(:,2),S(:,3),Vp,x(1,:),y(:,1),z(:,1));
        %         Wi = gridfit(S(:,1),S(:,2),S(:,3),Wp,x(1,:),y(:,1),z(:,1));
        F  = TriScatteredInterp(S(:,1),S(:,2),S(:,3),Up);
        Ui = reshape(F(x(:),y(:),z(:)),[n(1) n(2) n(3)]);
        F  = TriScatteredInterp(S(:,1),S(:,2),S(:,3),Vp);
        Vi = reshape(F(x(:),y(:),z(:)),[n(1) n(2) n(3)]);
        F  = TriScatteredInterp(S(:,1),S(:,2),S(:,3),Wp);
        Wi = reshape(F(x(:),y(:),z(:)),[n(1) n(2) n(3)]);
        
        U(:,:,:,i) = Ui;
        V(:,:,:,i) = Vi;
        W(:,:,:,i) = Wi;
        
    end
    
    PBad = 0;
    for j = 1:n(1)
        for k = 1:n(2)
            for m = 1:n(3)
                if sum(isnan(U(j,k,m,:))) == 0
                    try
                        [H.U,HX.U] = hist(permute(U(j,k,m,:),[4 1 2 3]),iter/2);
                        [H.V,HX.V] = hist(permute(V(j,k,m,:),[4 1 2 3]),iter/2);
                        [H.W,HX.W] = hist(permute(W(j,k,m,:),[4 1 2 3]),iter/2);
                        
                        modeU = HX.U(H.U==max(H.U));
                        modeV = HX.V(H.V==max(H.V));
                        modeW = HX.W(H.W==max(H.W));
                        
                        tU = abs((modeU - u(j,k,m))/modeU);
                        tV = abs((modeV - v(j,k,m))/modeV);
                        tW = abs((modeW - w(j,k,m))/modeW);
                        if tU > tol || tV > tol || tW > tol && Eval(j,k,m) ~= -1
                            u(j,k,m) = modeU(1);
                            v(j,k,m) = modeV(1);
                            w(j,k,m) = modeW(1);
                            Eval(j,k,m) = 200;
                            PBad = PBad+1;
                        end
                    catch%#ok
                        Ems=lasterror;%#ok
                        fprintf('\n\n')
                        fprintf(Ems.message)
                    end
                end
            end
        end
    end
    ktol = ktol + 1;
    tol = tol-(tol/(kmax-1))*(ktol-1);
end
U=u;
V=v;
W=w;



function [M1] = bootstrapping_dataremove(DSIZE,ENUM,MASK)
% --- Bootstrapping Data Removal ---

Nx = DSIZE(1);
Ny = DSIZE(2);
Nz = DSIZE(3);
Nt = 1;

M1   = zeros(DSIZE);
RMAT = rand(Nx,Ny,Nz,Nt);
EN   = 0;

while sum(M1(:))/(Nx*Ny*Nz) < ENUM && EN < 1
    M1 = RMAT<EN;
    M1(MASK<0) = 0;
    EN = EN + 0.005;
end
M1 = double(M1);



function params=parse_pv_pairs(params,pv_pairs)
% parse_pv_pairs: parses sets of property value pairs
% usage: params=parse_pv_pairs(default_params,pv_pairs)
%
% arguments: (input)
%  default_params - structure, with one field for every potential
%             property/value pair. Each field will contain the default
%             value for that property. If no default is supplied for a
%             given property, then that field must be empty.
%
%  pv_array - cell array of property/value pairs.
%             Case is ignored when comparing properties to the list
%             of field names. Also, any unambiguous shortening of a
%             field/property name is allowed.
%
% arguments: (output)
%  params   - parameter struct that reflects any updated property/value
%             pairs in the pv_array.

npv = length(pv_pairs);
n = npv/2;

if n~=floor(n)
    error 'Property/value pairs must come in PAIRS.'
end
if n<=0
    % just return the defaults
    return
end

if ~isstruct(params)
    error 'No structure for defaults was supplied'
end

% there was at least one pv pair. process any supplied
propnames = fieldnames(params);
lpropnames = lower(propnames);
for i=1:n
    pi = lower(pv_pairs{2*i-1});
    vi = pv_pairs{2*i};
    
    ind = strmatch(pi,lpropnames,'exact');
    if isempty(ind)
        ind = strmatch(pi,lpropnames);
        if isempty(ind)
            error(['No matching property found for: ',pv_pairs{2*i-1}])
        elseif length(ind)>1
            error(['Ambiguous property name: ',pv_pairs{2*i-1}])
        end
    end
    pi = propnames{ind};
    
    % override the corresponding default in params
    params = setfield(params,pi,vi);
    
end



function [X,Y,Z,U,V,W,Eval,C]=matrixform(x,y,z,u,v,w,eval,cc)
% --- Vector to Matrix Subfunction ---

%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
c=sort(unique(z));
N=length(x);

%initialize matrices
U=nan(length(b),length(a),length(c),size(u,2));
V=nan(length(b),length(a),length(c),size(v,2));
W=nan(length(b),length(a),length(c),size(w,2));
Eval=-1*ones(length(b),length(a),length(c),size(eval,2));

%generate grid matrix
[X,Y,Z]=meshgrid(a,b,c);

%generate variable matrices (nans where no data available)
for i=1:size(U,4)
    for n=1:N
        I=find(b==y(n));
        J=find(a==x(n));
        K=find(c==z(n));
        U(I,J,K,i) = u(n,i);
        V(I,J,K,i) = v(n,i);
        W(I,J,K,i) = w(n,i);
        Eval(I,J,K,i) = eval(n);
    end
end
if ~isempty(cc)
    C=nan(length(b),length(a),length(c),size(cc,2));
    for i=1:size(c,2)
        for n=1:N
            I= b==y(n);
            J= a==x(n);
            K= c==z(n);
            C(I,J,K,i)=cc(n,i);
        end
    end
else
    C=[];
end



function [u,v,w,eval,cc]=vectorform(x,y,z,U,V,W,Eval,C)
% --- Matrix to Vector Subfunction ---
x=x(:);y=y(:);z=z(:);
%find unique x and y grid points
a=sort(unique(x));
b=sort(unique(y));
c=sort(unique(z));
N=length(x(:));

%initialize vectors
S=size(x(:));
u    = zeros(S);
v    = zeros(S);
w    = zeros(S);
eval = zeros(S);
if ~isempty(C)
    cc = zeros(S);
else
    cc = [];
end

%generate data vectors where data is available
for n=1:N
    I=find(b==y(n));
    J=find(a==x(n));
    K=find(c==z(n));
    u(n)    = U(I,J,K);
    v(n)    = V(I,J,K);
    w(n)    = W(I,J,K);
    eval(n) = Eval(I,J,K);
    if ~isempty(C)
        cc(n)    = C(I,J,K);
    end
end



function []=write_dat_val_C(fname,X,Y,Z,U,V,W,Eval,C,strand,T,frametitle,t_opt)
% --- .dat Writer Subfunction ---

if nargin<13
    t_opt=[];
end

%find I,J for plt

% This size in each dimension is accessed individually in case any of the
% dimensions are singular
S(1)=size(U,1);
S(2)=size(U,2);
S(3)=size(U,3);

%generate text file
fid = fopen(fname,'w');
if fid==-1
    error(['error creating file ',fname])
end

varlist='"X" "Y" "Z" "U" "V" "W" "Eval"';

if ~isempty(C)
    varlist=[varlist,' "C"'];
    if size(U,4)>1
        for i=2:size(U,4)
            varlist=[varlist,' "U',num2str(i-1),'" "V',num2str(i-1),'" "W',num2str(i-1),'"'];
        end
    end
    if size(C,4)>1
        for i=2:size(C,4)
            varlist=[varlist,' "C',num2str(i-1),'"'];
        end
    end
end
if ~isempty(t_opt)
    varlist=[varlist,' "t_opt"'];
end

%header lines
fprintf(fid,['TITLE        = "' frametitle '"\n']);
fprintf(fid,['VARIABLES    = ',varlist,'\n']);
fprintf(fid,'ZONE T="Time=%0.6f" I=%i J=%i K=%i C=BLACK STRANDID=%i SOLUTIONTIME = %0.6f\n',T,S(2),S(1),S(3),strand,T);


%write data
for k = 1:S(3)
    for i=1:S(1)
        for j=1:S(2)
            if isnan(U(i,j,k)) || isnan(V(i,j,k)) || isnan(W(i,j,k))
                %second check to ensure no nans present
                fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j,k),Y(i,j,k),Z(i,j,k),0,0,0,-1);
            else
                %valid data points
                fprintf(fid,'%14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e %14.6e',X(i,j,k),Y(i,j,k),Z(i,j,k),U(i,j,k),V(i,j,k),W(i,j,k),Eval(i,j,k));
            end
            
            if ~isempty(C)
                if isnan(C(i,j,k,1))
                    fprintf(fid,' %14.6e',0);
                else
                    fprintf(fid,' %14.6e',C(i,j,k,1));
                end
                if size(U,4)>1
                    for m=2:size(U,4)
                        if isnan(U(i,j,k,m)) || isnan(V(i,j,k,m)) || isnan(W(i,j,k))
                            fprintf(fid,' %14.6e %14.6e %14.6e',0,0,0);
                        else
                            fprintf(fid,' %14.6e %14.6e %14.6e',U(i,j,k,m),V(i,j,k,m),W(i,j,k,m));
                        end
                    end
                end
                if size(C,4)>1
                    for m=2:size(C,4)
                        if isnan(C(i,j,k,m))
                            fprintf(fid,' %14.6e',0);
                        else
                            fprintf(fid,' %14.6e',C(i,j,k,m));
                        end
                    end
                end
            end
            
            if ~isempty(t_opt)
                keyboard
                if isnan(t_opt(i,j))
                    fprintf(fid,' %14.6e',0);
                else
                    fprintf(fid,' %14.6e',t_opt(i,j,k));
                end
            end
            
            fprintf(fid,'\n');
        end
    end
end

fclose(fid);



function [p]=findwidth(r)
% --- Window Size Interpolation Function ---

R = [0.0000 0.0051 0.0052 0.0053 0.0055 0.0056 0.0057 0.0059 0.0060 ...
    0.0063 0.0064 0.0066 0.0067 0.0069 0.0070 0.0072 0.0074 0.0076 ...
    0.0079 0.0081 0.0083 0.0085 0.0087 0.0089 0.0091 0.0093 0.0095 ...
    0.0100 0.0102 0.0104 0.0107 0.0109 0.0112 0.0114 0.0117 0.0120 ...
    0.0125 0.0128 0.0131 0.0134 0.0137 0.0141 0.0144 0.0147 0.0151 ...
    0.0158 0.0161 0.0165 0.0169 0.0173 0.0177 0.0181 0.0185 0.0190 ...
    0.0199 0.0203 0.0208 0.0213 0.0218 0.0223 0.0228 0.0233 0.0239 ...
    0.0250 0.0256 0.0262 0.0268 0.0274 0.0281 0.0287 0.0294 0.0301 ...
    0.0315 0.0322 0.0330 0.0337 0.0345 0.0353 0.0361 0.0370 0.0378 ...
    0.0396 0.0406 0.0415 0.0425 0.0435 0.0445 0.0455 0.0466 0.0476 ...
    0.0499 0.0511 0.0522 0.0535 0.0547 0.0560 0.0573 0.0586 0.0600 ...
    0.0628 0.0643 0.0658 0.0673 0.0689 0.0705 0.0721 0.0738 0.0755 ...
    0.0791 0.0809 0.0828 0.0847 0.0867 0.0887 0.0908 0.0929 0.0951 ...
    0.0996 0.1019 0.1042 0.1067 0.1092 0.1117 0.1143 0.1170 0.1197 ...
    0.1253 0.1283 0.1312 0.1343 0.1374 0.1406 0.1439 0.1473 0.1507 ...
    0.1578 0.1615 0.1652 0.1691 0.1730 0.1770 0.1812 0.1854 0.1897 ...
    0.1986 0.2033 0.2080 0.2128 0.2178 0.2229 0.2281 0.2334 0.2388 ...
    0.2501 0.2559 0.2619 0.2680 0.2742 0.2806 0.2871 0.2938 0.3006 ...
    0.3148 0.3221 0.3296 0.3373 0.3451 0.3531 0.3613 0.3696 0.3781 ...
    0.3957 0.4048 0.4140 0.4233 0.4329 0.4425 0.4524 0.4623 0.4724 ...
    0.4930 0.5034 0.5139 0.5244 0.5351 0.5457 0.5564 0.5672 0.5779 ...
    0.5992 0.6099 0.6204 0.6309 0.6414 0.6517 0.6619 0.6720 0.6819 ...
    0.7014 0.7109 0.7203 0.7295 0.7385 0.7473 0.7559 0.7643 0.7726 ...
    0.7884 0.7960 0.8035 0.8107 0.8177 0.8245 0.8311 0.8376 0.8438 ...
    0.8556 0.8613 0.8667 0.8720 0.8771 0.8820 0.8867 0.8913 0.8957 ...
    0.9041 0.9080 0.9118 0.9155 0.9190 0.9224 0.9256 0.9288 0.9318 ...
    0.9374 0.9401 0.9426 0.9451 0.9474 0.9497 0.9519 0.9539 0.9559 ...
    0.9597 0.9614 0.9631 0.9647 0.9662 0.9677 0.9691 0.9705 0.9718 ...
    0.9742 0.9753 0.9764 0.9775 0.9785 0.9794 0.9803 0.9812 0.9820 ...
    0.9836 0.9843 0.9850 0.9857 0.9863 0.9869 0.9875 0.9881 0.9886 ...
    0.9896 0.9900 0.9905 0.9909 0.9913 0.9917 0.9921 0.9924 0.9928 ...
    0.9934 0.9937 0.9940 0.9943 0.9945 0.9948 0.9950 1.0000]';

P = [500.0000 245.4709 239.8833 234.4229 229.0868 223.8721 218.7762 213.7962 208.9296 ...
    199.5262 194.9845 190.5461 186.2087 181.9701 177.8279 173.7801 169.8244 165.9587 ...
    158.4893 154.8817 151.3561 147.9108 144.5440 141.2538 138.0384 134.8963 131.8257 ...
    125.8925 123.0269 120.2264 117.4898 114.8154 112.2018 109.6478 107.1519 104.7129 ...
    100.0000  97.7237  95.4993  93.3254  91.2011  89.1251  87.0964  85.1138  83.1764 ...
    79.4328  77.6247  75.8578  74.1310  72.4436  70.7946  69.1831  67.6083  66.0693 ...
    63.0957  61.6595  60.2560  58.8844  57.5440  56.2341  54.9541  53.7032  52.4807 ...
    50.1187  48.9779  47.8630  46.7735  45.7088  44.6684  43.6516  42.6580  41.6869 ...
    39.8107  38.9045  38.0189  37.1535  36.3078  35.4813  34.6737  33.8844  33.1131 ...
    31.6228  30.9030  30.1995  29.5121  28.8403  28.1838  27.5423  26.9153  26.3027 ...
    25.1189  24.5471  23.9883  23.4423  22.9087  22.3872  21.8776  21.3796  20.8930 ...
    19.9526  19.4984  19.0546  18.6209  18.1970  17.7828  17.3780  16.9824  16.5959 ...
    15.8489  15.4882  15.1356  14.7911  14.4544  14.1254  13.8038  13.4896  13.1826 ...
    12.5893  12.3027  12.0226  11.7490  11.4815  11.2202  10.9648  10.7152  10.4713 ...
    10.0000   9.7724   9.5499   9.3325   9.1201   8.9125   8.7096   8.5114   8.3176 ...
    7.9433   7.7625   7.5858   7.4131   7.2444   7.0795   6.9183   6.7608   6.6069 ...
    6.3096   6.1660   6.0256   5.8884   5.7544   5.6234   5.4954   5.3703   5.2481 ...
    5.0119   4.8978   4.7863   4.6774   4.5709   4.4668   4.3652   4.2658   4.1687 ...
    3.9811   3.8905   3.8019   3.7154   3.6308   3.5481   3.4674   3.3884   3.3113 ...
    3.1623   3.0903   3.0200   2.9512   2.8840   2.8184   2.7542   2.6915   2.6303 ...
    2.5119   2.4547   2.3988   2.3442   2.2909   2.2387   2.1878   2.1380   2.0893 ...
    1.9953   1.9498   1.9055   1.8621   1.8197   1.7783   1.7378   1.6982   1.6596 ...
    1.5849   1.5488   1.5136   1.4791   1.4454   1.4125   1.3804   1.3490   1.3183 ...
    1.2589   1.2303   1.2023   1.1749   1.1482   1.1220   1.0965   1.0715   1.0471 ...
    1.0000   0.9772   0.9550   0.9333   0.9120   0.8913   0.8710   0.8511   0.8318 ...
    0.7943   0.7762   0.7586   0.7413   0.7244   0.7079   0.6918   0.6761   0.6607 ...
    0.6310   0.6166   0.6026   0.5888   0.5754   0.5623   0.5495   0.5370   0.5248 ...
    0.5012   0.4898   0.4786   0.4677   0.4571   0.4467   0.4365   0.4266   0.4169 ...
    0.3981   0.3890   0.3802   0.3715   0.3631   0.3548   0.3467   0.3388   0.3311 ...
    0.3162   0.3090   0.3020   0.2951   0.2884   0.2818   0.2754   0.2692   0.2630 ...
    0.2512   0.2455   0.2399   0.2344   0.2291   0.2239   0.2188   0.2138   0.2089 ...
    0.1995   0.1950   0.1905   0.1862   0.1820   0.1778   0.1738   0.0000]';

p=interp1q(R,P,r);



function [xh]=wlsq(y,H,W)
% --- Weighted Least Squares Fit for Phase Correlation ---
tempmat=sortrows([y',H',W'],2);
y=tempmat(:,1);
H=tempmat(:,2);
W=diag(tempmat(:,3));

% xh=inv(H'*W*H)*H'*W*y;
xh=(H'*W*H)\(H'*W*y);
