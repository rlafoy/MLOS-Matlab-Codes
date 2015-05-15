function pranaPIV3D_parallel_code_runner_03(piv_processing_parameters);
% This function is calls the prana3D code, but is designed to run in
% parallel on multiple computers.
%
% This version is similar to version 01 (which doesn't have a version
% number), except this version checks that all the frames are processed
% before exciting the loop in case any of the individual matlab workers are
% closed or the matlab workers crash or the computer crashes or the power
% goes out or the zombie apocolypse occurs or . . . well you get the point.
%
% Version 03 is designed to run using an input structure
% 'piv_processing_parameters' to control how the processing is ran.  The
% parameters defined in this structure will overwrite any parameters
% defined in this code, although if a parameter is not defined in
% piv_processing_parameters, the parameter defined in this code will be
% used for the processing.
%
% Authors: Rod La Foy, Sam Raben
% Version 03 Created On: 15 March 2013
% Last Edited: 21 June 2013

% This is the top directory to perform the processing in
top_write_directory=piv_processing_parameters.Data.outdirec;

% This is the filename of the data file for communicating between multiple
% computers that contains a list of frames that have begun processing.
% The frames listed in this file are not necessarily completed.  All
% parallel computers must be able to read/write this file.
processing_com_filename=[top_write_directory,'processing_com_file.dat'];

% This is the filename of the data file for communicating between multiple
% computers that contains a list of frames that have completed processing.
% The frames listed in this file have been fully processed.  All parallel 
% computers must be able to read/write this file.
completed_com_filename=[top_write_directory,'completed_com_file.dat'];

% This checks whether the communication files exist and if not creates them
if not(exist(processing_com_filename,'file'));
    % This opens the file for writing
    fid=fopen(processing_com_filename,'w');
    % This closes the file
    fclose(fid);
end;
if not(exist(completed_com_filename,'file'));
    % This opens the file for writing
    fid=fopen(completed_com_filename,'w');
    % This closes the file
    fclose(fid);
end;

% This is the range of frame numbers to process (this corresponds to the
% values that would ordinarily be assigned to
%  'Data.imfstep'
%  'Data.imfstart'
%  'Data.imfend'
% in the 'Data Construction' portion of the processing file.  These are
% saved as doubles to generate the processing vector and are later converted
% to strings for input to the 'Data' structure.
frame_pros_step=piv_processing_parameters.frame_correlation_step;
frame_pros_start=piv_processing_parameters.frame_processing_start;
frame_pros_end=piv_processing_parameters.frame_processing_end;

% This is a vector of frames to process
frame_pros_vect=frame_pros_start:frame_pros_step:frame_pros_end;

% This creates the data structure for processing
Data=generate_data_structure;

% This overwrites the 'Data' structure with the data defined in the
% 'piv_processing_parameters' structure
Data=overwrite_data_structure(Data,piv_processing_parameters);

tic;

% This iterates through the frame processing vector processing the frames
% in parallel on the several computers.
while true;

    % This opens the completed file for writing
    fid=fopen(completed_com_filename,'a+');
    % This reads the file data
    completed_com_data=fread(fid,inf,'uint16');
    % This checks whether all frames have been processed and if so exits
    % the loop, but first checks whether any frames have been written.
    if not(isempty(completed_com_data));
        % This sorts the completed data vector for comparison to the vector of
        % frames that should be processed
        completed_com_data=sort(completed_com_data,1,'ascend');
        % This removes any non-unique values in the vector. (These shouldn't
        % really exist, but there is a very small probability - less than 1
        % in 10 million - that two computers will simultaneously load the
        % vector and then write that they are processing the same vector back
        % to the file; in this case the second computer will likely overwrite
        % the data file of the first so the output data shouldn't be effected.)
        completed_com_data=unique(completed_com_data);
        % This checks if the length of this vector equals the length of the
        % frame processing vector
        if length(completed_com_data)==length(frame_pros_vect);
            % This is the sum of the differences between the vectors
            frame_vect_residual=sum((double(completed_com_data)-frame_pros_vect).^2);
            % This checks whether the vectors are equal (they should be if
            % they have the same length - but weird stuff happens)
            if frame_vect_residual<1e-14;
                % This closes the communication file
                fclose(fid);
                % This ends the processing loop
                break;
            else;
                % This closes the communication file
                fclose(fid);
                % This displays an error saying something weird is going
                % on in the processing
                error(['The processing file frame list has the same ',...
                    'length as the processing frame vector, but they are ',...
                    'not equal. This likely means that a bug exists ',...
                    'somewhere in the processing code or that there was a ',...
                    'write error.']);
            end;
        else;
            % This closes the communications file
            fclose(fid);
        end;
    else;
        % This closes the communications file
        fclose(fid);
    end;

    % Since the loop has gotten this far, it means that not all frames have
    % been processed.  In this case, this loads the list of frames being
    % processed and based upon this list calculates the next frame number
    % needing to be processed.
    
    % This opens the processing file for writing
    fid=fopen(processing_com_filename,'a+');
    % This reads the file data
    processing_com_data=fread(fid,inf,'uint16');
    % This checks which frames have begun processing, but first checks
    % whether any have begun processing
    if not(isempty(processing_com_data));
        % This is the maximum frame that has begun processing
        max_processing_frame=max(processing_com_data);
        % If the maximum frame does not equal the last frame then the
        % current frame is incremented, otherwise it is set as the last
        % frame
        if max_processing_frame<frame_pros_end;
            % This increments the current frame to process
            current_frame=max_processing_frame+frame_pros_step;
        else;
            % This sets the current processing frame to the last frame.
            % (This will result in at least one frame being redundantly
            % processed, which should be fixed, but there isn't an easy,
            % obvious way around this problem while still ensuring all
            % frames are processed.)
            current_frame=frame_pros_end;
        end;
    else;
        % Otherwise this sets the current frame to process to the first
        % frame to be processed
        current_frame=frame_pros_start;
    end;
    % This writes the current processing frame to the file
    fwrite(fid,current_frame,'uint16');
    % This closes the processing file
    fclose(fid);
    
    % This rewrites the values of 'Data.imfstart' and 'Data.imfend' for
    % calculating the current frame
    Data.imfstart=num2str(current_frame);
    Data.imfend=num2str(current_frame);
    
    % This processes the current frame
    pranaPIV3Dcode_mem_opt_04(Data);
    
    % This opens the completed com file for writing that the current frame
    % has finished processing
    fid=fopen(completed_com_filename,'a+');
    % This writes the currently completed frame to the file
    fwrite(fid,current_frame,'uint16');
    % This closes the completed com file
    fclose(fid);
    
end;

disp(['Prana3D Processing time: ',num2str(toc),' s.']);


function Data=overwrite_data_structure(Data,piv_processing_parameters);
% This uses the input parameter structure 'piv_processing_parameters' to 
% overwrite the fields defined within this code for the structure 'Data'.

% These are the fields within the Data structure of the input parameters
data_fieldnames=fieldnames(piv_processing_parameters.Data);
% This iterates through the fieldnames changing the values of the 'Data'
% structure
for n=1:length(data_fieldnames);
    % This is the current data fieldname
    current_fieldname=data_fieldnames{n};
    % This checks whether the current field name is a 'PIV' structure
    if isstruct(eval(['piv_processing_parameters.Data.',current_fieldname]))&&strcmp('PIV',current_fieldname(1:3));
        % These are the PIV fieldnames
        piv_fieldnames=fieldnames(eval(['piv_processing_parameters.Data.',current_fieldname]));
        % This iterates through the PIV fieldnames overwriting the values
        % already stored in the 'Data' structure
        for m=1:length(piv_fieldnames);
            % This is the current 'PIV' fieldname
            piv_current_fieldname=piv_fieldnames{m};
            % This is the string into the structure
            % 'piv_processing_parameters' to read the data from
            structure_read_string=['piv_processing_parameters.Data.',current_fieldname,'.',piv_current_fieldname];
            % This is the string into the structure 'Data' to write the
            % data to
            structure_write_string=['Data.',current_fieldname,'.',piv_current_fieldname];
            % This overwrites the value stored in the 'Data' structure
            eval([structure_write_string,'=',structure_read_string,';']);
        end;
    else;
        % This is the string into the structure
        % 'piv_processing_parameters' to read the data from
        structure_read_string=['piv_processing_parameters.Data.',current_fieldname];
        % This is the string into the structure 'Data' to write the
        % data to
        structure_write_string=['Data.',current_fieldname];
        % This overwrites the value stored in the 'Data' structure
        eval([structure_write_string,'=',structure_read_string,';']);
    end;
end;


function Data=generate_data_structure;
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Data Construction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.version                    = '4.3';
% Data.imbase                     = 'S';
Data.imzeros                    = '2';
Data.imext                      = 'mat';
Data.imcstep                    = '1';
Data.imfstep                    = '1';
Data.imfstart                   = '1';
Data.imfend                     = '1';
Data.wrmag                      = '1';
Data.wrsamp                     = '1';
Data.wrsep                      = '1';
Data.batchname                  = 'Test3D';
Data.datout                     = '1';
Data.multiplematout             = '1';
Data.exp_date                   = '';
Data.exp_L                      = '';
Data.exp_v0                     = '';
Data.exp_notes                  = {'Camera Description:'  ''  'Lens Description:'  ''  'Notes:'  ''};
Data.exp_density                = '1000';
Data.exp_viscosity              = '1.308e-3';
Data.exp_surfacetension         = '0.07197';
Data.exp_partD                  = '';
Data.exp_partdensity            = '';
Data.exp_wavelength             = '.532';
Data.exp_pixelsize              = '';
Data.exp_lensfocal              = '';
Data.exp_lensfnum               = '';
Data.exp_micro                  = '0';
Data.exp_NA                     = '';
Data.exp_n                      = '';
Data.exp_Re                     = '';
Data.exp_St                     = '';
Data.exp_M                      = '';
Data.exp_ROI                    = '';
Data.exp_diffractiondiameter    = '';
Data.exp_depthoffocus           = '';
Data.masktype                   = 'none';
Data.staticmaskname             = '';
Data.maskbase                   = 'maskfor_Img_';
Data.maskzeros                  = '4';
Data.maskext                    = 'tif';
Data.maskfstep                  = '1';
Data.maskfstart                 = '1';
Data.passes                     = '2';
Data.method                     = '2';
Data.velinterp                  = '3';
Data.iminterp                   = '1';
Data.framestep                  = '3';
Data.PIVerror                   = '0.1';
Data.splash                     = '1';
Data.par                        = '0';
Data.parprocessors              = '8';

Data.imdirec                    = '';
Data.maskdirec                  = '';
Data.outdirec                   = '';

Data.cpass                      = '1';

Data.imbase                     = 'frame_';
                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIV0
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.PIV0.winres                    = '32,32';
Data.PIV0.winsize                   = '64,64';
Data.PIV0.winauto                   = '1';
Data.PIV0.gridres                   = '8,8';
Data.PIV0.winoverlap                = '75,75';
Data.PIV0.gridtype                  = '1';
Data.PIV0.gridbuf                   = '8,8';
Data.PIV0.BWO                       = '0,0';
Data.PIV0.corr                      = '2';
Data.PIV0.RPCd                      = '2.8';
Data.PIV0.zeromean                  = '0';
Data.PIV0.peaklocator               = '1';
Data.PIV0.velsmooth                 = '0';
Data.PIV0.velsmoothfilt             = '2';
Data.PIV0.val                       = '0';
Data.PIV0.uod                       = '1';
Data.PIV0.bootstrap                 = '0';
Data.PIV0.thresh                    = '0';
Data.PIV0.uod_type                  = '2';
Data.PIV0.uod_window                = '3,3;3,3';
Data.PIV0.uod_thresh                = '3,2';
Data.PIV0.bootstrap_percentsampled  = '15';
Data.PIV0.bootstrap_iterations      = '700';
Data.PIV0.bootstrap_passes          = '12';
Data.PIV0.valuthresh                = '-16,16';
Data.PIV0.valvthresh                = '-16,16';
Data.PIV0.valwthresh                = '-16,16';
Data.PIV0.valextrapeaks             = '0';
Data.PIV0.savepeakinfo              = '0';
Data.PIV0.corrpeaknum               = '1';
Data.PIV0.savepeakmag               = '0';
Data.PIV0.savepeakvel               = '0';
Data.PIV0.outbase                   = 'PIV_';
Data.PIV0.write                     = '1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIV1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.PIV1.winres                    = '32,32,32'; % Effective resolution
Data.PIV1.winsize                   = '64,64,64'; % Actual resoltuion
Data.PIV1.winauto                   = '1';
Data.PIV1.gridres                   = '16,16,16';%'32,32,32'; % Distance between windows in voxels
Data.PIV1.winoverlap                = '25,25,25';
Data.PIV1.gridtype                  = '1';
Data.PIV1.gridbuf                   = '8,8,8';
Data.PIV1.BWO                       = '0,0,0';
Data.PIV1.corr                      = '2';  % corr = 1 gives SCC, corr = 2 gives RPC
Data.PIV1.RPCd                      = '2.8';
Data.PIV1.zeromean                  = '1';
Data.PIV1.peaklocator               = '1';
Data.PIV1.velsmooth                 = '1';
Data.PIV1.velsmoothfilt             = '2';
Data.PIV1.val                       = '1';
Data.PIV1.uod                       = '1';
Data.PIV1.bootstrap                 = '0';
Data.PIV1.thresh                    = '1';
Data.PIV1.uod_type                  = '2';
Data.PIV1.uod_window                = '3,3,3;3,3,3';
Data.PIV1.uod_thresh                = '3,2';
Data.PIV1.bootstrap_percentsampled  = '15';
Data.PIV1.bootstrap_iterations      = '7';
Data.PIV1.bootstrap_passes          = '12';
Data.PIV1.valuthresh                = '-16,16';
Data.PIV1.valvthresh                = '-16,16';
Data.PIV1.valwthresh                = '-16,16';
Data.PIV1.valextrapeaks             = '0';
Data.PIV1.savepeakinfo              = '0';
Data.PIV1.corrpeaknum               = '1';
Data.PIV1.savepeakmag               = '0';
Data.PIV1.savepeakvel               = '0';
Data.PIV1.outbase                   = 'PIV_Pass_1_';
Data.PIV1.write                     = '1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIV2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Data.PIV2.winres                    = '16,16,16';
Data.PIV2.winsize                   = '32,32,32';
Data.PIV2.winauto                   = '1';
Data.PIV2.gridres                   = '8,8,8';
Data.PIV2.winoverlap                = '25,25,25';
Data.PIV2.gridtype                  = '1';
Data.PIV2.gridbuf                   = '8,8,8';
Data.PIV2.BWO                       = '0,0,0';
Data.PIV2.corr                      = '2';
Data.PIV2.RPCd                      = '2.8';
Data.PIV2.zeromean                  = '0';
Data.PIV2.peaklocator               = '1';
Data.PIV2.velsmooth                 = '0';
Data.PIV2.velsmoothfilt             = '2';
Data.PIV2.val                       = '1';
Data.PIV2.uod                       = '1';
Data.PIV2.bootstrap                 = '0';
Data.PIV2.thresh                    = '0';
Data.PIV2.uod_type                  = '2';
Data.PIV2.uod_window                = '3,3,3;3,3,3';
Data.PIV2.uod_thresh                = '3,2';
Data.PIV2.bootstrap_percentsampled  = '15';
Data.PIV2.bootstrap_iterations      = '700';
Data.PIV2.bootstrap_passes          = '12';
Data.PIV2.valuthresh                = '-16,16';
Data.PIV2.valvthresh                = '-16,16';
Data.PIV2.valwthresh                = '-16,16';
Data.PIV2.valextrapeaks             = '0';
Data.PIV2.savepeakinfo              = '1';
Data.PIV2.corrpeaknum               = '2';
Data.PIV2.savepeakmag               = '1';
Data.PIV2.savepeakvel               = '0';
Data.PIV2.outbase                   = 'PIV_Pass_2_';
Data.PIV2.write                     = '1';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIV3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Data.PIV3.winres                    = '32,32,32';
% Data.PIV3.winsize                   = '64,64,64';
% Data.PIV3.winauto                   = '1';
% Data.PIV3.gridres                   = '8,8,8';
% Data.PIV3.winoverlap                = '75,75,75';
% Data.PIV3.gridtype                  = '1';
% Data.PIV3.gridbuf                   = '8,8,8';
% Data.PIV3.BWO                       = '0,0,0';
% Data.PIV3.corr                      = '2';
% Data.PIV3.RPCd                      = '2.8';
% Data.PIV3.zeromean                  = '0';
% Data.PIV3.peaklocator               = '1';
% Data.PIV3.velsmooth                 = '0';
% Data.PIV3.velsmoothfilt             = '2';
% Data.PIV3.val                       = '0';
% Data.PIV3.uod                       = '1';
% Data.PIV3.bootstrap                 = '0';
% Data.PIV3.thresh                    = '0';
% Data.PIV3.uod_type                  = '2';
% Data.PIV3.uod_window                = '3,3,3;3,3,3';
% Data.PIV3.uod_thresh                = '3,2';
% Data.PIV3.bootstrap_percentsampled  = '15';
% Data.PIV3.bootstrap_iterations      = '700';
% Data.PIV3.bootstrap_passes          = '12';
% Data.PIV3.valuthresh                = '-16,16';
% Data.PIV3.valvthresh                = '-16,16';
% Data.PIV3.valwthresh                = '-16,16';
% Data.PIV3.valextrapeaks             = '0';
% Data.PIV3.savepeakinfo              = '0';
% Data.PIV3.corrpeaknum               = '1';
% Data.PIV3.savepeakmag               = '0';
% Data.PIV3.savepeakvel               = '0';
% Data.PIV3.outbase                   = 'PIV_Pass3_';
% Data.PIV3.write                     = '1';

