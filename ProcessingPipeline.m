function ProcessingPipeline
% EEG processing pipeline  (Lauren Ostrowski, 06-21-2019
%    email lauren_ostrowski@brown.edu with any concerns)
% 
% Load EEG file(s) from a variety of collection systems and file types
%   as well as native EEGLab .set data structure
% Bandpass data at the range you specify below (default: 0.1 - 70 Hz)
% Notch filter data (default: 60 Hz)
% Re-reference the data to the average reference (optional)
% Interpolate bad channels identified by the FASTER algorithm, with manual
%   verification (bad channels in red, external channels in black)
% Run ICA and reject bad components identified by the FASTER algorithm, 
%   with manual verification (bad components in red)
% Segment data into task-specific runs (eyes closed, eyes open, etc.)
% Epoch data and reject very bad epochs (~70% sensitivity)
% Interpolate bad channels within epochs (identified as above)
% 
% Also save topological plots of the independent components (in the 
%  'Intermediate' directory) for your records.
% NOTE: the REJECTED components will have a RED BACKGROUND
%
%% VARIABLES
% Channel options
o.channel_options.eeg_chans=1:128; % Enter the channel indices of all EEG channels
o.channel_options.eog_chans=[8 25 126 127]; % If EOG collected, enter the EOG channel indices (expected that they are within the EEG channel set [not external])
o.channel_options.ext_chans=[]; % If external channel data collected (ex. EKG), enter channel index(ices)

o.channel_options.initialRef=17; % Enter the channel(s) to which input data will be initially referenced (suggest the nasion)
o.channel_options.do_avg_reref=0; % Set equal to 1 to re-reference to average
o.channel_options.mastoidRef=0; % Set equal to 1 to use mastoid reference
o.channel_options.mastoidChannels = [56,107];

% Epoch options
o.epoch_options.epoch_on=1; % Set equal to 1 to perform epoching
o.epoch_options.epoch_length=1; % Set length of epoch (in seconds)

% Filtering options - set equal to 1 to perform filtering
o.filter_options.hpf_on=1;
o.filter_options.hpf_freq=0.1; % Highpass frequency

o.filter_options.lpf_on=1;
o.filter_options.lpf_freq=70; % Lowpass frequency

o.filter_options.notch_on=1;
o.filter_options.notch_freq=60; % Notch frequency

% Data segmentation options - set equal to 1 to segment by event flags
o.epoch_options.segmentByEvent=1;
% Name event flags of interest to be saved in outputs
o.epoch_options.eventNames={'EyesClosed','EyesOpen',...
    'RHEyesOpen','LHEyesOpen','RHEyesClosed','LHEyesClosed'};
% Event flags of interest as stored in EEG data structure
o.epoch_options.eventFlags={{'EyCl','eyec','eycl'},{'EyOp','eyeo','eyop'},...
    {'RHEO','rheo'},{'LHEO','lheo'},{'RHEC','rhec'},{'LHEC','lhec'}};

%% Set up necessary files to track processing/catch errors
fprintf('**************************************\n');
fprintf('* Running EEG Preprocessing Pipeline *\n');
fprintf('**************************************\n');
fprintf('Make sure you have correctly set all customizable variables:\n');
fprintf([' - EEG channels:' sprintf(' %d', o.channel_options.eeg_chans) '\n']);
if ~isempty(o.channel_options.eog_chans)
    fprintf([' - EOG channels:' sprintf(' %d', o.channel_options.eog_chans) '\n']);
else
    fprintf(' - EOG channel data not collected.\n');
end
if ~isempty(o.channel_options.ext_chans)
    fprintf([' - External channel(s):' sprintf(' %d', o.channel_options.ext_chans) '\n']);
else
    fprintf(' - External channel data not collected.\n');
end
fprintf([' - Reference electrode:' sprintf(' %d', o.channel_options.initialRef) '\n'])
if o.channel_options.do_avg_reref
    fprintf(' - Average re-reference will be performed\n'); 
end
if o.filter_options.hpf_on
    fprintf([' - Data will be highpass filtered at ' num2str(o.filter_options.hpf_freq) ' Hz\n']);
end
if o.filter_options.lpf_on
    fprintf([' - Data will be lowpass filtered at ' num2str(o.filter_options.lpf_freq) ' Hz\n']);
end
if o.filter_options.notch_on
    fprintf([' - Data will be notch filtered at ' num2str(o.filter_options.notch_freq) ' Hz\n']);
end
if o.epoch_options.epoch_on
    fprintf([' - Data will be split into epochs of ' num2str(o.epoch_options.epoch_length) ' second(s)\n']);
end
if o.epoch_options.segmentByEvent
    fprintf([' - Data will be segmented by event flags: ' strjoin(o.epoch_options.eventNames,', ')]);
end; fprintf(newline);
fprintf('Select start directory...\n'); startDir=uigetdir(pwd);
if startDir==0
    fprintf('  Program terminated: You must select start directory.\n');
    fprintf('Please run program again and select a start directory if you wish to proceed.\n');
    return;
end
fprintf('Start directory selected.\n\n');
fprintf('Select output directory...\n'); outDir=uigetdir(pwd);
if outDir==0
    fprintf('  Program terminated: You must select output directory.\n'); 
    fprintf('Please run program again and select an output directory if you wish to proceed.\n');
    return;
end
fprintf('Output directory selected.\n');
eeglab;close;initDir=pwd; % set up path hierarchies by calling eeglab
Qname='ProcQ.eegQ';

% Search start directory for EEG files
[pSETlist, nSETlist] = extsearchc(startDir,'.set',0); % Native eeglab data structure
[pRAWlist, nRAWlist] = extsearchc(startDir,'.raw',0); % EGI (Electrical Geodesics Incorporated) continuous file
[pEDFlist, nEDFlist] = extsearchc(startDir,'.edf',0); % joint European 16-bit data format
[pBDFlist, nBDFlist] = extsearchc(startDir,'.bdf',0); % 24-bit variant of the EDF format used by EEG systems manufactured by BioSemi
[pSMAlist, nSMAlist] = extsearchc(startDir,'.sma',0); % Snapmaster file
[pCNTlist, nCNTlist] = extsearchc(startDir,'.cnt',0); % Neuroscan continuous file
% Note that there are frequently issues importing Neuroscan data, and it
% may be necessary to input manually by calling "eeglab", loading the data
% and saving as a .SET file before running the preprocessing pipeline

plist=[pSETlist pRAWlist pEDFlist pBDFlist pSMAlist pCNTlist];
nlist=[nSETlist nRAWlist nEDFlist nBDFlist nSMAlist nCNTlist];
oplist=cell(size(plist));
x=true(size(plist));

c=clock;
% Set up output hierarchy
filepath=outDir;
if (length(plist) == 1) && isempty(outDir)
    oplist{1}=outDir;
else
    for i=1:length(plist)
        oplist{i}=[filepath filesep nlist{i}(1:end-4) '_cleaned_' ...
            num2str(c(2)) '-' num2str(c(3)) '-' num2str(c(1))];
        if exist(oplist{i},'dir')
            promptMessage = sprintf('The output directory already exists:\n%s\nDo you want to overwrite it?', ...
                [nlist{i}(1:end-4) '_cleaned_' num2str(c(2)) '-' num2str(c(3)) '-' num2str(c(1))]);
            titleBarCaption = 'Overwrite?';
            buttonText = questdlg(promptMessage, titleBarCaption, 'Yes', 'No', 'Yes');
            if strcmpi(buttonText, 'No')
                options.Resize='on';options.WindowStyle='normal';options.Interpreter='tex';
                newName = inputdlg('Enter new name:',...
                    'New output folder name',[1 50],{''},options);
                if isempty(newName)
                    newName=['untitled_' num2str(c(2)) '-' num2str(c(3)) '-' num2str(c(1))]; 
                else
                    newName=newName{1};
                    if ~strcmp(newName(max(length(newName)-9,1):end),...
                            [num2str(c(2)) '-' num2str(c(3)) '-' num2str(c(1))])
                       newName=[newName '_' num2str(c(2)) '-' num2str(c(3)) '-' num2str(c(1))];
                    end
                end
                oplist{i}=[filepath filesep newName]; mkdir(oplist{i});
            end
        else
           mkdir(oplist{i});
        end
    end
end
plist={plist{x}};
nlist={nlist{x}};
oplist={oplist{x}};
o.file_options.plist=plist;
o.file_options.nlist=nlist;
o.file_options.oplist=oplist;


% Add vars to the 'ProcQ.eegQ' file
Q.plist=plist;
Q.nlist=nlist;
Q.plist_rel=cell(0);
Q.oplist_rel=cell(0);
for v=1:length(plist)
    Q.plist_rel{v}=find_relative_path(plist{v},startDir);
    try
        Q.oplist_rel{v}=find_relative_path(oplist{v},startDir);
    catch ME
        if (strcmp(ME.identifier,'MATLAB:cd:NonExistentDirectory'))
            if startsWith(outDir,filesep); mkdir(outDir)
            else; mkdir([initDir filesep outDir]); end
            Q.oplist_rel{v}=find_relative_path(oplist{v},startDir);
        else; rethrow(ME)
        end
    end
end
Q.outDir_rel=cell(0);
if ~isempty(outDir)
    Q.outDir_rel=find_relative_path(outDir,startDir);
end
if (~exist([startDir filesep 'Processing'],'dir'))
    mkdir([startDir filesep 'Processing']);
end
my_comp_num=1;
Q.comp_nums=1;
Q.finished=0;
Q.processed=zeros(size(Q.plist));
Q.errors=zeros(size(Q.plist));
Q.next_file=1;
save([startDir filesep Qname],'Q');

% Set up error checks and queue files for processing
all_errors=cell(0);
error_indices=zeros(size(plist));
first_file=1;
had_error=0;
my_proc_file=[];
make_processing_file();
EEG_state=[];
while 1
    L=load([startDir filesep Qname],'-mat');
    Q=L.Q;
    if ~first_file
        if ~had_error
            Q.processed(current_file)=1;
        else
            Q.errors(current_file)=1;
            if exist('m','var')
                if isempty(outDir)
                    if exist([startDir 'Preprocessing_errors.mat'],'file')
                        L=load([startDir 'Preprocessing_errors.mat'],'-mat');
                        all_errors=L.all_errors;
                    end
                    all_errors{end+1,1}=m;
                    all_errors{end,2}=EEG_state;
                    save([startDir filesep 'Preprocessing_errors.mat'],'all_errors','-mat');
                    error_indices(current_file)=size(all_errors,1);
                else
                    if exist([outDir filesep 'Preprocessing_errors.mat'],'file')
                        L=load([outDir filesep 'Preprocessing_errors.mat'],'-mat');
                        all_errors=L.all_errors;
                    end
                    all_errors{end+1,1}=m;
                    all_errors{end,2}=EEG_state;
                    save([outDir filesep 'Preprocessing_errors.mat'],'all_errors','-mat');
                    error_indices(current_file)=size(all_errors,1);
                end
            end
        end
    end
    had_error=0;
    if Q.next_file>length(Q.plist)
        Q.finished=1;
    end
    if Q.finished
        Q.comp_nums=setdiff(Q.comp_nums,my_comp_num);
        delete_processing_file();
        break;
    end

    % Iterate through list of files for processing
    current_file=Q.next_file;
    Q.next_file=Q.next_file+1;
    save([startDir filesep Qname],'Q');
    
    o.file_options.current_file = [plist{current_file} filesep nlist{current_file}];
    o.file_options.current_file_num=current_file;

    searchstring2=nlist{current_file};
    if ~isempty(strfind(nlist{current_file},searchstring2))
        t0=tic;
        fprintf('******************\n');
        fprintf('* File %.3d / %.3d *\n',current_file,length(nlist));
        fprintf('******************\n');

        log_file = fopen([oplist{current_file} filesep nlist{current_file}(1:end-4) '.log'],'a+');
        try
            o.file_options.output_folder_name=outDir;
            % RUN PREPROCESSING
            preprocess(o,log_file,t0);
        catch ME
            fprintf('\nError in function %s() at line %d:\nERROR: %s.\n', ...
                ME.stack(1).name, ME.stack(1).line, ME.message);
            try fclose(log_file); catch; end
            had_error=1;
        end
    else
        fprintf('Skipped file.\n');
    end

    % After processing
    first_file=0;
end

%%%%%%%%%%%%%%%%%%%
% Post processing %
%%%%%%%%%%%%%%%%%%%
D=dir([startDir filesep 'Processing']);
if length(D)>2
    fprintf('***************************\n');
    fprintf('* Pre-processing complete *\n');
    fprintf('***************************\n');
    fprintf('Finished processing all files.\n');
    return;
end

if isempty(outDir)
    top_log = fopen([startDir filesep 'Preprocessing.log'],'a');
    if exist([startDir filesep 'Preprocessing_errors.mat'],'file')
        L=load([startDir filesep 'Preprocessing_errors.mat'],'-mat');
        all_errors=L.all_errors;
    end
else
    top_log = fopen([outDir filesep 'Preprocessing.log'],'a');
    if exist([outDir filesep 'Preprocessing_errors.mat'],'file')
        L=load([outDir filesep 'Preprocessing_errors.mat'],'-mat');
        all_errors=L.all_errors;
    end
end

c=clock;
months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
fprintf(top_log,'\n%d/%s/%d %d:%d:%d\n',c(3),months{c(2)},c(1),c(4),c(5),round(c(6)));

for v=1:length(plist)
    fprintf(top_log,'%s%s%s:\n',plist{v},filesep,nlist{v});
end

delete([startDir filesep Qname]);
if length(dir([startDir filesep 'Processing']))==2
    rmdir([startDir filesep 'Processing']);
end

fprintf('**************************\n');
fprintf('* Preprocessing Finished *\n');
fprintf('*  %.3d processed         *\n',sum(Q.processed));
fprintf('*   %.3d errors           *\n',sum(Q.errors));
fprintf('*   %.3d skipped          *\n',length(plist)-sum(Q.processed)-sum(Q.errors));
fprintf('**************************\n');
fprintf(top_log,'\nFinished. %d processed, %d errors, %d skipped.\n',sum(Q.processed),sum(Q.errors),length(plist)-sum(Q.processed)-sum(Q.errors));
fclose(top_log);

%% File organization
function make_processing_file()
    my_proc_file=fullfile(startDir,'Processing',sprintf('%d',my_comp_num));
    fid=fopen(my_proc_file,'w');
    fclose(fid);
    assignin('caller','my_proc_file',my_proc_file);
end
function delete_processing_file()
    delete(my_proc_file);
    assignin('caller','my_proc_file',[]);
end
end

%% Main preprocessing pipeline   
function EEG=preprocess(o,log_file,tstart)
% Elements adapted from FASTER and MARA pipelines
try
    %%%%%%%%%%%%%%%%
    % File options %
    %%%%%%%%%%%%%%%%
    % 1 File name including full path (string)
    % 2 Reference channel (integer > 0)
    % 3 Number of data channels (integer > 0)
    % 4 Number of extra channels (integer > 0)
    % 5 Channel locations file including full path (string)
    % 6 Save options (cell)
    %%%%%%%%%%%%%%%%
    fullfilename = o.file_options.current_file;
    eeg_chans = o.channel_options.eeg_chans;
    ext_chans = o.channel_options.ext_chans;
    do_reref = o.channel_options.do_avg_reref;
    if do_reref; ref_chan = []; else; ref_chan=o.channel_options.initialRef; end
    [filepath,filename,extension] = fileparts(fullfilename);
    
    if o.epoch_options.segmentByEvent
    if length(o.epoch_options.eventNames) ~= length(o.epoch_options.eventFlags)
        EEG=[];
        fprintf('ERROR: The number of event names does not match the number of event flags.\n');
        fprintf(log_file,'%.2f - The number of event names does not match the number of event flags. Cannot process.\n',toc(tstart));
        return
    end
    end

    c=clock;
    months={'Jan' 'Feb' 'Mar' 'Apr' 'May' 'Jun' 'Jul' 'Aug' 'Sep' 'Oct' 'Nov' 'Dec'};
    fprintf(log_file,'\n%d/%s/%d %d:%d:%d\n',c(3),months{c(2)},c(1),c(4),c(5),round(c(6)));
    fprintf(log_file,'%.2f - Opened log file.\n',toc(tstart));
    
    
    %%%%%%%%%%%%%%
    % File setup %
    %%%%%%%%%%%%%%
    fprintf('Loading %s.\n',fullfilename);
    if strcmpi(extension,'.set')
        EEG = pop_loadset('filename',[filename '.set'],'filepath',filepath);
    elseif strcmpi(extension,'.raw')
        EEG = pop_readegi(fullfilename);
    elseif strcmpi(extension,'.edf') || strcmpi(extension,'.bdf')
        isActive = plugin_askandinstall('Biosig', 'sopen');
        if ~isActive; EEG = []; fprintf('Biosig plugin missing. Cannot process.\n');
            fprintf(log_file,'%.2f - Biosig plugin missing; cannot process without it. Use EEGLAB Extension Manager to install.\n',toc(tstart));
            return
        else
            EEG = pop_biosig(fullfilename,'ref',1);
        end
    elseif strcmpi(extension,'.sma')
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        gain = inputdlg('Enter relative gain (1/2^{12}*[V_{max}-V_{min}]*10^6/gain):',...
            'Enter the relative gain',[1 75],{'0'},options);
        [EEG,~] = pop_snapread(filename, gain);
    elseif strcmpi(extension,'.cnt')
        EEG = pop_loadcnt(fullfilename,'dataformat','auto');
    else % File extension not recognized
        EEG=[];
        fprintf('ERROR: Unknown file format.\n');
        fprintf(log_file,'%.2f - Unknown file format. Cannot process.\n',toc(tstart));
        return
    end
    fprintf(log_file,'%.2f - Loaded file %s.\n',toc(tstart),fullfilename);
    if ~isempty(o.file_options.output_folder_name)
        filepath=o.file_options.oplist{o.file_options.current_file_num};
        if ~exist([filepath filesep 'Intermediate'],'dir')
            mkdir([filepath filesep 'Intermediate']);
        end
    else
        filepath=o.file_options.oplist{o.file_options.current_file_num};            
        if ~exist([filepath filesep 'Intermediate'],'dir')
            mkdir([filepath filesep 'Intermediate']);
        end
        delete(fullfilename);
        if exist([fullfilename(1:end-4) '.fdt'],'file')
            delete([fullfilename(1:end-4) '.fdt']);
        end
        if exist([fullfilename(1:end-4) '.dat'],'file')
            delete([fullfilename(1:end-4) '.dat']);
        end
    end
    EEG.filename = [filename '.set'];
    EEG = eeg_checkset(EEG);
    if EEG.nbchan ~= length(eeg_chans) + length(ext_chans)
        eeg_chans=1:EEG.nbchan;
        options.Resize='on';
        options.WindowStyle='normal';
        options.Interpreter='tex';
        userExtChans = inputdlg([num2str(EEG.nbchan) ' channels were found. '...
            'If any EEG channels are external (EKG, nasion, etc.), enter their indices here:'],...
             'Indicate any external channels',[1 50],{''},options);
        if ~isempty(userExtChans); ext_chans=str2num(userExtChans{1}); end
    end
    % Check if channel locations exist, and if not load them from disk.
    if (~isfield(EEG.chanlocs,'X') || ~isfield(EEG.chanlocs,'Y') || ~isfield(EEG.chanlocs,'Z') || isempty(EEG.chanlocs)) || isempty([EEG.chanlocs(:).X]) || isempty([EEG.chanlocs(:).Y]) || isempty([EEG.chanlocs(:).Z])
        fprintf('Warning: Channel locations file not found\n');
        promptMessage = sprintf('Are you inputting data from the EGI GSN 128-channel EEG cap?');
        buttonText = questdlg(promptMessage, '', 'Yes', 'No', 'Yes');
        if strcmpi(buttonText, 'Yes')
            tempChanLocs=load('GSN128_chan_locs.mat');
            EEG.chanlocs=tempChanLocs.chanlocs; EEG.chaninfo=tempChanLocs.chaninfo; EEG.urchanlocs=tempChanLocs.urchanlocs;
            fprintf(log_file,'%.2f - GSN128 channel locations loaded.\n',toc(tstart));
        else % Not inputting data from the EGI GSN 128-channel EEG cap
            fprintf('Please select channel locations file\n');
            [channame, chanpath] = uigetfile2('*.sfp*;*.sph*;*.loc*;*.locs;*.ced;*.xyz*;*.asc*;*.polhemus*;*.besa*;*.chanedit;*.custom',...
                'Select channel locations file', 'multiselect', 'off');
            drawnow;
            if channame==0 || chanpath==0
                fprintf('ERROR: Channel locations file required.\n');
                fprintf(log_file,'%.2f - Channel locations file required. Cannot process.\n',toc(tstart));
                return
            else
                EEG.chanlocs = readlocs([chanpath channame]);
                fprintf(log_file,'%.2f - Loaded channel locations file from %s.\n',toc(tstart),[chanpath channame]);
            end
        end
        EEG = eeg_checkset(EEG);
        EEG.saved='no';
    end

    
    %%%%%%%%%%%%%%%%
    % Save options %
    %%%%%%%%%%%%%%%%
    EEG = pop_saveset(EEG,'filename',[filename '_original.set'],'filepath',filepath,'savemode','onefile');
    save_before_filter = 1;
    save_before_interp = 1;
    save_before_ica_rej = 1;
    save_before_segment = 1;

    if save_before_filter
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_filt_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['1_pre_filt_' filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end


    %%%%%%%%%%%%%
    % Filtering %
    %%%%%%%%%%%%%
    do_hipass=o.filter_options.hpf_on;
    do_lopass=o.filter_options.lpf_on;
    do_notch=o.filter_options.notch_on;


    if any(any(isnan(EEG.data)))
        fprintf('NaN in EEG data before filtering.\n');
    end
    
    % Initial reference before channel processing
    EEG = h_pop_reref(EEG,o.channel_options.initialRef,'keepref','on');
    EEG.ref = ['Channel ' num2str(EEG.chanlocs(o.channel_options.initialRef).labels)];
    
    if do_hipass
        w_h=o.filter_options.hpf_freq;
        t_h=0.5;  if t_h>w_h; t_h=w_h; end
        r_h=0.05;
        a_h=80;
        
        [m, wtpass, wtstop] = pop_firpmord([w_h-(t_h) w_h+(t_h)], [0 1], [10^(-1*abs(a_h)/20) (10^(r_h/20)-1)/(10^(r_h/20)+1)], EEG.srate);
        if mod(m,2);m=m+1;end
        EEG = pop_firpm(EEG, 'fcutoff', w_h, 'ftrans', t_h, 'ftype', 'highpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Highpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc(tstart),w_h,t_h,m);
    end

    if do_lopass
        w_l=o.filter_options.lpf_freq;
        t_l=2.5;
        r_l=0.01;
        a_l=40;

        [m, wtpass, wtstop] = pop_firpmord([w_l-(t_l) w_l+(t_l)], [1 0], [(10^(r_l/20)-1)/(10^(r_l/20)+1) 10^(-1*abs(a_l)/20)], EEG.srate);
        if mod(m,2);m=m+1;end
        EEG = pop_firpm(EEG, 'fcutoff', w_l, 'ftrans', t_l, 'ftype', 'lowpass', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
        EEG.saved='no';
        fprintf(log_file,'%.2f - Lowpass filter: %.3fHz, transition band: %.2f, order: %d.\n',toc(tstart),w_l,t_l,m);
    end

    if do_notch
        for n=1:length(o.filter_options.notch_freq)
            w_n=[o.filter_options.notch_freq(n)-1.5 o.filter_options.notch_freq(n)+1.5];
            t_n=1;
            r_n=0.05;
            a_n=80;

            [m, wtpass, wtstop] = pop_firpmord([w_n(1)-(t_n) w_n(1)+(t_n) w_n(2)-(t_n) w_n(2)+(t_n)], [0 1 0], [10^(-1*abs(a_n)/20) (10^(r_n/20)-1)/(10^(r_n/20)+1) 10^(-1*abs(a_n)/20)], EEG.srate);
            if mod(m,2);m=m+1;end
            EEG = pop_firpm(EEG, 'fcutoff', w_n, 'ftrans', t_n, 'ftype', 'bandstop', 'wtpass', wtpass, 'wtstop', wtstop, 'forder', m);
            EEG.saved='no';
            fprintf(log_file,'%.2f - Notch filter: %.3f to %.3fHz, transition band: %.2f, order: %d.\n',toc(tstart),w_n(1),w_n(2),t_n,m);
        end
    end

    if save_before_interp
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_interp_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['2_pre_interp_' filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%
    % Channel interpolation %
    %%%%%%%%%%%%%%%%%%%%%%%%%
    chans_to_interp=[];
    o.channel_options.rejection_options.measure=[1 1 1];
    o.channel_options.rejection_options.z=[3 3 3];

    if any(eeg_chans(:)==o.channel_options.initialRef)
        list_properties = channel_properties(EEG,eeg_chans,o.channel_options.initialRef); channel_callback=list_properties;
    elseif any(ext_chans(:)==o.channel_options.initialRef)
        list_properties = channel_properties(EEG,eeg_chans,[]); channel_callback=list_properties;
    else
        error('ERROR: Reference channel not found in neither the EEG nor external channel sets');
    end
    lengths = min_z(list_properties,o.channel_options.rejection_options);
    chans_to_interp = eeg_chans(logical(lengths));
    chans_to_interp = setdiff(chans_to_interp,ref_chan);
    EOG_chans=o.channel_options.eog_chans;
    chans_to_interp = setdiff(chans_to_interp,EOG_chans); % EOG channels necessary for IC algorithm
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Manually verify and search for bad channels %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    color_range=cell(length(eeg_chans)+length(ext_chans),1); color_range(:)={'b'};
    color_range(chans_to_interp)={'r'};
    color_range(ext_chans)={'k'};
    eegplot(EEG.data,'srate',EEG.srate,'eloc_file',EEG.chanlocs,...
    'dispchans',30,'spacing',50,'color',color_range,...
        'title','Manual inspection of EEG data for bad channels');
    options.Resize='on';
    options.WindowStyle='normal';
    options.Interpreter='tex';
    answer = inputdlg('Bad channels:',...
             'Visual data inspection',[1 50],{num2str(chans_to_interp)},options);
    if isempty(answer); chans_to_interp = []; else; chans_to_interp = str2num(answer{1}); end
    while any(chans_to_interp > (length(eeg_chans)+length(ext_chans))) || any(chans_to_interp <= 0)
        answer = inputdlg(sprintf(['Bad channels:\nChannels must be greater than 0 and less than ' num2str(length(eeg_chans)+length(ext_chans))]),...
            'Visual data inspection',[1 50],{num2str(chans_to_interp)},options);
        chans_to_interp = str2num(answer{1});
    end
    close(gcf);
    fprintf('Resuming processing...\n');
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Pre-epoch data for ICA processing %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    EEGpreEpoch=EEG;
    if length(size(EEG.data)) < 3
    oldname = EEG.setname;
    EEG = make_epochs(EEG,1);
    EEG.setname = oldname;
    EEG.saved='no'; 
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Epoch rejection %
    %%%%%%%%%%%%%%%%%%%
    EEGtemp = h_pop_reref(EEG, [], 'exclude', ext_chans, 'refstate', ref_chan);
    o.epoch_options.rejection_options.measure=[1 1 1];
    o.epoch_options.rejection_options.z=[3 3 3];
    if size(EEGtemp.data,3) > 1
        list_properties = epoch_properties(EEGtemp,setdiff(eeg_chans,chans_to_interp));
        [lengths] = min_z(list_properties,o.epoch_options.rejection_options);
        EEG=pop_rejepoch(EEG, find(lengths),0);
        fprintf(log_file,['%.2f - Rejected %d epoch(s) at positions(s) ' ...
            regexprep(num2str(find(lengths)'),'\s+',', ') '.\n'],toc(tstart),length(find(lengths)));
        EEG.saved='no';
    end
    
    %%%%%%%%%%%%%%%%%%%%%%
    % Selected reference %
    %%%%%%%%%%%%%%%%%%%%%%
    if do_reref
        EEG = h_pop_reref(EEG, [], 'exclude',[ext_chans chans_to_interp], 'refstate', ref_chan);
    end
    if o.channel_options.mastoidRef
        mastoidchans = o.channel_options.mastoidChannels;
        mastoidRef = mean(EEG.data(mastoidchans,:),1);
        EEG.data = EEG.data - mastoidRef;
        EEG.data = EEG.data(setdiff(1:EEG.nbchan,mastoidchans),:);
        EEG.chanlocs = EEG.chanlocs(setdiff(1:EEG.nbchan,mastoidchans));
        EEG.nbchan = EEG.nbchan - 2;
        EEG.ref = 'Mastoid';
        EEG = eeg_checkset(EEG);
    end
    
    %%%%%%%%%%%%%%%
    % ICA options %
    %%%%%%%%%%%%%%%
    do_ica = 1;
    k_value = 25;
    do_component_rejection = 1;
    ica_chans = eeg_chans;
    o.ica_options.rejection_options.measure=[1 1 1 1 1];
    o.ica_options.rejection_options.z=[3 3 3 3 3];
    o.ica_options.IC_images=1;

    %%%%%%%%%%%
    % Run ICA %
    %%%%%%%%%%%
    if do_ica && isempty(EEG.icaweights)
        num_pca = min(floor(sqrt(size(EEG.data(:,:),2) / k_value)),(size(EEG.data,1) - length(chans_to_interp) - 1));
        num_pca = min(num_pca,length(setdiff(ica_chans,chans_to_interp)));
        
        ica_chans=intersect(setdiff(ica_chans,chans_to_interp),union(eeg_chans,ext_chans));
        EEG = pop_runica(EEG,'dataset',1, 'chanind',setdiff(ica_chans,chans_to_interp),'options',{'extended',1,'pca',num_pca});    
        
        EEG.saved='no';
        fprintf(log_file,'%.2f - Ran ICA.\n',toc(tstart));
    end

    if save_before_ica_rej
        EEGBAK=EEG;
        EEGBAK.setname = ['pre_comp_rej_' EEG.setname];
        pop_saveset(EEGBAK,'filename',['3_pre_comp_rej_' filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
        clear EEGBAK;
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    Component rejection     %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if do_component_rejection && ~isempty(EEG.icaweights)
        EEG = eeg_checkset(EEG);
        original_name=EEG.setname;
        if do_lopass
            list_properties = component_properties(EEG,EOG_chans,[w_l-(t_l/2) w_l+(t_l/2)]);
        else
            list_properties = component_properties(EEG,EOG_chans);
            o.ica_options.rejection_options.measure(2)=0;
        end
        [lengths] = min_z(list_properties,o.ica_options.rejection_options);
        
        % Plot ICs with suggested rejections for manual inspection
        EEG.reject.gcompreject=lengths;
        [EEG,bad_comps]=pop_selectcomps_integrated(EEG,1:length(lengths));
        
        % Plot components
        if (o.ica_options.IC_images)
            if ~exist([filepath filesep 'Component maps'],'dir')
                mkdir([filepath filesep 'Component maps']);
            end
            p=1;
            activations=eeg_getica(EEG);
            perc_vars = var(activations(:,:),[],2);
            perc_vars = 100*perc_vars./sum(perc_vars);
            for u=1:size(EEG.icawinv,2)
                if ~mod(u-1,16)
                    if (u~=1)
                        saveas(h,sprintf('%s%sComponent maps%sComponents_%d.png',filepath,filesep,filesep,p));
                        p=p+1;
                        close(h);
                    end
                    h=figure;
                end
                subplot(4,4,1+mod(u-1,16));
                topoplot(EEG.icawinv(:,u),EEG.chanlocs(EEG.icachansind));
                title(sprintf('Component %d\n%.1f%% variance',u,perc_vars(u)));
                if ismember(u,bad_comps)
                    c=get(h,'Children');
                    c2=get(c(1),'Children');
                    set(c2(5),'FaceColor',[0.6 0 0]);
                    x=get(c2(5),'XData');
                    x(1:end/2)=1.5*(x(1:end/2));
                    set(c2(5),'XData',x);
                    y=get(c2(5),'YData');
                    y(1:end/2)=1.5*(y(1:end/2));
                    set(c2(5),'YData',y);
                end
            end
            saveas(h,sprintf('%s%sComponent maps%sComponents_%d.png',filepath,filesep,filesep,p));
            if ~isempty(h)
                close(h);
            end
        end
    elseif ~isempty(EEG.icawinv) && o.ica_options.IC_images
        activations=eeg_getica(EEG);
        perc_vars = var(activations(:,:),[],2);
        perc_vars = 100*perc_vars./sum(perc_vars);
        p=1;
        for u=1:size(EEG.icawinv,2)
            if ~mod(u-1,16)
                if (u~=1)
                    saveas(h,sprintf('%s%sComponent maps%sComponents_%d.png',filepath,filesep,filesep,p));
                    p=p+1;
                    close(h);
                end
                h=figure;
            end
            subplot(4,4,1+mod(u-1,16));
            topoplot(EEG.icawinv(:,u),EEG.chanlocs);
            title(sprintf('Component %d\n%.1f%% variance',u,perc_vars(u)));
        end
        saveas(h,sprintf('%s%sComponent maps%sComponents_%d.png',filepath,filesep,filesep,p));
        if ~isempty(h)
            close(h);
        end
    end
    
    % Revert data to match event tags
    EEGtemp=EEG; EEG=EEGpreEpoch;
    EEG.icachansind=EEGtemp.icachansind; EEG.icawinv=EEGtemp.icawinv;
    EEG.icasphere=EEGtemp.icasphere; EEG.icaweights=EEGtemp.icaweights;
    clear EEGtemp;
    % Reject bad components
    if do_component_rejection && ~isempty(EEG.icaweights)
        if ~isempty(bad_comps)
            fprintf('Rejecting components');
            fprintf(' %d',bad_comps);
            fprintf('.\n');
            EEG = pop_subcomp(EEG, bad_comps, 0);
            fprintf(log_file,['%.2f - Rejected %d component(s): ' ...
                regexprep(num2str(bad_comps'),'\s+',', ') '.\n'],toc(tstart),length(bad_comps));
        else
            fprintf('Rejected no components.\n');
            fprintf(log_file,'%.2f - Rejected no components.\n',toc(tstart));
        end
        EEG.setname=original_name;
        EEG.saved='no';
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Complete channel interpolation %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~isempty(chans_to_interp)
        EEG = h_eeg_interp_spl(EEG,chans_to_interp,ext_chans);
        EEG.saved='no';
        fprintf(log_file,['%.2f - Interpolated channel(s) ' ...
            regexprep(num2str(chans_to_interp),'\s+',', ') '.\n'],toc(tstart));
        zs=channel_callback-repmat(mean(channel_callback,1),size(channel_callback,1),1);
        zs=zs./repmat(std(zs,[],1),size(channel_callback,1),1);
        for l=1:length(chans_to_interp)
            cha=chans_to_interp(l);
            chaCorr=zs(cha,1); chaVar=zs(cha,2); chaHurst=zs(cha,3);
            if chaCorr < -3; reasoning='Low mean correlation';
            elseif chaCorr > 3; reasoning='High mean correlation'; % should never happen
            elseif chaVar < -3; reasoning='Low variance';
            elseif chaVar > 3; reasoning='High variance';
            elseif chaHurst < -3; reasoning='Low Hurst exponent';
            elseif chaHurst > 3; reasoning='High Hurst exponent';
            else; reasoning='Manually marked for rejection';
            end
            fprintf(log_file,['   %d: ',reasoning,'.\n'],chans_to_interp(l));
        end
    end
    
    
    %%%%%%%%%%%%%%%%%%%%%%%
    % Unsegmented dataset %
    %%%%%%%%%%%%%%%%%%%%%%%
    tStart=EEG.xmin; tEnd=EEG.xmax;
    fprintf(log_file,['%.2f - ' filename '_full_preprocessed.set'],toc(tstart));
    fprintf(log_file,'.\n');
    EEGBAK=pop_select(EEG,'time',[tStart tEnd]);
    EEGBAK = eeg_checkset(EEGBAK);
    EEGBAK.setname = ['full_preprocessed_' EEG.setname];
    if o.epoch_options.epoch_on
        oldname = EEGBAK.setname;
        if EEGBAK.xmax > o.epoch_options.epoch_length
            EEGBAK = make_epochs(EEGBAK,o.epoch_options.epoch_length);
            fprintf(log_file,'%.2f - Epoched data every %.2f seconds.\n',toc(tstart),o.epoch_options.epoch_length); 
        else
            fprintf(log_file,'%.2f - Data not epoched (segment length less than %.2f seconds).\n',toc(tstart),o.epoch_options.epoch_length); 
        end
        EEGBAK.setname = oldname;
    end
    
    if size(EEGBAK.data,3) > 1
        EEGtemp = EEGBAK;
        if do_reref
            EEGtemp = h_pop_reref(EEGtemp, [], 'exclude',ext_chans, 'refstate', ref_chan);
        end
        fprintf(log_file,'Initial baseline variance: %.2f.\n',median(var(mean(EEGtemp.data(:,1:round(EEGtemp.srate-1*EEGtemp.xmin),:),3),[],2)));
        clear EEGtemp;
    end
    
    %%%%%%%%%%%%%%%%%%%
    % Epoch rejection %
    %%%%%%%%%%%%%%%%%%%
    o.epoch_options.rejection_options.measure=[1 1 1];
    o.epoch_options.rejection_options.z=[3 3 3];
    EEGtemp = h_pop_reref(EEGBAK, [], 'exclude', ext_chans, 'refstate', ref_chan);
    if size(EEGtemp.data,3) > 1
        list_properties = epoch_properties(EEGtemp,setdiff(eeg_chans,chans_to_interp));
        [lengths] = min_z(list_properties,o.epoch_options.rejection_options);
        EEGBAK=pop_rejepoch(EEGBAK, find(lengths),0);     
        fprintf(log_file,['   %.2f - Rejected %d epoch(s) at positions(s) ' ...
            regexprep(num2str(find(lengths)'),'\s+',', ') '.\n'],toc(tstart),length(find(lengths)));
        EEGBAK.saved='no';
    end
    clear EEGtemp

    %%%%%%%%%%%%%%%%%%%%%%%
    % Epoch interpolation %
    %%%%%%%%%%%%%%%%%%%%%%%
    do_epoch_interp=1;
    o.epoch_interp_options.rejection_options.measure=[1 1 1 1];
    o.epoch_interp_options.rejection_options.z=[3 3 3 3];

    if do_epoch_interp && length(size(EEGBAK.data)) > 2
        status = '';
        lengths_ep=cell(1,size(EEGBAK.data,3));
        for v=1:size(EEGBAK.data,3)
            list_properties = single_epoch_channel_properties(EEGBAK,v,eeg_chans);
            lengths_ep{v}=eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options)));
            status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) sprintf('\n')];
        end
        EEGBAK=h_epoch_interp_spl(EEGBAK,lengths_ep,ext_chans);
        EEGBAK.saved='no';
        if ~exist([filepath filesep 'Channel interpolations by epoch'],'dir')
            mkdir([filepath filesep 'Channel interpolations by epoch']);
        end % Directory for epoch interpolation text files
        epoch_interps_log_file=fopen([filepath filesep ...
            'Channel interpolations by epoch' filesep filename ...
            '_full_preprocessed_channel_interpolations_by_epoch.txt'],'a');
        fprintf(epoch_interps_log_file,'%s',status);
        fclose(epoch_interps_log_file);
        fprintf(log_file,'   %.2f - Did per-epoch interpolation cleanup.\n',toc(tstart));
        fprintf(log_file,['   See ' filename '_full_preprocessed' ...
            '_channel_interpolations_by_epoch.txt for details.\n']);
    end
    EEGBAK.urevent=[];
    pop_saveset(EEGBAK,'filename',[filename '_full_preprocessed'],...
        'filepath',filepath,'savemode','onefile');
    clear EEGBAK;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%
    % Segment by event tag %
    %%%%%%%%%%%%%%%%%%%%%%%%
    if o.epoch_options.segmentByEvent
        if isempty(EEG.event)
            fprintf('WARNING: User specified segmentation by task, but no event flags were found in the data.\n');
            fprintf('         All EEG data saved under ''_full_preprocessed.set'' file.\n');
        else
            % Find all matching event flags in EEG data
            eventflags=o.epoch_options.eventFlags;
            allflags=horzcat(eventflags{:});
            allrecorded = {EEG.event.type};
            toDelete=[];
            for i=1:length(allrecorded)
                if any(strcmp(allrecorded(i),allflags))
                else; toDelete = [toDelete i];
                end
            end
            % Delete all unmatched event tags
            EEG = pop_editeventvals(EEG,'delete',toDelete);
            % Check if any events remain
            if isempty(EEG.event)
                fprintf('WARNING: User specified segmentation by task, but no event flags were found in the data.\n');
                fprintf('         All EEG data saved under ''_full_preprocessed.set'' file.\n');
            else % Segment data by event flag
                if save_before_segment
                    EEGBAK=EEG;
                    EEGBAK.setname = ['pre_segment_' EEG.setname];
                    pop_saveset(EEGBAK,'filename',['4_pre_segment_' filename],'filepath',[filepath filesep 'Intermediate'],'savemode','onefile');
                    clear EEGBAK;
                end
                fprintf(log_file,'%.2f - Segmented data by task flags.\n',toc(tstart));
                events=cell(size(eventflags));
                for n=1:length(eventflags)
                    flagVersions = eventflags{n};
                    flagger = [];
                    for m=1:length(flagVersions)
                        flagger = [flagger, ...
                            find(strcmp({EEG.event.type},flagVersions{m}))];
                    end
                    events{n} = sort(flagger);
                end
                eventnames=o.epoch_options.eventNames;
                eventnames=[{'pre-events'},eventnames];
                EEG = pop_editeventvals(EEG,'insert',{1 [] [] []},'changefield',...
                    {1 'type' 'pre-events'},'changefield',{1 'latency' 0});
                for i=1:length(events)
                    events{i}=events{i}+1;
                end
                events=[{1},events];

                % Segment data based on consecutive appearance of flags
                lastEvent=1;
                for i=1:length(eventnames)
                    eventnames{i} = strrep(eventnames{i},' ','_');
                    if ~isempty(events{i}); lastEvent=max(max(events{i}),lastEvent); end
                end
                for i=1:size(events,2)
                    thisEvents=events{i};
                    if isempty(thisEvents); disp([eventnames{i} ' not recorded.']);
                       fprintf(log_file,[eventnames{i} ' data not recorded']);
                       fprintf(log_file,'.\n');
                    else % Events found
                        for j=1:length(thisEvents)
                            tStart=cell2mat({EEG.event(thisEvents(j)).latency})/EEG.srate;
                            if thisEvents(j)==lastEvent; tEnd=EEG.xmax;
                            else; tEnd=cell2mat({EEG.event(thisEvents(j)+1).latency})/EEG.srate;
                            end
                            fprintf(log_file,['%.2f - ' filename '_' eventnames{i} '_' num2str(j) '.set'],toc(tstart));
                            fprintf(log_file,'.\n'); warning off % Disable event latency warnings -- eeglab issue
                            EEGBAK = pop_select(EEG,'time',[tStart tEnd]); w = warning('query','last'); warning on
                            EEGBAK = eeg_checkset(EEGBAK);
                            EEGBAK.setname = [eventnames{i} '_' num2str(j) '_' EEG.setname];
                            if o.epoch_options.epoch_on
                                oldname = EEGBAK.setname;
                                if EEGBAK.xmax > o.epoch_options.epoch_length
                                    EEGBAK = make_epochs(EEGBAK,o.epoch_options.epoch_length);
                                    fprintf(log_file,'%.2f - Epoched data every %.2f seconds.\n',toc(tstart),o.epoch_options.epoch_length); 
                                else
                                    fprintf(log_file,'%.2f - Data not epoched (segment length less than %.2f seconds).\n',toc(tstart),o.epoch_options.epoch_length); 
                                end
                                EEGBAK.setname = oldname;
                            end

                            %%%%%%%%%%%%%%%%%%%
                            % Epoch rejection %
                            %%%%%%%%%%%%%%%%%%%
                            o.epoch_options.rejection_options.measure=[1 1 1];
                            o.epoch_options.rejection_options.z=[3 3 3];
                            EEGtemp = h_pop_reref(EEGBAK, [], 'exclude', ext_chans, 'refstate', ref_chan);
                            if size(EEGtemp.data,3) > 1
                                list_properties = epoch_properties(EEGtemp,setdiff(eeg_chans,chans_to_interp));
                                [lengths] = min_z(list_properties,o.epoch_options.rejection_options);
                                EEGBAK=pop_rejepoch(EEGBAK, find(lengths),0);     
                                fprintf(log_file,['   %.2f - Rejected %d epoch(s) at positions(s) ' ...
                                    regexprep(num2str(find(lengths)'),'\s+',', ') '.\n'],toc(tstart),length(find(lengths)));
                                EEGBAK.saved='no';
                            end

                            %%%%%%%%%%%%%%%%%%%%%%%
                            % Epoch interpolation %
                            %%%%%%%%%%%%%%%%%%%%%%%
                            do_epoch_interp=1;
                            o.epoch_interp_options.rejection_options.measure=[1 1 1 1];
                            o.epoch_interp_options.rejection_options.z=[3 3 3 3];

                            if do_epoch_interp && length(size(EEGBAK.data)) > 2
                                if ~exist([filepath filesep 'Channel interpolations by epoch'],'dir')
                                    mkdir([filepath filesep 'Channel interpolations by epoch']);
                                end % Directory for epoch interpolation text files
                                status = '';
                                lengths_ep=cell(1,size(EEGBAK.data,3));
                                for v=1:size(EEGBAK.data,3)
                                    list_properties = single_epoch_channel_properties(EEGBAK,v,eeg_chans);
                                    lengths_ep{v}=eeg_chans(logical(min_z(list_properties,o.epoch_interp_options.rejection_options)));
                                    status = [status sprintf('%d: ',v) sprintf('%d ',lengths_ep{v}) sprintf('\n')];
                                end
                                EEGBAK=h_epoch_interp_spl(EEGBAK,lengths_ep,ext_chans);
                                EEGBAK.saved='no';
                                epoch_interps_log_file=fopen([filepath filesep ...
                                    'Channel interpolations by epoch' filesep filename ...
                                    ['_' eventnames{i} '_' num2str(j) '_channel_interpolations_by_epoch.txt']],'a');
                                fprintf(epoch_interps_log_file,'%s',status);
                                fclose(epoch_interps_log_file);
                                fprintf(log_file,'   %.2f - Did per-epoch interpolation cleanup.\n',toc(tstart));
                                fprintf(log_file,['   See ' filename '_' eventnames{i} '_' ...
                                    num2str(j) '_channel_interpolations_by_epoch.txt for details.\n']);
                            end
                            remainingFlags = {EEGBAK.event.type};
                            toDelete=[];
                            for a=1:length(remainingFlags)
                                if strcmp(remainingFlags(a),'X')
                                else; toDelete = [toDelete a];
                                end
                            end
                            EEGBAK = pop_editeventvals(EEGBAK,'delete',toDelete); EEGBAK.urevent=[];
                            pop_saveset(EEGBAK,'filename',[filename '_' eventnames{i} '_' num2str(j) '_preprocessed'],...
                                'filepath',filepath,'savemode','onefile');
                            clear EEGBAK;
                        end % Iterate through all events within event type
                    end % Events found
                end % Iterate through all event types
            end % ~isempty(EEG.event) after discarding extraneous flags
        end % ~isempty(EEG.event)
    end % if o.epoch_options.segmentByEvent
    
    fprintf('Done with file %s.\nProcessing time: %d seconds.\n',[filepath filesep filename extension],toc(tstart));
    fprintf(log_file,'%.2f - Finished.\n',toc(tstart));
    fclose(log_file);
    
catch
    m=lasterror;
    EEG_state{1}=evalc('disp(EEG)');
    try
        if ~isempty(fopen(log_file))
            frewind(log_file);
            EEG_state{2}=fscanf(log_file,'%c',inf);
            
            try fclose(log_file); catch; end
        end
    catch
    end
    EEG_state{3}=o;
    EEG_state{4}=builtin('version');
    if exist('eeg_getversion','file')
        EEG_state{5}=eeg_getversion;
    else
        EEG_state{5}=which('eeglab');
    end
    
    assignin('caller','EEG_state',EEG_state);
    rethrow(m);
end

%% Display components and label them for artifact rejection
function [EEG,bad_comps] = pop_selectcomps_integrated(EEG, compnum)

PLOTPERFIG = 35;
bad_comps=[];

for index = 1:PLOTPERFIG:length(compnum)
    [~,b]=pop_selectcomps_helper(EEG,compnum(index:min(length(compnum),index+PLOTPERFIG-1)),index);
    bad_comps=[bad_comps;b];
end
end

function [EEG,bad_comps] = pop_selectcomps_helper(EEG,compnum,index)

COLREJ = '[1 0.6 0.6]';
COLACC = '[0.75 1 0.75]';
    
try
    icadefs
catch
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end

% set up the figure
column =ceil(sqrt( length(compnum) ))+1;
rows = ceil(length(compnum)/column);
if ~exist('fig','var')
	figure('name', [ 'Reject components by map (dataset: ' EEG.setname ')'], 'tag', 'ADEC - Plot', ...
		   'numbertitle', 'off', 'color', BACKCOLOR);
	set(gcf,'MenuBar', 'none');
	pos = get(gcf,'Position');
	set(gcf,'Position', [pos(1) 20 800/7*column 600/5*rows]);
    incx = 120;
    incy = 110;
    sizewx = 100/column;
    if rows > 2
        sizewy = 90/rows;
	else 
        sizewy = 80/rows;
    end
    pos = get(gca,'position'); % plot relative to current axes
	hh = gca;
	q = [pos(1) pos(2) 0 0];
	s = [pos(3) pos(4) pos(3) pos(4)]./100;
	axis off;
end

% figure rows and columns
if EEG.nbchan > 64
    plotelec = 0;
else
    plotelec = 1;
end
count = 1;
for ri = compnum
	if exist('fig','var')
        button = findobj('parent', fig, 'tag', ['comp' num2str(ri)]);
        if isempty(button) 
            error('ERROR: figure does not contain the component button');
        end
    else
		button = [];
    end	
		 
	if isempty( button )
		% compute coordinates
		X = mod(count-1, column)/column * incx-10;  
        Y = (rows-floor((count-1)/column))/rows * incy - sizewy*1.3;  

		% plot the head
        if ~strcmp(get(gcf, 'tag'), 'ADEC - Plot')
            figure(findobj('tag', 'ADEC - Plot'));
        end
		ha = axes('Units','Normalized', 'Position',[X Y sizewx sizewy].*s+q);
        if plotelec
            topoplot(EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                     'off', 'style' , 'fill', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        else
            topoplot(EEG.icawinv(:,ri), EEG.chanlocs, 'verbose', ...
                     'off', 'style' , 'fill','electrodes','off', 'chaninfo', EEG.chaninfo, 'numcontour', 8);
        end
		axis square;

		% plot the button
        if ~strcmp(get(gcf, 'tag'), 'ADEC - Plot')
            figure(findobj('tag', 'ADEC - Plot'));
        end
        guidata(gcf,EEG);
		button = uicontrol(gcf, 'Style', 'pushbutton', 'Units','Normalized', 'Position',...
                   [X Y+sizewy sizewx sizewy*0.25].*s+q, 'tag', ['comp' num2str(ri)]);
		set(button, 'Callback', {@callback_compbutton, ri});
    end
	set( button, 'backgroundcolor', eval(fastif(EEG.reject.gcompreject(ri), COLREJ,COLACC)), 'string', int2str(ri));
    drawnow;
	count = count +1;
end

function callback_compbutton(hObject,eventdata, pos) 
compstatus=popup_comp_prop(guidata(gcf), 0, pos, gcbo, { 'freqrange', [1 50] });
EEG.reject.gcompreject(pos) = compstatus;
%save changes
openplots = findall(0, 'tag', 'ADEC - Plot');
for k = 1: length(openplots)
    guidata(openplots(k), EEG);
end
end


% draw the bottom button
if ~exist('fig','var')
    if ~strcmp(get(gcf, 'tag'), 'ADEC - Plot')
        figure(findobj('tag', 'ADEC - Plot'));
    end
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Cancel', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[-10 -10  15 sizewy*0.25].*s+q, 'callback', 'close(gcf); fprintf(''Operation cancelled\n'')' );
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Set threhsolds', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[10 -10  15 sizewy*0.25].*s+q, 'callback', 'pop_icathresh(EEG); pop_selectcomps(EEG, gcbf);' );
	if isempty(EEG.stats.compenta), set(hh, 'enable', 'off'); end
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See comp. stats', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[30 -10  15 sizewy*0.25].*s+q, 'callback',  ' ' );
	if isempty(EEG.stats.compenta), set(hh, 'enable', 'off'); end
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'See projection', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[50 -10  15 sizewy*0.25].*s+q, 'callback', ' ', 'enable', 'off'  );
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'Help', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[70 -10  15 sizewy*0.25].*s+q, 'callback', 'pophelp(''pop_selectcomps'');' );
    command = sprintf('close(gcf);');
	hh = uicontrol(gcf, 'Style', 'pushbutton', 'string', 'OK', 'Units','Normalized', 'backgroundcolor', GUIBUTTONCOLOR, ...
			'Position',[90 -10  15 sizewy*0.25].*s+q, 'callback', command);%'Callback', {@callback_OK, compnum});
        
end
uiwait(gcf);
bad_comps = find(EEG.reject.gcompreject(compnum))+index-1;
end

%% Epoch data at regular intervals
function EEG = make_epochs(EEG, recurrence)

if recurrence < 0 || recurrence > EEG.xmax
  error('recurrence interval out of bounds');
end
limits = [0 recurrence];

if length(limits) ~= 2 || limits(2) <= limits(1) 
   error('epoch limits must be a 2-vector [minsec maxsec]')
end

% Calculate number of events to add
nu = floor(EEG.xmax/recurrence);

if nu < 1
  error('specified recurrence interval too long')
end

eplength = limits(2)-limits(1);
fprintf('The input dataset will be split into %d epochs of %g s\n',nu,eplength);
fprintf('Epochs will overlap by%2.0f%%.\n',(eplength-recurrence)/eplength*100);

nevents = length(EEG.event);
nurevents = length(EEG.urevent);
for k = 1:nu
   if rem(k,40)
      fprintf('.')
   else
      fprintf('%d',k)
   end
   if k==40 || ( k>40 && ~rem(k-40,70))
     fprintf('\n');
   end

   EEG.event(nevents+k).type = 'X';
   EEG.event(nevents+k).latency = recurrence*(k-1)*EEG.srate+1;

   EEG.urevent(nurevents+k).type = 'X';
   EEG.urevent(nurevents+k).latency = recurrence*(k-1)*EEG.srate+1;
   EEG.event(nevents+k).urevent = nurevents+k;
end
fprintf('\n');

% Sort the events based on their latency
EEG = pop_editeventvals( EEG, 'sort', {'latency' 0}); 

% Split the dataset into epochs
setname = sprintf('%s - %g-s epochs', EEG.setname, recurrence);
EEG = pop_epoch( EEG, { 'X' }, limits, 'newname', ...
                                  setname, 'epochinfo', 'yes');
end

%% Install missing toolboxes
function installRes = plugin_askandinstall(pluginName, pluginFunc, forceInstall)

if nargin < 3, forceInstall = false; end
if nargin < 2 || ~exist(pluginFunc)    
    if ~forceInstall
        installRes = 0;
        % check if deactivated
        try PLUGINLIST = evalin('base', 'PLUGINLIST'); catch, PLUGINLIST = []; end;
        if ~isempty(PLUGINLIST) && isfield(PLUGINLIST, 'plugin')
            indPlugin = find(strcmpi(pluginName,{PLUGINLIST.plugin}));
            if ~isempty(indPlugin) && strcmpi(PLUGINLIST(indPlugin(1)).status, 'deactivated')
                res = questdlg2( [ pluginName ' extension is de-activated. Do you want to reactivate it now?' ], [ pluginName ' extension installation' ], 'No', 'Yes', 'Yes' );
                if strcmpi(res, 'no'), return, end
                plugin_reactivate(PLUGINLIST(indPlugin(1)).foldername);
                evalin('base', 'eeglab rebuild');
                installRes = 1;
                return
            end
        end
        % check if needs to be installed
        res = questdlg2([pluginName ' extension is not installed. Do you want to download it now?'],[pluginName ' extension installation'],'No','Yes','Yes');
    else
        res = 'yes';
    end
    
    if strcmpi(res, 'no'); return; end
    plugins = plugin_getweb('import', []);
    indPlugin = find(strcmpi(pluginName,{plugins.name}));
    if isempty(indPlugin)
        plugins = plugin_getweb('process', []);
        indPlugin = find(strcmpi(pluginName,{plugins.name}));
        if isempty(indPlugin)
            error([ pluginName ' extension not found' ]); 
        end
    end
    result = plugin_install(plugins(indPlugin(1)).zip, plugins(indPlugin(1)).name, plugins(indPlugin(1)).version, forceInstall);
    if result == 1, installRes = 1; end
    evalin('base', 'eeglab rebuild'); close;
else
    installRes = 1;
end
end

end