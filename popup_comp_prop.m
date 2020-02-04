function outstat = popup_comp_prop(EEG, typecomp, chanorcomp, winhandle, spec_opt)

if nargin < 1
	help popup_comp_prop;
	return;   
end
if nargin < 5
	spec_opt = {};
end
if nargin == 1
	typecomp = 1;    % defaults
        chanorcomp = 1;
end
if typecomp == 0 && isempty(EEG.icaweights)
   error('No ICA weights recorded for this dataset -- first run ICA on it');
end   
if nargin == 2
	promptstr    = { fastif(typecomp,'Channel index(ices) to plot:','Component index(ices) to plot:') ...
                     'Spectral options (see spectopo() help):' };
	inistr       = { '1' '''freqrange'', [2 50]' };
	result       = inputdlg2( promptstr, 'Component properties - popup_comp_prop()', 1,  inistr, 'popup_comp_prop');
	if size( result, 1 ) == 0; return; end
   
	chanorcomp   = eval( [ '[' result{1} ']' ] );
    spec_opt     = eval( [ '{' result{2} '}' ] );
end

% plotting several component properties
% -------------------------------------
if length(chanorcomp) > 1
    for index = chanorcomp
        popup_comp_prop(EEG, typecomp, index, 0, spec_opt);  % call recursively for each chanorcomp
    end
    return;
end

if chanorcomp < 1 || chanorcomp > EEG.nbchan
   error('Component index out of range');
end   

% assumed input is chanorcomp
try 
    icadefs
catch 
	BACKCOLOR = [0.8 0.8 0.8];
	GUIBUTTONCOLOR   = [0.8 0.8 0.8]; 
end
basename = [fastif(typecomp,'Channel ', 'Component ') int2str(chanorcomp) ];

fhandle = figure('name', ['popup_comp_prop() - ' basename ' properties'], ...
    'color', BACKCOLOR, 'numbertitle', 'off', 'visible', 'off');
pos     = get(fhandle,'Position');
set(fhandle,'Position', [pos(1) pos(2)-500+pos(4) 500 500], 'visible', 'on');
hh = axes('parent',fhandle);
pos      = get(hh,'position'); % plot relative to current axes
q = [pos(1) pos(2) 0 0];
s = [pos(3) pos(4) pos(3) pos(4)]./100;
axis(hh,'off');

% plot topoplot
h = axes('parent',fhandle,'Units','Normalized', 'Position',[-10 60 40 42].*s+q);

if isfield(EEG.chanlocs, 'theta')
    if typecomp == 1 % plot single channel locations
        topoplot( chanorcomp, EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                 'electrodes','off', 'style', 'blank', 'emarkersize1chan', 12); axis square;
    else             % plot component map
        topoplot( EEG.icawinv(:,chanorcomp), EEG.chanlocs, 'chaninfo', EEG.chaninfo, ...
                 'shading', 'interp', 'numcontour', 3); axis square;
    end
else
    axis(h,'off');
end
basename = [fastif(typecomp,'Channel ', 'IC') int2str(chanorcomp) ];
title(basename, 'fontsize', 14);

% plot time series data & spectogram
hhh = axes('Parent', fhandle,'Units','Normalized', 'Position',[45 62 48 38].*s+q);
eeglab_options; 
if EEG.trials > 1
    % put title at top of time series data & spectogram
    axis(hhh,'off');
    hh = axes('Parent', fhandle,'Units','Normalized', 'Position',[45 62 48 38].*s+q);
    EEG.times = linspace(EEG.xmin, EEG.xmax, EEG.pnts);
    if EEG.trials < 6
      ei_smooth = 1;
    else
      ei_smooth = 3;
    end
    if typecomp == 1 % plot channel
           printf('time_series_image() received channel data rather than component data');
    else % plot component
         icaacttmp  = eeg_getdatact(EEG, 'component', chanorcomp);
         offset     = nan_mean(icaacttmp(:));
%          era        = nan_mean(squeeze(icaacttmp)')-offset;
%          era_limits = get_era_limits(era);
         gcapos=get(gca,'Position');
         time_series_image(icaacttmp-offset, EEG.times*1000, gcapos);%, '', ...
             %'caxis', 2/3, 'cbar','erp','eegylabel','','vltg_ticks',era_limits);   
    end
    axes(hhh);
    title(sprintf('%s activity \\fontsize{10}', basename), 'fontsize', 14);
else
      if typecomp == 1 % plot channel
           printf('time_series_image() received channel data rather than component data');  
      else % plot component
            printf('time_series_image() received unepoched data');  
      end
end	

% plot spectrum
if ~exist('winhandle')
    winhandle = NaN;
end
if ishandle(winhandle) 
	h = axes('Parent', fhandle,'units','normalized', 'position',[5 10 95 35].*s+q);
else
	h = axes('Parent', fhandle,'units','normalized', 'position',[5 0 95 40].*s+q);
end

try
	eeglab_options; 
	if typecomp == 1
		[spectra freqs] = spectopo( EEG.data(chanorcomp,:), EEG.pnts, EEG.srate, spec_opt{:} );
	else 
		if option_computeica  
			[spectra freqs] = spectopo( EEG.icaact(chanorcomp,:), EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
        else
    		icaacttmp = (EEG.icaweights(chanorcomp,:)*EEG.icasphere)*reshape(EEG.data(EEG.icachansind,:,:), length(EEG.icachansind), EEG.trials*EEG.pnts); 
			[spectra freqs] = spectopo( icaacttmp, EEG.pnts, EEG.srate, 'mapnorm', EEG.icawinv(:,chanorcomp), spec_opt{:} );
		end
	end
    % set up new limits
	set( get(h, 'ylabel'), 'string', 'Power 10*log_{10}(\muV^{2}/Hz)', 'fontsize', 14); 
	set( get(h, 'xlabel'), 'string', 'Frequency (Hz)', 'fontsize', 14); 
	title('Average power spectrum', 'fontsize', 14); 
catch err
	axis off;
    text(0.1, 0.3, [ 'Error: no spectrum plotted' 10 err.message 10]);
end	
	
% display buttons
if ishandle(winhandle)
	COLREJ = '[1 0.6 0.6]';
	COLACC = '[0.75 1 0.75]';
	% CANCEL button
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Cancel', 'Units','Normalized','Position',[-10 -10 15 6].*s+q, 'callback', 'close(gcf);');

	% VALUE button
	hval  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'Values', 'Units','Normalized', 'Position', [15 -10 15 6].*s+q);

	% REJECT/ACCEPT button
    if ~isempty(EEG.reject.gcompreject)
    	status = EEG.reject.gcompreject(chanorcomp);
    else
        status = 0;
    end
    guidata(gcf,status);
	hr = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', eval(fastif(status,COLREJ,COLACC)), ...
				'string', fastif(status, 'REJECT', 'ACCEPT'), 'Units','Normalized', 'Position', [40 -10 15 6].*s+q, 'userdata', status, 'tag', 'rejstatus');
	command = [ 'set(gcbo, ''userdata'', ~get(gcbo, ''userdata''));' ...
				'if get(gcbo, ''userdata''),' ...
				'     set( gcbo, ''backgroundcolor'',' COLREJ ', ''string'', ''REJECT'');' ...
				'else ' ...
				'     set( gcbo, ''backgroundcolor'',' COLACC ', ''string'', ''ACCEPT'');' ...
				'end' ];					
	set(hr, 'callback', command); 

	% HELP button
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'backgroundcolor', GUIBUTTONCOLOR, 'string', 'HELP', 'Units','Normalized', 'Position', [65 -10 15 6].*s+q, 'callback', 'pophelp(''popup_comp_prop'');');

	% OK button
 	command = [ 'global EEG;' ...
 				'tmpstatus = get( findobj(''parent'', gcbf, ''tag'', ''rejstatus''), ''userdata'');' ...
 				'tmpfile = fopen(''tmp.txt'',''a''); fprintf(tmpfile,''%d'',tmpstatus);'];
            
	if winhandle ~= 0
        if VERS < 8.04
            command = [ command ...
                sprintf('if tmpstatus set(%3.15f, ''backgroundcolor'', %s); else set(%3.15f, ''backgroundcolor'', %s); end; ', ...
                winhandle, COLREJ, winhandle, COLACC)];
        elseif VERS >= 8.04
            command = [ command ...
                sprintf('if tmpstatus set(findobj(''Tag'', ''%s''), ''backgroundcolor'', %s); else set(findobj(''Tag'',''%s''), ''backgroundcolor'', %s); end; ', ...
                winhandle.Tag, COLREJ, winhandle.Tag, COLACC)];
        end
    end
	command = [ command ' fclose(tmpfile); close(gcf); clear tmpstatus' ];
	h  = uicontrol(fhandle, 'Style', 'pushbutton', 'string', 'OK', 'backgroundcolor', GUIBUTTONCOLOR, 'Units','Normalized', 'Position',[90 -10 15 6].*s+q, 'callback', command);

	% draw the figure for statistical values
	index = num2str( chanorcomp );
	command = [ ...
		'figure(''MenuBar'', ''none'', ''name'', ''Statistics of the component'', ''numbertitle'', ''off'');' ...
		'' ...
		'pos = get(gcf,''Position'');' ...
		'set(gcf,''Position'', [pos(1) pos(2) 340 340]);' ...
		'pos = get(gca,''position'');' ...
		'q = [pos(1) pos(2) 0 0];' ...
		's = [pos(3) pos(4) pos(3) pos(4)]./100;' ...
		'axis off;' ...
		''  ...
		'txt1 = sprintf(''(\n' ...
						'Entropy of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ' AND                 \t\t\t----\n\n' ...
					    'Kurtosis of component activity\t\t%2.2f\n' ...
					    '> Rejection threshold \t\t%2.2f\n\n' ...
					    ') OR                  \t\t\t----\n\n' ...
					    'Kurtosis distibution \t\t\t%2.2f\n' ...
					    '> Rejection threhold\t\t\t%2.2f\n\n' ...
					    '\n' ...
					    'Current thesholds sujest to %s the component\n\n' ...
					    '(after manually accepting/rejecting the component, you may recalibrate thresholds for future automatic rejection on other datasets)'',' ...
						'EEG.stats.compenta(' index '), EEG.reject.threshentropy, EEG.stats.compkurta(' index '), ' ...
						'EEG.reject.threshkurtact, EEG.stats.compkurtdist(' index '), EEG.reject.threshkurtdist, fastif(EEG.reject.gcompreject(' index '), ''REJECT'', ''ACCEPT''));' ...
		'' ...				
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-11 4 117 100].*s+q, ''Style'', ''frame'' );' ...
		'uicontrol(gcf, ''Units'',''Normalized'', ''Position'',[-5 5 100 95].*s+q, ''String'', txt1, ''Style'',''text'', ''HorizontalAlignment'', ''left'' );' ...
		'h = uicontrol(gcf, ''Style'', ''pushbutton'', ''string'', ''Close'', ''Units'',''Normalized'', ''Position'', [35 -10 25 10].*s+q, ''callback'', ''close(gcf);'');' ...
		'clear txt1 q s h pos;' ];
	set( hval, 'callback', command); 
	if isempty( EEG.stats.compenta ); set(hval, 'enable', 'off'); end
else
end
uiwait(gcf);
fID=fopen('tmp.txt','r'); fSpec='%d'; 
outstat=fscanf(fID,fSpec);
fclose(fID); delete('tmp.txt');
return;

function out = nan_mean(in)

    nans = find(isnan(in));
    in(nans) = 0;
    sums = sum(in);
    nonnans = ones(size(in));
    nonnans(nans) = 0;
    nonnans = sum(nonnans);
    nononnans = find(nonnans==0);
    nonnans(nononnans) = 1;
    out = sum(in)./nonnans;
    out(nononnans) = NaN;

    
function era_limits=get_era_limits(era)
%function era_limits=get_era_limits(era)
%
% Returns the minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)
%
% Inputs:
% era - [vector] Event related activation or potential
%
% Output:
% era_limits - [min max] minimum and maximum value of an event-related
% activation/potential waveform (after rounding according to the order of
% magnitude of the ERA/ERP)

mn=min(era);
mx=max(era);
mn=orderofmag(mn)*round(mn/orderofmag(mn));
mx=orderofmag(mx)*round(mx/orderofmag(mx));
era_limits=[mn mx];


function ord=orderofmag(val)
%function ord=orderofmag(val)
%
% Returns the order of magnitude of the value of 'val' in multiples of 10
% (e.g., 10^-1, 10^0, 10^1, 10^2, etc ...)
% used for computing erpimage trial axis tick labels as an alternative for
% plotting sorting variable

val=abs(val);
if val>=1
    ord=1;
    val=floor(val/10);
    while val>=1,
        ord=ord*10;
        val=floor(val/10);
    end
    return;
else
    ord=1/10;
    val=val*10;
    while val<1,
        ord=ord/10;
        val=val*10;
    end
    return;
end