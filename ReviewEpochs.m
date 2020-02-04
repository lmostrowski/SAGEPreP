% EEG processing pipeline  (Lauren Ostrowski, 12-19-18
%    email lauren_ostrowski@brown.edu with any concerns)
eeglab; close
fprintf('\nSelect appropriate ".set" file to review/reject epochs ...\n');
EEG=pop_loadset;
if ~isempty(EEG)
    if EEG.trials > 1
        macrorej = 'EEG.reject.rejmanual';
        macrorejE = 'EEG.reject.rejmanualE';
        elecrange = (1:EEG.nbchan);
        colrej = EEG.reject.rejmanualcol;
        rej = eval(macrorej);
        rejE = eval(macrorejE);
        icacomp = 1; superpose = 1; reject = 1;
        eeg_rejmacro; % script macro for generating command and old rejection arrays
    end
    eegplotoptions = [eegplotoptions 'submean' 'off'];
    eegplot( EEG.data, 'srate', EEG.srate, 'title', 'Scroll channel activities', ...
			  'limits', [EEG.xmin EEG.xmax]*1000 , 'command', command, ...
              'spacing', 40, eegplotoptions{:});
    uiwait(gcf);
    fprintf(['\nSave the dataset under the "Epochs reviewed" folder '...
        'in the subject folder ...\n']); close
    if ~exist([EEG.filepath 'Epochs reviewed/'],'dir')
        mkdir([EEG.filepath 'Epochs reviewed/']);
    end
    [file,path]=uiputfile('*.set','Save final EEG',[EEG.filepath 'Epochs reviewed/' EEG.filename(1:end-4) '_epochs_reviewed']);
    if ( ischar(file) && ischar(path) )
        EEG = pop_saveset(EEG,'filename',file,'filepath',path,'savemode','onefile');
    end
end; clear