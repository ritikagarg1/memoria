function BIDS_step04_badchannelreject(project_name, sub, ses, run, main_dir, ieegsource_dir, event_dir,  diagnostics_dir,pathstr, use_interactive, badchan_file )
    %% Branch 4 - bad channel rejection
    % [refChan, badChan, epiChan, emptyChan] = GetMarkedChans(sub);
     manualVar.refChan = [];
     manualVar.epiChan = [];
     manualVar.emptyChan = []; 
        
    % BadChanRejectCAR_SSD(sub, project_name, block_names, dirs,pathstr)%matlab must be started from anaconda command prompt
    % use_interactive = false;
    % badchan_file = 'nofile'; % or add path .tsv
    BIDS_BadChanRejectCAR_generic(sub, ses, run, project_name, main_dir, ieegsource_dir, ...
    event_dir, diagnostics_dir, manualVar, pathstr, use_interactive, badchan_file)

end