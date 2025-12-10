function BIDS_step04_manual_badchannelappend(fig_path, globalvar_path, output_tsv_path)
    % 1. Load and show interactive figure
    fprintf('Opening figure: %s\n', fig_path);
    % open(fig_path);
    openfig(fig_path, 'new', 'visible');

    figure(gcf); % bring to front
    zoom on;     % enable zoom
    datacursormode on; % enable data cursor for manual inspection

    % 2. Load globalVar struct
    fprintf('Loading globalVar from: %s\n', globalvar_path);
    tmp = load(globalvar_path);
    if ~isfield(tmp, 'globalVar')
        error('globalVar not found in %s', globalvar_path);
    end
    globalVar = tmp.globalVar;

    % 3. Prompt for bad channels
    fprintf('Inspect the plot and input bad channel numbers (e.g., [1 5 19]):\n');
    if exist('promptBadChanSpec', 'file')
        bad_chan_spec = promptBadChanSpec;
    else
        bad_chan_spec = input('Enter bad channel indices: ');
    end

    badChanIndex = bad_chan_spec(:);
    % badChanLabel = globalVar.channame(badChanIndex);
    % reason = repmat({''}, size(badChanIndex));

    T = table(badChanIndex) %, badChanLabel, reason);

    % 5. Write to .tsv
    fprintf('Saving bad channels to TSV: %s\n', output_tsv_path);
    writetable(T, output_tsv_path, 'FileType', 'text', 'Delimiter', '\t');
    disp('TSV saved. You can open and edit reasons later if needed.');

    % % 4. Append and deduplicate
    % globalVar.badChan = unique([globalVar.badChan, bad_chan_spec]);
    % globalVar.badChan_specs.manual = bad_chan_spec;

    % % 5. Save updated globalVar
    % fprintf('Saving updated globalVar to: %s\n', globalvar_path);
    % save(globalvar_path, 'globalVar', '-v7.3');
    % disp('Done.');
end
