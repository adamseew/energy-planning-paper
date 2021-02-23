
% MAIN


%% selection

 
[indx] = listdlg ...
    (   'PromptString', {'Select a plan (refer to the figures)'}, ...
        'SelectionMode', 'single', 'ListString', ...
        { ...
            'Plan i (dynamic)', ...
            'Plan ii (dynamic)', ...
            'Plan I (static)', ...
            'Plan II (static)'
        }, ... 
        'ListSize', [220 80] ...
    );

if indx == 1
    
    SIM8_revised
    
elseif indx == 2
    
    SIM9_revised
    
elseif indx == 3
    
    SIM11
    
elseif indx == 4
    
    SIM10
    
end

