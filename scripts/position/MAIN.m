%MAIN
% Simulation of the Energy-Aware Dynamic Mission Planning Algorithm


%% Position sim

disp('[P    ] Position sim');

selection = 0;
           

%% Selection of subroutine
 
% Runs if one runs just MAIN instead of selectively running sections...

[indx] = listdlg ...
    (   'PromptString', {'Select a subroutine'}, ...
        'SelectionMode', 'single', 'ListString', ...
        { ...
            'OP1: given TEE generate GDN, plot, animate', ...
            'OP2: ', ...
        }, ... 
        'ListSize', [350 150] ...
    );
selection = 1;

% Now comes these subroutines:


%% Position sim, OP1

% Vector field simulation from a defined TEE, plot and possible animation

if or(selection == 0, indx == 1)

    OP1;
end


%% Position sim, OP2

% 

if or(selection == 0, indx == 2)

    OP2;
end

