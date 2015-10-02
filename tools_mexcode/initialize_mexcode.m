% initialize_mexcode.m
%
% Compiles .c files in the mexcode subdirectory 

func_names = {'fastinterp2', 'spikeconv_mex'}; % Function names to mex

for j =1:length(func_names)
    if ~exist(func_names{j})
        fprintf(1, 'mexing %s.c\n', func_names{j});
        mex([func_names{j} '.c']);
    end
end

clear func_names;
