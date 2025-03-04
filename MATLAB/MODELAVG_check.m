function [method,options] = MODELAVG_check(method,D,y,options)
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LL        AAAA   VV   VV GGGGGGG %%
%% MMM     MM OOOOOOOO DDDDDDDD EEEEEEE LL       AA  AA  VV   VV GG   GG %%
%% MMMM  MMMM OO    OO DD    DD EE      LL       AA  AA  VV   VV GG   GG %%
%% MM MMMM MM OO    OO DD    DD EEEE    LL       AA  AA  VV   VV GGGGGGG %%
%% MM  MM  MM OO    OO DD    DD EEEE    LL       AAAAAA   VV VV  GGGGGGG %%
%% MM      MM OO    OO DD    DD EE      LL      AA    AA   VVV        GG %%
%% MM      MM OOOOOOOO DDDDDDDD EEEEEEE LLLLLLL AA    AA   VVV       GGG %%
%% MM      MM  OOOOOO  DDDDDDD  EEEEEEE LLLLLLL AA    AA    V    GGGGGGG %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%% Checks the input defined by the user                                  %%
%%                                                                       %%
%% SYNOPSIS: [method,options] = MODELAVG_check(method,D,y,options)       %%
%%  where                                                                %%
%%   method    [input] String with model averaging used                  %%
%%   D         [input] nxK matrix forecasts ensemble members             %%
%%   y         [input] nx1 vector with verifying observations            %%
%%   options   [input] structure BMA algorithmic variables               %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%
%%                                                                       %%
%%   (c) Written by Jasper A. Vrugt, Feb 2012                            %%
%%   University of California Irvine                                     %%
%%                                                                       %%
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%

% Derive current time and set deadline
deadline = datenum('28-Feb-2025');

% Now check whether this is a trial version or not
if ( deadline - now ) < 0
    % ERROR -- trial version ended
    error('MODELAVG ERROR: Trial version of MODELAVG V1.1 has ended');
end

% open an output file with warnings
fid = fopen('warning_file.txt','w+');
fprintf(fid,'-------------- MODELAVG warning file --------------\n');

% Check MODELAVG input data structures
if isempty(method)
    error(['MODELAVG ERROR: Input argument method should not be empty ' ...
        'but a string enclosed between quotes with model averaging ' ...
        'method to use']);
elseif ~ischar(method)
    error(['MODELAVG ERROR: Input argument method should be a string ' ...
        'enclosed between quotes']);
end

% Check MODELAVG input data structures
if isempty(D)
    error(['MODELAVG ERROR: Input argument D should not be empty but ' ...
        'a n x K matrix with the n forecasts of K models']);
elseif ~isnumeric(D)
    error(['MODELAVG ERROR: Input argument D should be a n x K matrix ' ...
        '( = numerical values )']);
end

% Check MODELAVG input data structures
if isempty(y)
    error(['MODELAVG ERROR: Input argument y should not be empty but a ' ...
        'n x 1 vector with the verifying observations of the training ' ...
        'data set']);
elseif ~isnumeric(y)
    error(['MODELAVG ERROR: Input argument y should be a n x 1-vector ' ...
        '( = numerical values )']);
end

% Check MODELAVG input data structures
if ~isstruct(options)
    error(['MODELAVG ERROR: Input argument options should be a structure ' ...
        'with fields']);
end

% How many data points and models?
[n,K] = size ( D ); [m,K2] = size( y );

% Now check this
if ( K > n )
    error(['MODELAVG ERROR: Ill-determined problem as you have more ' ...
        'models in the ensemble than their corresponding ' ...
        'predictions/training observations; K > n']);
end

% Now check this
if ( K == 1 )
    error(['MODELAVG ERROR: No point to use model averaging if only ' ...
        'a single model is used ( matrix D only has one column ) ']);
end

if ( n == 1 )
    error(['MODELAVG ERROR: No point to use model averaging if matrix ' ...
        'D stores only a single prediction of each of the K models, ' ...
        'that is n = 1']);
end

% Now check this
if ( K2 > 1 )
    if ( m > 1 )
        error(['MODELAVG ERROR: Input argument y should be a n x 1 ' ...
            'vector with training data - not a matrix']);
    elseif ( m == 1 )
        % Now write to screen
        evalstr = char(['MODELAVG WARNING: Vector y with verifying ' ...
            'observations of training data set should be column vector; ' ...
            'transposed in code \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        % Now tranpose so it is always vertical
        y = y(:); m = K2;
    end
end

% Now check this
if ( m ~= n )
    error(['MODELAVG ERROR: Number of rows, n, of n x K matrix D ' ...
        'and number of elements ( = training data points ) of ' ...
        'vector y do not match']);
end

% First use lower case letters
method = lower ( method );

% Check content method 
if ~any(strcmp(method,{'ewa','bga','aica','bica','gra','bma', ...
        'mma','mma-s'}))
    error(['MODELAVG ERROR: Unknown model averaging method: method ' ...
        'should be ''ewa''/''bga''/''aica''/''bica''/''gra''' ...
        '/''bma''/''mma''/''mma-s'' '])
end

% Which fields has the user specified?
field_names = fieldnames(options);

% How many field_names?
Nf = numel ( field_names );

% If field_names is empty --> nothing has to change
if (Nf == 0)
    % Do nothing
else
    % PDF and VAR should be capital; alpha, p and print lower values
    for i = 1 : Nf
        if strcmpi(field_names(i),'PDF') || strcmpi(field_names(i),'VAR') ...
                || strcmpi(field_names(i),'TAU') ...
                || strcmpi(field_names(i),'CPU')
            % make sure it is upper case (capitals) for PDF and VAR + 
            % CPU (May, 2017) - and lower case for their settings
            eval(char(strcat('temp.',upper(field_names(i)), ...
                ' = lower(options.',field_names(i),');')));
        elseif strcmpi(field_names(i),'p') || strcmpi(field_names(i), ...
                'print') || strcmpi(field_names(i),'alpha') || ...
                strcmpi(field_names(i),'postproc')
            % make sure it is lower case for p, print, postproc and alpha 
            % and lower case for their settings
            eval(char(strcat('temp.',lower(field_names(i)), ...
                ' = lower(options.',field_names(i),');')));
        end
    end
    % Now assign options to be equal to temp
    options = temp; clear temp;
end

% Now check p
if any(strcmp(method,{'aica','bica','mma','mma-s'})) % sum(strncmp(method,{'aica','bica','mma','mma-s'},inf))
    if ~isfield(options,'p')
        error(char(strcat(['MODELAVG ERROR: Unknown complexity of each ' ...
            'model for'],{' '},upper(method),{' '},['method --> Define ' ...
            'field ''p'' of structure options'])));
    end
    if isfield(options,'p')
        if isempty(options.p)
            error(['MODELAVG ERROR: Field ''p'' of structure options ' ...
                'should not be empty but a ( K x 1 )-vector with number ' ...
                'of "free" parameters of each model ( = proxy for model ' ...
                'complexity )']);
        end
        if ~isnumeric(options.p)
            error(['MODELAVG ERROR: Field ''p'' of structure options ' ...
                'should be a ( K x 1 )-vector with number of "free" ' ...
                'parameters of each model ( = proxy for model complexity )']);
        end
        % Determine size of options.p ( -> must be numerical values )
        [r2,d2] = size(options.p);
        % Now check this
        if ( d2 > 1 )
            if ( r2 > 1 )
                error(['MODELAVG ERROR: Input argument options.p should ' ...
                    'not be a matrix but a K x 1 vector with the ' ...
                    'complexity ( = number of "free" parameters ) ' ...
                    'of each model ']);
            end
        end
        if ( d2*r2 ~= K )
            error(['MODELAVG ERROR: Number of elements of field ''p'' ' ...
                'of structure options does not match K, the number of ' ...
                'models of the ensemble in matrix D']);
        end
        % For mma and mma-s -> options.p is vertical; for aica and bica options.p is horizontal
        %if strmatch(method,{'aica','bica'}) % sum(strncmp(method,{'aica','bica'},inf))
        if any(strcmp(method,{'aica','bica'}))        
            options.p = options.p(:)';
        elseif any(strcmp(method,{'mma','mma-s'})) % sum(strncmp(method,{'mma','mma-s'},inf))
            options.p = options.p(:);
        end
    end
elseif any(strcmp(method,{'ewa','bga','gra','bma'})) % sum(strncmp(method,{'ewa','bga','gra','bma'},inf))
    % Now check p again
    if isfield(options,'p')
        % Now write to screen
        evalstr = char(strcat(['MODELAVG WARNING: Field ''p'' of structure ' ...
            'options defined but not used by '],{' '}, ...
            upper(method),{' '},'method \n'));
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
        %         % Now remove field p of options as method does not use it
        %         options = rmfield(options,'p');
    end
end

% Now check whether to print to screen
if ~isfield(options,'print')
    % Now write to screen
    evalstr = char(['MODELAVG WARNING: Field ''print'' of structure ' ...
        'options not defined as ''yes'' or ''no'' ( default of ' ...
        'options.print = ''yes'' used ) \n']);
    % Now print warning to screen and to file
    fprintf(evalstr); fprintf(fid,evalstr);
end
if isfield(options,'print')
    if isempty(options.print)
        error(['MODELAVG ERROR: Field ''print'' of structure options ' ...
            'should not be empty but should contain a string ' ...
            '(''yes'' or ''no'') \n']);
    elseif ~ischar(options.print)
        error(['MODELAVG ERROR: Field ''print'' of structure options ' ...
            'should be a string (content equal to ''yes'' or ''no'') \n']);
    elseif ~any(strcmp(options.print,{'yes','no'})) % ~sum(strncmp(options.print,{'yes','no'},inf))
        error(['MODELAVG ERROR: Field ''print'' of structure options ' ...
            'should equal ''yes'' or ''no'' ']);
    end
end

% Now check whether BMA method is used
if strcmp(method,'bma')
    % Which percentage of ensemble forecasts is negative
    pct = 100 * sum(sum(D<0))/numel(D);
    % Which percentage of training data (y) is negative
    pct_2 = 100 * sum(y<0)/numel(y);
    
    % warning if BMA conditional forecast PDF is not defined
    if ~isfield(options,'PDF')
        error(['MODELAVG ERROR: BMA method is used but conditional ' ...
            'distribution is not defined in field ''PDF'' of structure ' ...
            'options']);
    end
    if isfield(options,'PDF')
        if isempty(options.PDF)
            error(['MODELAVG ERROR: Field ''PDF'' of structure options ' ...
                'should not be empty but should be a string (options: ' ...
                '''normal'' or ''gamma'' or ''gen_normal'') ']);
        end
        if ~ischar(options.PDF)
            error(['MODELAVG ERROR: Field ''PDF'' of structure options ' ...
                'should be a string (options: ''normal'' or ''gamma'' ' ...
                'or ''gen_normal'' ) ']);
        end
        if ~any(strcmp(options.PDF,{'normal','gamma','gen_normal', ...
                'lognormal','weibull','power_law','tnormal','gev','gpareto'}))
            error(['MODELAVG ERROR: Unknown conditonal PDF of BMA method ' ...
                '(use ''normal'' or ''lognormal'' or ''gen_normal'' ' ...
                'or ''gamma'' or ''weibull'' or ''gev'' or ''gpareto'' ']);
        end
% %         if strcmp(options.PDF,'normal')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of normal ' ...
% %                 'conditional distribution is evaluated with built-in ' ...
% %                 'PDF function -> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %         end
% %         if strcmp(options.PDF,'lognormal')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of lognormal ' ...
% %                 'conditional distribution is evaluated with built-in ' ...
% %                 'PDF function -> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key ' ...
% %                 'consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %         end
% %         if strcmp(options.PDF,'weibull')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of Weibull ' ...
% %                 'distribution is evaluated with built-in PDF function ' ...
% %                 '-> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %         end
% %         if strcmp(options.PDF,'gev')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of GEV ' ...
% %                 'distribution is evaluated with built-in PDF function ' ...
% %                 '-> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %         end
% %         if strcmp(options.PDF,'gpareto')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of generalized pareto ' ...
% %                 'distribution is evaluated with built-in PDF function ' ...
% %                 '-> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %         end
        if strcmp(options.PDF,'gamma')
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Density of gamma ' ...
% %                 'conditional distribution is evaluated with built-in ' ...
% %                 'PDF function -> this is not particularly CPU-efficient\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);
% %             % Print that code is slow because built-in PDF function is used
% %             evalstr = char(strcat(['MODELAVG NOTE: Contact author ' ...
% %                 '(Email: jasper@uci.edu) if efficiency is a key consideration\n']));
% %             % Now print warning to screen and to file
% %             fprintf(evalstr); fprintf(fid,evalstr);            
            % Provide warning
            if pct > 0
                evalstr = char(strcat(['MODELAVG WARNING: Gamma ' ...
                    'distribution is used but '],{' '}, ...
                    num2str(round(100*pct)/100),{' '},...
                    '%% of the forecasts of the ensemble are negative \n'));
                % Now print warning to screen and to file
                fprintf(evalstr); fprintf(fid,evalstr);
            end
            % Provide warning
            if pct_2 > 0
                evalstr = char(strcat(['MODELAVG WARNING: Gamma ' ...
                    'distribution is used but '],{' '}, ...
                    num2str(round(100*pct_2)/100),{' '},...
                    '%% of the verifying observations are negative \n'));
                % Now print warning to screen and to file
                fprintf(evalstr); fprintf(fid,evalstr);
            end
        end
        if strcmp(options.PDF,{'gen_normal'}) || strcmp(options.PDF,{'gev'})
            if ~isfield(options,'TAU')
                error(['MODELAVG ERROR: Field ''TAU'' of structure ' ...
                    'options should be specified for generalized ' ...
                    'normal conditional PDF (options: ''1'' or ''2'') ']);
            end
            if isempty(options.TAU)
                error(['MODELAVG ERROR: Field ''TAU'' of structure ' ...
                    'options should not be empty but contain a ' ...
                    'string (options: ''1'' or ''2'') ']);
            end
            if ~ischar(options.VAR)
                error(['MODELAVG ERROR: Field ''TAU'' of structure ' ...
                    'options should be a string (options: ''1'' or ''2'') ']);
            end
            if ~any(strcmp(options.TAU,{'1','2'})) % ~sum(strncmp(options.TAU,{'1','2','3','4'},inf))
                error(['MODELAVG ERROR: Unknown value of field ' ...
                    '''TAU'' of structure options (use ''1''or ''2'') ']);
            end
        end
    end    
    if ~isfield(options,'VAR')
        evalstr = char(['MODELAVG WARNING: BMA method is used but ' ...
            'field ''VAR'' of structure options not defined ' ...
            '( default of options.VAR = ''1'' is used ) \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
    if isfield(options,'VAR')
        if isempty(options.VAR)
            error(['MODELAVG ERROR: Field ''VAR'' of structure options ' ...
                'should not be empty but contain a string (options: ' ...
                '''1'' or ''2'' or ''3'' or ''4'') ']);
        end
        if ~ischar(options.VAR)
            error(['MODELAVG ERROR: Field ''VAR'' of structure options ' ...
                'should be a string (options: ''1'' or ''2'' or ''3'' ' ...
                'or ''4'') ']);
        end
        if ~any(strcmp(options.VAR,{'1','2','3','4'})) % ~sum(strncmp(options.VAR,{'1','2','3','4'},inf))
            error(['MODELAVG ERROR: Unknown value of field ''VAR'' ' ...
                'of structure options (use ''1''or ''2'' or ''3'' ' ...
                'or ''4'') ']);
        end
    end
    % Now check sigma ahead of time
    if any(strcmp(options.VAR,{'3','4'})) % sum(strncmp(options.VAR,{'3','4'},inf))
        % Provide warning
        if pct > 0
            evalstr = char(strcat('MODELAVG WARNING: Variance option', ...
                {' '},options.VAR,{' '},'used but',{' '},...
                num2str(round(100*pct)/100),{' '},['%% of forecasts ' ...
                'of the ensemble are negative \n']));
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        if strcmp(options.VAR,'3')
            evalstr = char(['MODELAVG WARNING: As "sigma(:,k) = ' ...
                'c*D(:,k)" ( c > 0 ) then sigma(:,k) can be < 0 --> ' ...
                'code uses "sigma(:,k) = c*abs(D(:,k))"; k in [1...K] \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
        if strcmp(options.VAR,'4')
            evalstr = char(['MODELAVG WARNING: As "sigma(:,k) = ' ...
                'c(k)*D(:,k)" ( c(k) > 0 ) then sigma(:,k) can be < 0 ' ...
                '-->  code uses "sigma(:,k) = c(k)*abs(D(:,k))"; ' ...
                'k in [1...K] \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end
    if ~isfield(options,'alpha')
        % Now write to screen
        evalstr = char(['MODELAVG WARNING: Field ''alpha'' of ' ...
            'structure options not defined -> default of ' ...
            'alpha = 0.05 is used \n']);
        % Now print warning to screen and to file
        fprintf(evalstr); fprintf(fid,evalstr);
    end
    if isfield(options,'alpha')
        if isempty(options.alpha)
            error(['MODELAVG ERROR: Field ''alpha'' of structure ' ...
                'options should not be empty but should be ' ...
                'a numeric value (options: (0-1)) \n']);
        end
        if ~isnumeric(options.alpha)
            error(['MODELAVG ERROR: Field ''alpha'' of structure ' ...
                'options should be a numeric value (options: (0-1)) ']);
        end
        if options.alpha < 0
            error(['MODELAVG ERROR: Field ''alpha'' of structure ' ...
                'options cannot be negative (e.g. ' ...
                'use 0.10 or 0.05 or 0.01)']);
        end
        if options.alpha > 1
            error(['MODELAVG ERROR: Field ''alpha'' of structure ' ...
                'options cannot be larger than one (e.g. ' ...
                'use 0.10 or 0.05 or 0.01)']);
        end
        if options.alpha > 0.75
            evalstr = char(['MODELAVG WARNING: Field ''alpha'' of ' ...
                'structure options set rather unusual (e.g. ' ...
                'use 0.10 or 0.05 or 0.01) \n']);
            % Now print warning to screen and to file
            fprintf(evalstr); fprintf(fid,evalstr);
        end
    end

end

if ~strcmp(method,'bma')
    % Remove PDF field
    if isfield(options,'PDF')
        options = rmfield(options,'PDF');
    end
    % Remove VAR field
    if isfield(options,'VAR')
        options = rmfield(options,'VAR');
    end
    % Remove TAU field
    if isfield(options,'TAU')
        options = rmfield(options,'TAU');
    end
    % Remove alpha field
%     if isfield(options,'alpha')
%         options = rmfield(options,'alpha');
%     end
end

if strcmp(method,'bma')
    % Remove TAU field
    if ~ sum(strcmp(options.PDF,{'gen_normal','power_law','gev','gpareto'})) 
        if isfield(options,'TAU')
            options = rmfield(options,'TAU');
        end
    end
end

% Now check this one
if isfield(options,'CPU')
    if isempty(options.CPU)
        error(['MODELAVG ERROR: Field ''CPU'' of structure options ' ...
            'should not be empty but should contain a string ' ...
            '(''yes'' or ''no'') \n']);
    elseif ~ischar(options.CPU)
        error(['MODELAVG ERROR: Field ''CPU'' of structure options ' ...
            'should be a string (content equal to ''yes'' or ''no'') \n']);
    elseif ~any(strcmp(options.CPU,{'yes','no'})) % ~sum(strncmp(options.CPU,{'yes','no'},inf))
        error(['MODELAVG ERROR: Field ''CPU'' of structure options ' ...
            'should equal ''yes'' or ''no'' ']);
    end
else % --> set options.CPU = no;
    options.CPU = 'no';
end

if ~isfield(options,'postproc')
    options.postproc = 'yes';
end

% Now close warning_file.txt file
fclose(fid);

end
