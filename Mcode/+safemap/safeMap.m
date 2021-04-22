function varargout = safeMap(function_handle, inputs, config)
%safeMap stores on disk partial results freeing memory and avoiding data loss
%
% varargout = safemap.safeMap(function_handle, inputs, config)
%
% The passed one parameter function is fed with each input cell and the outputs
% are saved in a file as they are generated; a cell matrix containing the outputs,
% with the same structure of the input one, is returned.
%
% If something interrupts the execution, `safeMap` will automatically resume once 
% executed again.
%
%
% Arguments
%   - function_handle: a function handle (@nameOfYourFunction or @(x) expression)
%   - inputs: indicates the `function_handle` inputs, can be a cell array, a
%       matrix or a scalar (`n` means `1:n`)
%   - config: a struct with some configuration
%       - filePath: absolute or relative path of the file where the output data  
%           is saved. NB if using a file already existing, make sure the mat 
%           version is 7.3 to allow partial writing [default is `output.mat`]
%       - tempFilePath: absolute or relative path of the file where the current
%           state is saved (the index of the cycle) [default is `progress.safemap.mat`]
%       - variableName: the names of the variables under which the output data is
%           saved inside the file specified in filePath. Can be a single string,
%           a cell of strings if multiple outputs. If less names than outputs
%           "data_i" will be used for the missing names [default is `data`]
%       - returnData: true/false. If true, a cell containing all the outputs 
%           will be returned and the output file will be deleted. 
%           If false the output file handler is returned [default is `false`]
%       - joinUniformOutput: like `cellfun` uniformOutput option, instead of 
%           returning a cell matrix of individual outputs, tries to join them 
%           together in a single matrix. Only numeric, logical, cell and struct 
%           outputs can be joined [default is `true`]
%       - numberOfFunctionOutputs: number of outputs to expect from the passed 
%           function [by default it is inferred from `function_handle`]
%
% Outputs
% By default, safeMap returns a matlab.io.MatFile instance, that is a file handler 
% that allow access to its variables or to parts of them (e.g. `file.data(i,:)`)
%
% If `config.returnData` is `true`, safeMap will return as many cell matrices as 
% the number of outputs, each with the same dimensions of `inputs` (or more with 
% `config.joinUniformOutput`).
%
% In case `config.joinUniformOutput` is `true` (which is the default), each output 
% will be equal to the function outputs stacked along the leading dimensions 
% (e.g. if inputs is of size `[3, 3]` and each output is of size `[4, 2]` the 
% final output will be of size `[3, 3, 4, 2]`).
%

  %% sanitize inputs and eventually use defaults
  % function_handle
  if ~exist('function_handle', 'var') || ~isa(function_handle, 'function_handle')
    error('The first argument must be a function handle. Check "help safeMap" for more info.');
  end
  % inputs
  if ~exist('inputs', 'var')
    error('No input.');
  elseif ~isa(inputs, 'cell')
    % if a scalar n is passed it is interpreted as 1:n
    if isscalar(inputs)
      inputs = 1:inputs;
    end
    % a matrix is transformed element-wise into a cell matrix
    inputs = num2cell(inputs);
  end
  % assert(~isscalar(inputs), 'Inputs should have more than one value');
  % config
  if ~exist('config', 'var')
    config = struct();
  end
  if ~isfield(config, 'filePath')
    config.filePath = 'output.mat';
  end
  if ~isfield(config, 'tempFilePath')
    config.tempFilePath = 'progress.safemap.mat';
  end
  if ~isfield(config, 'returnData')
    config.returnData = false;
  end
  if ~isfield(config, 'joinUniformOutput')
    config.joinUniformOutput = true;
  end
  if isfield(config, 'variableName') && ~isa(config.variableName, 'cell')
    config.variableName = {config.variableName};
  end

  %% if not provided, infer `config.numberOfFunctionOutputs`
  if ~isfield(config, 'numberOfFunctionOutputs')
    if config.returnData
      % if returning data instead of the handle we can infer outputs from safeMap `nargout`
      config.numberOfFunctionOutputs = nargout;
    else
      if isfield(config, 'variableName')
        % otherwise, if passed from the user, use number of strings in `config.variableName`
        config.numberOfFunctionOutputs = length(config.variableName);
      elseif nargout(function_handle) == 0
        error('The passed function declares no output');
      elseif nargout(function_handle) < 0 
        % if infer from `function_handle` nargout is not possible (< 0 means that 
        % there is a varargout) warn the user and use matlab guess (can be bigger)
        config.numberOfFunctionOutputs = abs(nargout(function_handle));
        warning('Cannot determine exactly the number of wanted outputs for the given function, if it is different from %d please provide the `config.numberOfFunctionOutputs` or the `config.variableName` parameter.', config.numberOfFunctionOutputs);
      else
        % infer from `function_handle` nargout
        config.numberOfFunctionOutputs = nargout(function_handle);
      end
    end
  end

  %% check or generate variable names
  if ~isfield(config, 'variableName')
    if config.numberOfFunctionOutputs == 1
      config.variableName = {'data'};
    else
      % variable names: data_1, data_2, ..., data_N with N being `config.numberOfFunctionOutputs`
      config.variableName = arrayfun(@(n) sprintf('data_%d', n), 1:config.numberOfFunctionOutputs, 'un', 0);
    end
  end
  namesNumber_ = length(config.variableName);
  if namesNumber_ < config.numberOfFunctionOutputs
    % fill names with 'data_i', ..., 'data_N' if there are less than outputs
    config.variableName(namesNumber_+1:config.numberOfFunctionOutputs) = ...
      arrayfun(@(n) sprintf('data_%d', n), namesNumber_+1:config.numberOfFunctionOutputs, 'un', 0);
  end

  %% files handling, checks and initializations
  % here is when the current index info is kept. The file will be deleted when finished
  temp = matfile(config.tempFilePath, 'Writable', true);
  resumeIdxFound_ = any(strcmp('i', who(temp)));
  % % hide temp file on windows
  % if ispc 
  %   fileattrib(config.tempFilePath, '+h');
  % end
  if ~isfile(config.filePath)
    % if no data file already found, initialise it
    % make sure mat file is 7.3 version to enable partial reading/writing
    save(config.filePath, 'config', '-v7.3');
    m = matfile(config.filePath, 'Writable', true);
    if resumeIdxFound_
      warning('Resume index found in temporary file, but no data file found, ignoring it and starting from scratch.');
    end
    initNeeded = true;
    startIndex = 1;
  else
    % otherwise verify data file contents
    m = matfile(config.filePath, 'Writable', true);
    mFileVariables_ = who(m);
    variablesInFile_ = cellfun(@(name) any(strcmp(name, mFileVariables_)), config.variableName);
    if resumeIdxFound_
      % check we can resume
      % the file should already contain all the variables where we'll save data
      assert(all(variablesInFile_), 'Tried to resume but not all the output variables have been found in the data file');
      if config.joinUniformOutput
        % make sure file variable size begins with the input size 
        if ndims(inputs) > 2
          checkSizeFn_ = @(name) startsWith(char(size(m, name)), char(size(inputs))); % length(strfind(size(m, name), size(inputs))) && strfind(size(m, name), size(inputs))==1
        else
          % we'll just check later cause we don't know exact size until seeing first output
          checkSizeFn_ = @(name) true;
        end
      else
        % make sure file variable size matches the input size 
        checkSizeFn_ = @(name) isequal(size(m, name), size(inputs));
      end
      sizeCorrespond_ = cellfun(checkSizeFn_, config.variableName);
      assert(all(sizeCorrespond_), 'The data file contains one or more variables with the same name but uncoherent dimensions');
      assert(isscalar(temp.i), 'Temporary file index is not scalar');
      initNeeded = false;
      startIndex = temp.i + 1;
      fprintf('Resuming previous work, from %d-th input\n', startIndex);
    else
      % starting from scratch
      % the file should not contain the variables where we'll save data
      assert(~any(variablesInFile_), 'The data file already contains one or more variables with the given name');
      initNeeded = true;
      startIndex = 1;
    end
  end
  % initialise variables inside the data file 
  % it will be done later when the output sizes are known if joinUniformOutput is true
  if initNeeded && ~config.joinUniformOutput
    for j = 1:config.numberOfFunctionOutputs
      m.(config.variableName{j}) = cell(size(inputs));
    end
  end
  % initialize cells to store the correct number of indexs for the file variables
  inputIndexes = cell(1, ndims(inputs));
  % initialize cells to store the correct number of outputs from function_handle
  outputChunks = cell(1, config.numberOfFunctionOutputs);

  %% execution cycle
  % print intial state  
  stateTextLength = fprintf('Progress: %d / %d (%.1f%%)\n', startIndex-1, numel(inputs), 100 * (startIndex-1) / numel(inputs));
  iterationTimes = nan(numel(inputs), 1);
  for i = startIndex:numel(inputs)
    timer = tic;
    % execute the function
    [outputChunks{:}] = function_handle(inputs{i});
    % compute all input dimensions indexs from i
    [inputIndexes{:}] = ind2sub(size(inputs), i);
    
    if config.joinUniformOutput
      % only first iteration
      if i == startIndex
        inputDims = zeros(1, config.numberOfFunctionOutputs);
        outputIndexes = cell(1, config.numberOfFunctionOutputs);
        initOutputIndexes = cell(1, config.numberOfFunctionOutputs);
        for j = 1:config.numberOfFunctionOutputs
          [variableSize, inputDims(j)] = uniformVariableSize(inputs, outputChunks{j});
          if initNeeded
            % check output data type is ok
            if ~iscell(outputChunks{j}) && ~isnumeric(outputChunks{j}) && ~islogical(outputChunks{j}) && ~isstruct(outputChunks{j})
              error('You enabled uniform output, but the %d-th output is not a cell nor a matrix', j);
            end
            % compute initial matrix size (partially empty)
            initSize = [variableSize(1:inputDims(j)) zeros(1, length(variableSize)-inputDims(j))];
            % initialize file storage
            m.(config.variableName{j}) = uniformInit(outputChunks{j}, initSize);
          else
            % if resuming, just check the sizes of variables in storage is [inputSize outputSize]
            if ~isequal(size(m, config.variableName{j}), variableSize)
              error('The %d-th output size is [%s], which is not coherent with variable %s in data file of size [%s]', j, num2str(variableSize), config.variableName{j}, num2str(size(m, config.variableName{j})));
            end
          end
          % precompute the indexs to slice the stored variable 
          [initOutputIndexes{j}, outputIndexes{j}] = outputIndexing(outputChunks{j});
        end
      end

      % save output
      for j = 1:config.numberOfFunctionOutputs
        if initNeeded && i == startIndex
          idxs = [inputIndexes(1:inputDims(j)) initOutputIndexes{j}];
        else
          idxs = [inputIndexes(1:inputDims(j)) outputIndexes{j}];
        end
        leadingDims = inputDims(j) - ndims(outputChunks{j}) + length(outputIndexes{j});
        m.(config.variableName{j})(idxs{:}) = shiftdim(outputChunks{j}, -leadingDims); 
      end
    else
      % save output
      for j = 1:config.numberOfFunctionOutputs
        m.(config.variableName{j})(inputIndexes{:}) = outputChunks(j);
      end
    end
    
    % save the current state
    temp.i = i;

    % print state info
    iterationTimes(i) = toc(timer);
    cancelText_ = repmat('\b', 1, stateTextLength);
    stateText_ = sprintf('Progress: %d / %d (%.1f%%)\nEstimated remaining time: %.0fs\n', ...
      i, numel(inputs), 100 * i / numel(inputs), ...
      mean(iterationTimes, 'omitnan') * (numel(inputs) - i)...
    );
    fprintf([cancelText_ '%s'], stateText_);
    stateTextLength = length(stateText_);
  end

  % delete the temporary file
  delete(temp.Properties.Source);
  % prepare return
  if config.returnData
    varargout = arrayfun(@(index) m.(config.variableName{index}), 1:config.numberOfFunctionOutputs, 'UniformOutput', false);
    % delete data file
    delete(m.Properties.Source);
  else
    varargout = {m};
  end
end

function compactSize = computeCompactSize(x)
  if isscalar(x)
    compactSize = [];
  elseif isvector(x)
    compactSize = length(x);
  else
    compactSize = size(x);
  end
end

function [variableSize, inputDims] = uniformVariableSize(inputs, outputChunk)
  inputCompactSize = computeCompactSize(inputs);
  outputCompactSize = computeCompactSize(outputChunk);
  totalLen = length(inputCompactSize) + length(outputCompactSize);
  if totalLen < 2
    % edge cases
    if size(inputs, 2) > 1
      % [1 n] [1 1] -> [1 n], 2
      variableSize = [1 inputCompactSize];
    else
      % [n 1] [1 1] -> [n 1], 2       % [1 1] [1 1] -> [1 1], 2
      % [1 1] [m 1] -> [1 m], 1       % [1 1] [1 m] -> [1 m], 1
      variableSize = [inputCompactSize ones(1, 2 - totalLen) outputCompactSize];
    end
    inputDims = 2 - length(outputCompactSize);
  else
    variableSize = [inputCompactSize outputCompactSize];
    inputDims = length(inputCompactSize);
  end
end

function [initOutputIndexes, outputIndexes] = outputIndexing(outputChunk)
  % create text cell with colons to index correct matrix location while saving
  % and array cell with ranges to index it the first time
  compactSize = computeCompactSize(outputChunk);
  initOutputIndexes = cell(1, length(compactSize));
  for i = 1:length(compactSize)
    initOutputIndexes{i} = 1:compactSize(i);
  end
  outputIndexes = repmat({':'}, 1, length(compactSize));
end

function [initData] = uniformInit(outputChunk, initSize)
  if isnumeric(outputChunk)
    initData = zeros(initSize, class(outputChunk));
    if ~isreal(outputChunk)
      initData = complex(initData);
    end
  elseif islogical(outputChunk)
    initData = false(initSize);
  elseif iscell(outputChunk)
    initData = cell(initSize);
  elseif isstruct(outputChunk)
    fieldnames_ = fieldnames(outputChunk);
    initData = cell2struct(cell([length(fieldnames_) initSize]), fieldnames_);
  else
    error('Cannot handle data of type %s with `config.joinUniformOutput` set, try disabling it.', class(outputChunk));
  end
end