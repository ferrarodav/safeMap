# Memory and interruption safe map <!-- omit in toc -->

## Alternative for Matlab® `cellfun` <!-- omit in toc -->

[![View safeMap on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://it.mathworks.com/matlabcentral/fileexchange/90807-safemap)

`safeMap` saves on disk data generated by long and heavy computations so that you don't have to worry about memory problems and unexpected interruptions.

The passed function is fed with each input cell and the outputs are saved in a file as they are generated;`safeMap` can either return the file handle or the whole data (with `config.returnData`, which also deletes the whole file after returning the data). If something interrupts the execution, `safeMap` will automatically resume once executed again.

Each output of the function can be saved in a cell matrix with the same structure of the input one, or joined (with `config.joinUniformOutput`). In the latter case, which is the default, each output matrix will be equal to the function outputs stacked along the leading dimensions (e.g. if inputs is of size `[3, 3]` and each output is of size `[4, 2]` the final output will be of size `[3, 3, 4, 2]`).

Moreover, the used file paths and the name under which each output variable is stored can be specified in the configuration.
<!--When needed, the number of outputs can be manually specified (`config.numberOfFunctionOutputs`) otherwise it will be inferred.-->

```text
>> output = safemap.safeMap(@(x) rand(1,300)*x^5*rand(300,1), repmat({rand(300)}, 100, 100));
Progress: 731 / 10000 (7.3%)
Estimated remaining time: 143s
🔴(Ctrl+C) Operation terminated by user during ...
>> output = safemap.safeMap(@(x) rand(1,300)*x^5*rand(300,1), repmat({rand(300)}, 100, 100));
Resuming previous work, from 732-th input
Progress: 10000 / 10000 (100.0%)
Estimated remaining time: 0s
```

## Limitations

`safeMap` writes data on the disk after every output computed, so it's inherently slower than `cellfun`. Moreover, currently it accepts only one input, though this can be easily overcome like in the example below.

## Table of contents <!-- omit in toc -->

- [Limitations](#limitations)
- [Installation](#installation)
- [Arguments](#arguments)
- [Outputs](#outputs)
- [Examples](#examples)
  - [Explore a grid of hyperparameters for a simulation or training a Neural Network](#explore-a-grid-of-hyperparameters-for-a-simulation-or-training-a-neural-network)
  - [Setting a configuration](#setting-a-configuration)
- [Contributing](#contributing)
- [License](#license)

## Installation

Clone the [repo](https://github.com/ferrarodav/safeMap/) (or [download a copy of the latest code](https://github.com/ferrarodav/safeMap/archive/refs/heads/main.zip)) and add the `Mcode` directory to your MATLAB path with `addpath`.

```bash
git clone https://github.com/ferrarodav/safeMap
```

## Arguments

`safemap.safeMap(function_handle, inputs, config)`
- `function_handle`: a function handle (@nameOfYourFunction or @(x) expression)
- `inputs`: indicates the `function_handle` inputs, can be a cell array, a
    matrix or a scalar (`n` means `1:n`)
- `config`: a struct with some configuration
  - `filePath`: absolute or relative path of the file where the output data is saved. NB if using a file already existing, make sure the mat version is 7.3 to allow partial writing [default is `output.mat`]
  - `tempFilePath`: absolute or relative path of the file where the current state is saved (the index of the cycle) [default is `progress.safemap.mat`]
  - `variableName`: the names of the variables under which the output data is saved inside the file specified in filePath. Can be a single string, a cell of strings if multiple outputs. If less names than outputs "data_i" will be used for the missing names [default is `data`]
  - `returnData`: true/false. If true, a cell containing all the outputs will be returned and **the output file will be deleted**. If false the output file handler is returned [default is `false`]
  - `joinUniformOutput`: like `cellfun` uniformOutput option, instead of returning a cell matrix of individual outputs, tries to join them together in a single matrix. Only numeric, logical, cell and struct outputs can be joined [default is `true`]
  - `numberOfFunctionOutputs`: number of outputs to expect from the passed function [by default it is inferred from `function_handle`]

## Outputs

By default, safeMap returns a matlab.io.MatFile instance, that is a file handler that allow access to its variables or to parts of them (e.g. `file.data(i,:)`).

If `config.returnData` is `true`, safeMap will return as many cell matrices as the number of outputs, each with the same dimensions of `inputs` (or more with `config.joinUniformOutput`).

In case `config.joinUniformOutput` is `true` (which is the default), each output will be equal to the function outputs stacked along the leading dimensions (e.g. if inputs is of size `[3, 3]` and each output is of size `[4, 2]` the final output will be of size `[3, 3, 4, 2]`).

## Examples

### Explore a grid of hyperparameters for a simulation or training a Neural Network

```matlab
data = 1;

% list of hps values 
hps = {{i, 1}, num2cell(1:9), {1, 10, 100}};

% make grid of all possible combinations
all_params = cell(1, length(hps));
[all_params{:}] = ndgrid(hps{:});

% trick to pass multiple inputs since it is currently not supported by safeMap
ins = cellfun(@(a,b,c,d) {a, b, c, d}, all_params{:}, 'un', 0);
output_file = safemap.safeMap(@(vars) data * vars{1} * vars{2} * vars{3}, ins);
% desiderable: safemap.safeMap(@(a,b,c,d) a * b * c, all_params{:});

% read the output corresponding to hps: 1, 4, 100
disp(output_file.data(1,4,3,1)); % 400

% read the output corresponding to hps: i, 3, 10
disp(output_file.data(2,3,2,1)); % 30i
```

### Setting a configuration

```matlab
>> data = num2cell(rand(100, 1000), 2);
>> size(data) 
ans =
    100    1
>> size(data{1}) 
ans =
    1      1000    
>> maxValues, maxIndexs = safemap.safeMap(@(x) maxk(x, 2), data, struct('returnData', true));
>> size(maxValues) 
ans =
    100    2
>> size(maxIndexs) 
ans =
    100    2
>> output_file = safemap.safeMap(@(x) maxk(x, 2), data, struct('filePath', 'max.mat', 'variableNames', {{'values', 'indexs'}}));
>> who -file max.mat

Your variables are:

values  indexs  

```

## Contributing

Pull requests are welcome. For major changes, please open an issue first to discuss what you would like to change.

Here is some contribution idea:

- allow multiple inputs with varargin
- a parameter to save progress to disk less frequently in order to balance between safety, memory occupation and speed of execution
- parallel execution with sharding
- use MATLAB datastore API

To run the test use matlab command:

```matlab
results = runtests('safemap.safeMapTest')
```

## License

[MIT License](https://choosealicense.com/licenses/mit/)
