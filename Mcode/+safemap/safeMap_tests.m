% Tests for safeMap
%
% To run this, call safemap.safeMap_tests.

% test triangles

% inputs
vector_scalar_inputs = 1:100;
dim3_scalar_inputs = ones(5, 5, 5);
vector_vector_inputs = num2cell(rand(20, 20), 2);
% size(vector_vector_inputs): [20 1]
% size(vector_vector_inputs{1}): [1 20]

% utils
allInFile = @(h, names) all(cellfun(@(n) any(strcmp(n, who(h))), names));
delete progress.safemap.mat

%% Test 1: vector of scalar inputs, scalar output
delete 1.mat
expSize_ = size(vector_scalar_inputs);
config_ = struct('filePath', '1.mat');
output_file = safemap.safeMap(@(x) 1, vector_scalar_inputs, config_);
sz_ = size(output_file, 'data');
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));

%% Test 2: vector of scalar inputs, scalar output without joinUniformOutput
delete 2.mat
expSize_ = size(vector_scalar_inputs);
expChildSize_ = [1 1];
config_ = struct('joinUniformOutput', false, 'filePath', '2.mat');
output_file = safemap.safeMap(@(x) 1, vector_scalar_inputs, config_);
child_ = output_file.data(1, 1);
sz_ = size(output_file, 'data'); 
szc_ = size(child_{1}); 
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));
assert(isequal(szc_, expChildSize_), 'Output cell size wrong: obtained [%s] expected [%s]', num2str(szc_), num2str(expChildSize_));

%% Test 3: tensor of scalar inputs, scalar output
delete 3.mat
expSize_ = size(dim3_scalar_inputs);
config_ = struct('joinUniformOutput', false, 'filePath', '3.mat');
output_file = safemap.safeMap(@(x) 1, dim3_scalar_inputs, config_);
sz_ = size(output_file, 'data');
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));

%% Test 4: tensor of scalar inputs, vector output
delete 4.mat
expSize_ = [size(dim3_scalar_inputs) 10 10 2];
config_ = struct('filePath', '4.mat');
output_file = safemap.safeMap(@(x) repmat(magic(10), 1, 1, 2), dim3_scalar_inputs, config_);
sz_ = size(output_file, 'data');
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));

%% Test 5: tensor of scalar inputs, vector output without joinUniformOutput
delete 5.mat
expSize_ = size(dim3_scalar_inputs);
expChildSize_ = [10 10 2];
config_ = struct('joinUniformOutput', false, 'filePath', '5.mat');
output_file = safemap.safeMap(@(x) repmat(magic(10), 1, 1, 2), dim3_scalar_inputs, config_);
child_ = output_file.data(1,1,1);
sz_ = size(output_file, 'data');
szc_ = size(child_{1});
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));
assert(isequal(szc_, expChildSize_), 'Output cell size wrong: obtained [%s] expected [%s]', num2str(szc_), num2str(expChildSize_));

%% Test 6: vector of vector inputs, 2 scalar outputs, return data
delete 6.mat
expSize_ = size(vector_vector_inputs);
config_ = struct('returnData', true, 'filePath', '6.mat');
[maxValues, maxIndexs] = safemap.safeMap(@max, vector_vector_inputs, config_);
sz1_ = size(maxValues);
sz2_ = size(maxIndexs);
assert(isequal(sz1_, expSize_) && isequal(sz2_, expSize_), 'Output sizes wrong: obtained [%s] and [%s], expected [%s]', num2str(sz1_), num2str(sz2_), num2str(expSize_));

%% Test 7: vector of vector inputs, 2 vector outputs
delete 7.mat
expSize_ = [length(vector_vector_inputs) 5];
config_ = struct('variableName', {{'values', 'indexs'}}, 'filePath', '7.mat');
output_file = safemap.safeMap(@(x) maxk(x, 5), vector_vector_inputs, config_);
sz1_ = size(output_file, 'values');
sz2_ = size(output_file, 'indexs');
assert(allInFile(output_file, {'values', 'indexs'}), 'Variables not found in output file');
assert(isequal(sz1_, expSize_) && isequal(sz2_, expSize_), 'Output sizes wrong: obtained [%s] and [%s] expected [%s]', num2str(sz1_), num2str(sz2_), num2str(expSize_));

%% Test 8: vector of vector inputs, 2 vector outputs without joinUniformOutput, return only first output
delete 8.mat
expSize_ = size(vector_vector_inputs);
expChildSize_ = [1 5];
config_ = struct('returnData', true, 'joinUniformOutput', false, 'filePath', '8.mat');
maxValues = safemap.safeMap(@(x) maxk(x, 5), vector_vector_inputs, config_);
sz_ = size(maxValues);
szc_ = size(maxValues{1});
assert(isequal(sz_, expSize_), 'Output size wrong: obtained [%s] expected [%s]', num2str(sz_), num2str(expSize_));
assert(isequal(szc_, expChildSize_), 'Output cell size wrong: obtained [%s] expected [%s]', num2str(szc_), num2str(expChildSize_));

%% Test 9: error interruption
delete 9.mat
config_ = struct('filePath', '9.mat');
data = true(size(vector_scalar_inputs));
data(50) = false;
try
    output_file = safemap.safeMap(@(i) data(i) || {}, vector_scalar_inputs, config_);
catch
    output_file = safemap.safeMap(@(i) false, vector_scalar_inputs, config_);
end
assert(isequal(output_file.data(1:49), true(1, 49)), '');
assert(isequal(output_file.data(50:end), false(1, length(vector_scalar_inputs)-49)), '');