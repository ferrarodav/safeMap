% Tests for safeMap
%
% To run this, call runtests('safemap.safeMapTest').

function tests = safeMapTest
    tests = functiontests(localfunctions);
end

function setupOnce(testCase)
    deleteTempFiles;
end

function setup(testCase) 
    % utility function
    testCase.TestData.allInFile = @(h, names) all(cellfun(@(n) any(strcmp(n, who(h))), names));
    % base config
    testCase.TestData.config = struct('filePath', 'safeMapTestData.mat', 'tempFilePath', 'safeMapTestProgressFile.mat');
    % input data
    testCase.TestData.vector_scalar_inputs = 1:100;
    testCase.TestData.dim3_scalar_inputs = ones(5, 5, 5);
    testCase.TestData.vector_vector_inputs = num2cell(rand(20, 20), 2);
    % size of vector_vector_inputs: [20 1]
    % size of vector_vector_inputs{1}: [1 20]
end

function teardown(testCase)
    deleteTempFiles;
end

function testOne(testCase)
% vector of scalar inputs, scalar output
    input = testCase.TestData.vector_scalar_inputs;
    output_file = safemap.safeMap(@scalarFun, input, testCase.TestData.config);
    verifyEqual(testCase, size(output_file, 'data'), size(input));
end

function testTwo(testCase)
% vector of scalar inputs, scalar output without joinUniformOutput
    input = testCase.TestData.vector_scalar_inputs;
    testCase.TestData.config.joinUniformOutput = false;
	output_file = safemap.safeMap(@scalarFun, input, testCase.TestData.config);
    dataElement = output_file.data(1, 1);
	verifyEqual(testCase, size(output_file, 'data'), size(input));
    verifyEqual(testCase, size(dataElement{1}), [1 1]);
end

function testThree(testCase)
% tensor of scalar inputs, scalar output
    input = testCase.TestData.dim3_scalar_inputs;
    output_file = safemap.safeMap(@scalarFun, input, testCase.TestData.config);
    verifyEqual(testCase, size(output_file, 'data'), size(input));
end

function testFour(testCase)
% tensor of scalar inputs, vector output
    input = testCase.TestData.dim3_scalar_inputs;
    output_file = safemap.safeMap(@tensorFun, input, testCase.TestData.config);
    verifyEqual(testCase, size(output_file, 'data'), [size(input) 10 10 2]);
end

function testFive(testCase)
% tensor of scalar inputs, vector output without joinUniformOutput
    input = testCase.TestData.dim3_scalar_inputs;
    testCase.TestData.config.joinUniformOutput = false;
    output_file = safemap.safeMap(@tensorFun, input, testCase.TestData.config);
    dataElement = output_file.data(1,1,1);
    verifyEqual(testCase, size(output_file, 'data'), size(input));
    verifyEqual(testCase, size(dataElement{1}), [10 10 2]);
end

function testSix(testCase)
% vector of vector inputs, 2 scalar outputs, return data
    input = testCase.TestData.vector_vector_inputs;
    testCase.TestData.config.returnData = true;
    [maxValues, maxIndexs] = safemap.safeMap(@max, input, testCase.TestData.config);
    verifyEqual(testCase, size(maxValues), size(input));
    verifyEqual(testCase, size(maxIndexs), size(input));
end

function testSeven(testCase)
% vector of vector inputs, 2 vector outputs
    input = testCase.TestData.vector_vector_inputs;
    testCase.TestData.config.variableName = {'values', 'indexs'};
    output_file = safemap.safeMap(@max5, input, testCase.TestData.config);
    verifyEqual(testCase, size(output_file, 'values'), [length(input) 5]);
    verifyEqual(testCase, size(output_file, 'indexs'), [length(input) 5]);
end

function testEight(testCase)
% vector of vector inputs, 2 vector outputs without joinUniformOutput, return only first output
    input = testCase.TestData.vector_vector_inputs;
    testCase.TestData.config.returnData = true;
    testCase.TestData.config.joinUniformOutput = false;
    maxValues = safemap.safeMap(@max5, input, testCase.TestData.config);
    verifyEqual(testCase, size(maxValues), size(input));
    verifyEqual(testCase, size(maxValues{1}), [1 5]);
end

function testNine(testCase)
% error interruption
    input = testCase.TestData.vector_scalar_inputs;
    data = true([1, length(input)]);
    false_idx = round(length(input) / 2);
    data(false_idx) = false;
    try
        output_file = safemap.safeMap(@(i) data(i) || {}, input, testCase.TestData.config);
    catch
        output_file = safemap.safeMap(@(i) false, input, testCase.TestData.config);
    end
    verifyEqual(testCase, output_file.data(1,1:false_idx-1), true(1, 49));
    verifyEqual(testCase, output_file.data(1,false_idx:end), false(1, length(input) - false_idx + 1));
end

function y = scalarFun(x)
    y = 1;
end

function y = tensorFun(x)
    y = repmat(magic(10), 1, 1, 2);
end

function [values, indexs] = max5(x)
    [values, indexs] = maxk(x, 5);
end

function deleteTempFiles
    if isfile('safeMapTestProgressFile.mat')
        delete safeMapTestProgressFile.mat;
    end
    if isfile('safeMapTestData.mat')
        delete safeMapTestData.mat;
    end
end