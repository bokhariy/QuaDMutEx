clear;

% load dataset
load test_data/GBM.mat 
% mySparseMatrix is the binary matrix of mutations (patients in rows, genes in columns)
% GenesNames is a cell matrix with gene names
% selectionCDFn is a non-uniform cumulative discribution function for random selection of new genes to be evaluated (selectionCDFu is uniform CDF)

% set parameters
C=0.5; % select larger value of C for larger (more genes) datasets
k=1; 

% run QuaDMutEx
[objective_function,selectedGenes, namesOfSelectedGenes]=QuadMutEx(mySparseMatrix,GenesNames,10000,30,C,k,selectionCDFn)

% calculate solution quality metrics
solutionMetrics=QuadMutExMetricsStruct(selectedGenes,mySparseMatrix,C,k)

% EXPECTED RESULTS: 
% note that QuaDMutEx is a randomized algorithm,
% so the results may not be identical
%
% objective_function =
% 
%     18
% 
% 
% selectedGenes =
% 
%      1
%     14
%     28
%     34
%     74
%     84
%    118
%    122
%    128
%    138
%    146
%    157
% 
% 
% namesOfSelectedGenes =
% 
%   12×1 cell array
% 
%     'CDKN2B'
%     'FAM119B'
%     'RB1'
%     'ERBB2'
%     'ITGB3'
%     'TRIM2'
%     'WEE1'
%     'CHD5'
%     'MARK4'
%     'CES3'
%     'SHH'
%     'IQGAP1'
% 
% 
% solutionMetrics = 
% 
%   struct with fields:
% 
%        qObj: 18
%        lObj: 72
%        covR: 0.9286
%     excessR: 0.0769
%        mCnt: 84
%        sCnt: 84
%        gCnt: 12
%      covTot: 78
%      covOvr: 6
