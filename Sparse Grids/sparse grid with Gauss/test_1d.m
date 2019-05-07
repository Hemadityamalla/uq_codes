%This script was written to compute gPC approximations in 1d/2d and verify
%my own functions 'gpc_coeffs' and 'gpc_polyval'


clc
close all
clearvars

% add path to UQLab
addpath(genpath('~/Documents/UQLabCore_Rel1/'));


%% start uqlab
uqlab

%% set up model
Model.mFile = 'test1D';               % name of Matlab file representing the model
myModel     = uq_createModel(Model);  % create and add the model to UQLab   


%% set up inputs

% Description of input marginals
% for ii = 1 : 10
%     Input.Marginals(ii).Type = 'Uniform'; 
%     Input.Marginals(ii).Parameters = [0 1];
% end
Input.Marginals.Type = 'Uniform';
Input.Marginals.Parameters = [-1,1];
myInput = uq_createInput(Input);

% show information about input:
uq_print(myInput);
% graphical display input:
% uq_display(myInput);


%% set up Polynomial Chaos Expansion as Metamodel
metamodelQuad.FullModel = myModel;  
metamodelQuad.Input     = myInput;
metamodelQuad.Type      = 'Metamodel';
metamodelQuad.MetaType  = 'PCE';
metamodelQuad.Method          = 'quadrature'; % using quadrature to compute coefficients
metamodelQuad.Degree = 1;

myPCE_Quad           = uq_createModel(metamodelQuad);
uq_print(myPCE_Quad)
myPCE_Quad.PCE.Coefficients

%Evaluating the gPC at some samples
xtest = linspace(-1,1,500);
Y = uq_evalModel(xtest');
plot(xtest,Y,'b.');
hold on;
plot(xtest,test1D(xtest),'r-');

%% some analysis of PCE:

% NsamplesPCE = myPCE_Quad.ExpDesign.NSamples
% mean_PCE    = myPCE_Quad.PCE.Moments.Mean
% std_PCE     = sqrt(myPCE_Quad.PCE.Moments.Var)
