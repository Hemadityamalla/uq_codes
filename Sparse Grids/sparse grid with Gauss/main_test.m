%Script to construct gPC approximation of an r.v. Y=f(X) using Galerkin
%projection.

clc
close all
clearvars

% add path to UQLab (this needs to be modified accordingly)
addpath(genpath('~/Documents/UQLabCore_Rel1/'));


%% start uqlab
uqlab

%% set up model
Model.mFile = 'test1D';               % name of Matlab file representing the model
myModel     = uq_createModel(Model);  % create and add the model to UQLab   


%% set up inputs

% Description of input marginals (uniform distribution in 1D, in the interval [-1, 1])
Input.Marginals.Type = 'Uniform';
Input.Marginals.Parameters = [-1,1];
myInput = uq_createInput(Input);

% show information about input:
uq_print(myInput);

%% set up Polynomial Chaos Expansion as Metamodel
metamodelQuad.FullModel = myModel;  
metamodelQuad.Input     = myInput;
metamodelQuad.Type      = 'Metamodel';
metamodelQuad.MetaType  = 'PCE';
metamodelQuad.Method    = 'quadrature'; % using quadrature to compute coefficients
maxPolyDegree = 2;
metamodelQuad.Degree = maxPolyDegree;

%I am using the '..Quadrature.Level' command below with the understanding
%that it will give me a quadrature rule of degree higher than the default
%setting of '(maxPolyDegree + 1)'. But this gives me an error (matrix dimension exceeded). 
metamodelQuad.Quadrature.Level = 5;

myPCE_Quad           = uq_createModel(metamodelQuad);
uq_print(myPCE_Quad)
myPCE_Quad.PCE.Coefficients %prints the coefficients of the gPC expansion

%Evaluating the gPC at some samples
interval = Input.Marginals.Parameters;
xtest = linspace(interval(1),interval(2),500);
Y = uq_evalModel(xtest');
plot(xtest,Y,'b.');
hold on;
plot(xtest,test1D(xtest),'r-');
