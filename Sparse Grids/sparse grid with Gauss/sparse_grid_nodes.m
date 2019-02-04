clc
close all
clearvars

% add path to UQLab
addpath(genpath('~/Documents/UQLabCore_Rel1/'));


%% start uqlab
uqlab

%% set up model
Model.mFile = 'test2D';               % name of Matlab file representing the model
myModel     = uq_createModel(Model);  % create and add the model to UQLab   


%% set up inputs

% Description of input marginals
for ii = 1 : 10
    Input.Marginals(ii).Type = 'Uniform'; 
    Input.Marginals(ii).Parameters = [-1 1];
end
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
metamodelQuad.Method          = 'Quadrature'; % quadrature to compute coefficients
metamodelQuad.Quadrature.Type = 'Smolyak'; % tensor grid
metamodelQuad.Degree = 5;

myPCE_Quad           = uq_createModel(metamodelQuad);

figure
plot(myPCE_Quad.ExpDesign.X(:,1),myPCE_Quad.ExpDesign.X(:,2),'s')
hold on;
%% some analysis of PCE:

NsamplesPCE = myPCE_Quad.ExpDesign.NSamples
mean_PCE    = myPCE_Quad.PCE.Moments.Mean
std_PCE     = sqrt(myPCE_Quad.PCE.Moments.Var)
