close all;clear;format long;
set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

%points used for piecewise linear interpolation
Kmax = 5e2;
rng(1,'twister'); %Seeding for reproducibility.

xpos = 500;ypos = 500; width = 1000; height = 800;


for testFn = 1:6

    lmax = 8;
    D = 4; N = 5;
    numavg = 50;
    
    efiqavg = zeros(lmax,1)'; egiqavg = zeros(lmax,1)';
    for iavg = 1:numavg
        error_giq = [];
        error_fiq = [];
        nnodes = [];
        samples = rand(Kmax,1);
        %Generating nested QR sequence 
        a = rand(1); u = rand(1); 
        f = genz_fns(0:0.01:1, a, u, testFn);
        [QR] = cappedaccuracy_quad_nested([D,N], lmax,f, samples);
        for ii=1:lmax
            xa = QR{ii}.nodes; wa = QR{ii}.weights;
            nnodes(end+1) = length(xa);
            %Generating the FIQR (Gives a QR of degree similar to the ii th QR of the nested sequence)
            [xf,wf] = fixed_implict_quad(N-2 + ii, samples);
            error_fiq(end+1) = abs(sum(f(xf).*wf) - mean(f(samples)));
            
            error_giq(end+1) = abs(sum(f(xa).*wa) - mean(f(samples)));
        end
        efiqavg = efiqavg + error_fiq; egiqavg = egiqavg + error_giq;
    end
    figure(testFn);
    loglog(nnodes,egiqavg/numavg,'bo-',nnodes,efiqavg/numavg,'r+-');
    xlabel('Number of nodes/levels');ylabel('Error');
    legend('adaptive','fixed');
    grid on;set(gcf,'Position',[xpos ypos width height]); box on;
end