close all;clear;clc;format long;

%points used for piecewise linear interpolation
Kmax = 2.5e2;


for testFn = 1:6

    lmax = 10;
    D = 2; N = 5;
    numavg = 25;
    efiqavg = zeros(lmax,1)'; egiqavg = zeros(lmax,1)';
    for iavg = 1:numavg
        error_giq = [];
        error_fiq = [];
        samples = rand(Kmax,1);
        %Generating nested QR sequence 
        a = rand(1); u = rand(1); 
        f = genz_fns(0:0.01:1, a, u, testFn);
        [QR] = cappedaccuracy_quad_nested([D,N], lmax,f, samples);
        for ii=1:lmax
            %Generating the FIQR (Gives a QR of degree similar to the ii th QR of the nested sequence)
            [xf,wf] = fixed_implict_quad((N-2)+ii, samples);
            error_fiq(end+1) = abs(sum(f(xf).*wf) - mean(f(samples)));
            xa = QR{ii}.nodes; wa = QR{ii}.weights;
            error_giq(end+1) = abs(sum(f(xa).*wa) - mean(f(samples)));
        end
        efiqavg = efiqavg + error_fiq; egiqavg = egiqavg + error_giq;
    end
    figure(testFn);
    semilogy(1:lmax,egiqavg/numavg,'bo-',1:lmax,efiqavg/numavg,'r+-');
    xlabel('Number of nodes/levels');ylabel('Error');
    legend('adaptive','fixed');
end