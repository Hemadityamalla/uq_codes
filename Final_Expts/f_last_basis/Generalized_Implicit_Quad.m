%clear;format long;
%set(0,'DefaultAxesFontSize',16,'DefaultAxesFontWeight','bold','DefaultLineLineWidth',2,'DefaultLineMarkerSize',8);

function main(Qrule, testFn)

    Kmax = 1e2;

    
%Degree/ cardinality of the basis
    for N = [2:2:30]
        %Basis definition
        switch Qrule
        case 1
        %log(x)
        f_deets.fn = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(log(x).*x.^(k/2)); f_deets.coeffs = 1:N;
        case 2
        %sin(x)
        f_deets.fn = @(x,k) (mod(k-1,2)==0).*(x.^((k-1)/2)) + (mod(k-1,2)==1).*(sin(x).*x.^(k/2)); f_deets.coeffs = 1:N;
        case 3
         %x^k
        f_deets.fn = @(x,k) x.^(k-1); f_deets.coeffs = 1:N; 
        end
        
        mean_error = 0.0;
        Navg = 10;
        for iter = 1:Navg
            %Almost Implicit Quadrature rule
            [x,w,y] = almost_fixed_implicit_quad(N, f_deets.fn, Kmax);
            %Integrand
            a = rand(1); u = rand(1);
            [fn_quad, ~] = genz(x,a,u, testFn);
            [fn_samples, ~] = genz(y,a,u, testFn);
            exact = mean(fn_samples);
            quadrature = sum(fn_quad.*w);
            mean_error = mean_error + abs(quadrature - exact);
        end
        mean_error = mean_error/Navg;
        fname = strcat('Errors_quad',num2str(Qrule),'_fn',num2str(testFn),'.dat');
        if N == 5
            dlmwrite(fname,mean_error);
        else
            dlmwrite(fname,mean_error,'-append');
        end
    end

end
