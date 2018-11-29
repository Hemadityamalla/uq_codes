function [weights, startid] = gpweights(maxlev)
% GPWEIGHTS   Get the 1D integration weights for the
% Gauss-Patterson grid.

% Determine size of weights vector
n = 0;
for lev = 0:maxlev
	n = n + 2^lev;
end
weights = zeros(n, 1);
startid = zeros(maxlev+1,1);
wid = 0;

for lev = 0:maxlev
	w = getweights(lev);
	nw = length(w);
	weights(wid+1:wid+nw) = w;
	startid(lev+1) = wid + 1;
	wid = wid + nw;
	if lev > 0
		weights(wid+1:wid+nw) = w(end:-1:1);
		wid = wid + nw;
	end
end

% ------------------------------------------------------------
function w = getweights(level);

switch level
 case 0
	w = 2;
 case 1
	w = 0.555555555555555555556e-0;
 case 2
	w = [0.104656226026467265194e-0 0.401397414775962222905e-0];
 case 3
	w = [0.170017196299402603390e-1 0.929271953151245376859e-1 ...
       0.171511909136391380787e-0 0.219156858401587496404e-0];
 case 4
	w = [0.254478079156187441540e-2 0.164460498543878109338e-1 ...
       0.359571033071293220968e-1 0.569795094941233574122e-1 ...
       0.768796204990035310427e-1 0.936271099812644736167e-1 ...
       0.105669893580234809744e-0 0.111956873020953456880e-0];
 case 5
	w = [0.363221481845530659694e-3 0.257904979468568827243e-2 ...
       0.611550682211724633968e-2 0.104982469096213218983e-1 ...
       0.154067504665594978021e-1 0.205942339159127111492e-1 ...
       0.258696793272147469108e-1 0.310735511116879648799e-1 ...
       0.360644327807825726401e-1 0.407155101169443189339e-1 ...
       0.449145316536321974143e-1 0.485643304066731987159e-1 ...
       0.515832539520484587768e-1 0.539054993352660639269e-1 ...
       0.554814043565593639878e-1 0.562776998312543012726e-1];
 case 6
	w = [0.505360952078625176247e-4 0.377746646326984660274e-3 ...
       0.938369848542381500794e-3 0.168114286542146990631e-2 ...
       0.256876494379402037313e-2 0.357289278351729964938e-2 ...
       0.467105037211432174741e-2 0.584344987583563950756e-2 ...
       0.707248999543355546805e-2 0.834283875396815770558e-2 ...
       0.964117772970253669530e-2 0.109557333878379016480e-1 ...
       0.122758305600827700870e-1 0.135915710097655467896e-1 ...
       0.148936416648151820348e-1 0.161732187295777199419e-1 ...
       0.174219301594641737472e-1 0.186318482561387901863e-1 ...
       0.197954950480974994880e-1 0.209058514458120238522e-1 ...
       0.219563663053178249393e-1 0.229409642293877487608e-1 ...
       0.238540521060385400804e-1 0.246905247444876769091e-1 ...
       0.254457699654647658126e-1 0.261156733767060976805e-1 ...
       0.266966229274503599062e-1 0.271855132296247918192e-1 ...
       0.275797495664818730349e-1 0.278772514766137016085e-1 ...
       0.280764557938172466068e-1 0.281763190330166021307e-1];
 otherwise
	warning('MATLAB:spinterp:unsupportedWeights',['Maximum depth (level 6) ' ... 
	  'exceeded for precomputed weights of Gauss-Patterson grid. Using all-zero ' ...
		'weights for levels > 6.']);
	w = zeros(2^(level-2),1);
end

% Divide by 2 since range is [ 0, 1] instead of [-1, 1].
w = w/2;
