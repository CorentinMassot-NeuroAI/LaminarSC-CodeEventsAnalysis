function g = filter_FR(spk,type,sampling,param1,param2)

%function g = filter_FR(spk,type,sampling,param1,param2)
%
% type='gauss': smooth signal spk with gaussian of standard deviation sigma and remove
% spurious values at the borders
%
% type='epsp': use a model of epsp temporal dynamic to filter the spiking
% activity and avoid possible problem of non-causal filtering introduced by
% the gausian filtering.
%
% Corentin University of Pittsburgh, 11/12/2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

debug=0;

spk = spk(:);

if debug
    figsig=figure;subplot(1,2,2);hold on
    plot(spk*100)
end;

switch type
    case 'gauss'
        %gaussian filtering
        %gaussian window
        
        sigma=param1*sampling/1000;
        len = 2*round(sigma*4);
        f = gausswin(len,len/sigma);
        f = f/sum(f);
        if debug
            figure(figsig);subplot(1,2,1);hold on
            plot(f)
        end
        %other definition of gaussian window
        %l=[-len/2+1:len/2];
        %sigma=sigma/2;
        %f = (1/(sigma*sqrt(2*pi)))*exp((-l.^2)/(2*sigma^2));
        
        
        
    case 'epsp'
        %epsp
        %tau1 = 2;
        %tau2 = 10/9;
        %%epspf = exp(-f(f>=0)/tau1)-exp(-f(f>=0)/tau2);
        
        t=[1:150];
        tau_g=param1;
        tau_d=param2;
        f = (1-exp(-t/tau_g)).*exp(-t/tau_d);
        f = f/sum(f);
        %centering f
        f2 = [zeros(1,floor(1*length(f))),f];
        f=f2';
        len=length(f);
        
        if debug
            figure(figsig);subplot(1,2,1);hold on;
            plot(1:length(f),f)
        end
        
end

%convolution
g = conv(spk,f);

% %taking care of borders %does not work
% c = cumsum(f);
% g(1:length(c)) = g(1:length(c))./c;
% g(end-length(c)+1:end) = g(end-length(c)+1:end) ./ flipud(c);

%taking central part of convolved signal
g = g(len/2+1:end-len/2+1);

%sampling
g=sampling * g;

if debug
    figure(figsig);subplot(1,2,2);hold on;
    plot(g)
    pause
end
