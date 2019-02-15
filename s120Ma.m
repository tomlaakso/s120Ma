
%%% Model for 'Volcanic controls on seawater sulfate over the
%%% past 120 million years' by T. Laakso, A. Waldeck, D. Johnston
%%% and F. Macdonald.

%%% Written and tested in Matlab version R2016a by T. Laakso.

%%% Input: In order to run this code, the user must supply the
%%%        d34S barite record from:
%%%          A. Paytan, M. Kastner, D. Campbell, M. Thiemens, 
%%%          Seawater Sulfur Isotope Fluctuations in the Cretaceous. 
%%%          Science 304, 1663-1665 (2004).
%%%          DOI: 10.1126/science.1095258
%%%         This data should be configured as a comma-separated file
%%%         giving the age of the sample in Ma in the first column,
%%%         and the delta34S of the sample in permil in the second
%%%         column. The full file name (including path) must
%%%         be saved as a string in the 'nomen' variable, below:
            nomen = '';
%%% Output: d0: a vector of delta34S values of seawater sulfate
%%%             over the simulation (permil)
%%%         M0: a vector of the sulfate content of seawater over
%%%             the simulation (mol S)
%%%         t0: a vector giving the ages for the elements of
%%%             d0 and M0 in millions of years before present.

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% User specifications %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%These are really used for exploration of the parameter space using
%index terms i1, i2, etc. To run a single integration, assign the SET
%terms below single values and fix the index looping terms {i} to 1.
%{
iset1 = ;   %Amplitude of 100 Ma event (multiplicative)
iset2 = ;   %Amplitude of 50 Ma event  (multiplicative)
iset3 = ;   %d34S of 100 Ma input (additive, baseline +2)
iset4 = ;   %d34S of 50 Ma input  (additive, bseline +2)
iset5 = ;   %Peak of 100 Ma event (additive, baseline 101 Ma)
iset6 = ;   %Peak of 50 Ma event  (additive, baseline 53 Ma)
iset7 = ;   %st.dev. of 100 Ma event (mult., baseline 5 My)
iset8 = ;   %st.dev. of 100 Ma event (mult., baseline 5 My)
iset9 = 0;  %Unused
i1 = 1;i2 = 1;i3 = 1;i4 = 1;i5 = 1;i6 = 1;i7 = 1;i8 = 1;i9 = 1;
inclmod = 0; %Include modern?
%}

%Current best fit
inclmod = 0;
lax = 10;
iset1 = linspace(.5,1.5,lax);
iset2 = linspace(.5,1.5,lax);
 iset3 = -1:2;
 iset4 = -1:2;
iset5 = linspace(-.5,1.5,lax); 
iset6 = linspace(-.5,1.5,lax);
iset7 = linspace(.5,1.5,lax);
iset8 = 1;
iset9 = 0;
i1=6;i2=6;i3=2;i4=2;i5=6;i6=5;i7=4;i8=1;i9=1;

warning('off','MATLAB:singularMatrix')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Smooth data for fit evaluation or presentation %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Using loess
tcut = 120;
a = dlmread(nomen);
[tAd,I] = sort(a(:,1));
dAd = a(:,2);dAd = dAd(I);
dAd = dAd(tAd<tcut);tAd = tAd(tAd<tcut);
w = 15; %smoothing window width in Ma (Paytan)
dts = 1; %smooth stepping
dt = 1;
ts = 0:dt:120;
ds = [];
for tau = ts
    z = tAd<(tau+w/2) & tAd>(tau-w/2);
    wt = (1-abs((tau-tAd)/w).^3).^3; %tricube
    wt(~z) = 0;
    wt = diag(wt);
    x = [ones(length(tAd),1) tAd tAd.^2];
    y = (dAd);
    bet = inv(x'*wt*x)*x'*wt*y;
    ds = [ds bet(1)+bet(2)*tau+bet(3)*tau^2];
end
ts = ts(end:-1:1); %smoothed time vector
ds = ds(end:-1:1); %smoothed d34S

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Basic parameter set-up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dA = 16.75;                          %initial steady-state d34S       
bd = (25.2+iset9(i9)-dA)/120/10^6;   %slope of steady-state d34S drift
Mdiv = 5;                 
MA = 28/Mdiv*10^-3*10^21;  %initial steady state [S]
bM = (Mdiv-1)*MA/120/10^6; %slope of steady-state [S] drift
Fw = MA*Mdiv/10/10^6;      %S weathering input
dw = 0;                    %S weathering input d34S

%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Perturbation set-up %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
megaI = 3;                 %S perturbation amplitude relative to initial [S]
ampM = MA*megaI*2;         %baseline perturbation amplitude
pertdelta = 2;             %d34S of S perturbation flux
%Event-specific tuning of perturbation amplitude & isotopic composition
ampM = [.5*ampM*iset1(i1) 0.6*ampM*iset2(i2)];
pertdelta = [pertdelta+iset3(i3) pertdelta+iset4(i4)];
%Event-specific onset
pt = [(19+iset5(i5))*10^6 (67+iset6(i6))*10^6];
%Event-specific duration
pl = 5*10^6;pl = [pl*iset7(i7) pl*iset8(i8)];

%Perturbation vectors:
 %pM is the rate of perturbed S input at each time step
 %pertd is isotopic composition of that input
dt = 10^5;
t = 0:dt:1.2*10^8;
pM = 0*t;
pertd = 0*t;
for i = 1:length(pt)
 Mbit = ampM(i)/sqrt(2*pi*pl(i)^2)*exp(-(t-pt(i)).^2/(2*pl(i)^2));
 Mbit(t>pt(i))=0;         %truncate in time
 Mbit((t-pt(i))<3*pl(i)); %truncate in time
 dbit = repmat(pertdelta(i),1,length(pertd));
 dbit(Mbit==0) = 0;
 dbit(pertd~=0) = pertd(pertd~=0);
 pM = pM+Mbit;
 pertd = dbit;
end

%%%%%%%%%%%%%%%%%
%%% Parameter %%%
%%%%%%%%%%%%%%%%%

%These terms are used later to evaluate goodness of fit
keys = [0 0 0 0 0 0 0 0 0 10^10];
zkey = 10^10;

%If you wish to loop over the perturbation parameters, initialize here, e.g
%for i1 = 1:length(iset1)
%for i2 = 1:length(iset2)
%...etc

%%%%%%%%%%%%%%%%%%%%%%%%
%%% Initialize ocean %%%
%%%%%%%%%%%%%%%%%%%%%%%%
M0 = MA+bM*t;  %Steady-state [S]  at each time step
d0 = dA+bd*t;  %Steady-state d34S at each time step
Mnow = M0(1);  %Current [S]
dnow = d0(1);  %Current d34S
%Note: it is helpful to spin-up the model, as it is never truly
%at equilibrium due to the drift in the steady-state values; it
%goes through an initial oscillation as it relaxes toward smooth
%behavior
spint = -10*10^6:dt:0; %spin-up time vector
M0spin = MA+bM*spint;  %steady-state [S] at each spin-up time step
d0spin = dA+bd*spint;  %steady-state d34S at each spin-up time step
%Spin-up (Euler integration)
for i = 1:length(spint);
    dM = Fw*(1-Mnow/M0spin(i));
    dd = (dw-dnow)*Fw/Mnow-(dw-d0spin(i))*Fw/M0spin(i);
    Mnow = Mnow+dM*dt;
    dnow = dnow+dd*dt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Primary integration (Eulerian) %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = Mnow;
d = dnow;
for i = 1:length(t)
    dM = Fw*(1-Mnow/M0(i));
    dd = (dw-dnow)*Fw/Mnow-(dw-d0(i))*Fw/M0(i);
    pd = (pertd(i)-dnow)*pM(i)/Mnow;
    Mnow = Mnow+dM*dt+pM(i)*dt;
    dnow = dnow+dd*dt+pd(1)*dt;
    M = [M Mnow];
    d = [d dnow];
end
M = M(1:end-1);
d = d(1:end-1);
%Final time, d34S, and [S] vectors of integration
t0 = 120-t/10^6;
d0 = d;
M0 = M;

%%%%%%%%%%%%%%%%%%%%%%
%%% Fit evaluation %%%
%%%%%%%%%%%%%%%%%%%%%%
%Trimming: fit only for certain ages
d0z = d0(t0<tcut);t0z = t0(t0<tcut);
dsz = ds(ts<tcut);tsz = ts(ts<tcut);
%Metric is root sum of square residuals relative to full data.
 ykey = 0;
 for i = 1:length(dAd)
     mk = 1:length(t0);
     mk = mk(abs(tAd(i)-t0z)==min(abs(tAd(i)-t0z)));
     ykey = ykey+mean((dAd(i)-d0z(mk)).^2);
 end
 keynorm = length(dAd);
 ykey = sqrt(ykey);
%For parameter exploration, preserve the current best fit parameters
%and plot current best run
if ykey<zkey
    zkey = ykey;
    ikey = [i1 i2 i3 i4 i5 i6 i7 i8 i9];
    try
        close(1)
    end
    plot(tAd,dAd,'.');hold on
    plot(ts,ds)
    plot(t0,d0);
    legend('data','loess','simulation','Location','Northeast')
    xlabel('time (Ma)')
    ylabel(['\delta^{34}S (' char(8240) ')'])
    drawnow
end
%Save top 10 parameter sets
if (ykey/sqrt(keynorm))<max(keys(:,end));
    keys = [keys;[i1 i2 i3 i4 i5 i6 i7 i8 i9 ykey/sqrt(keynorm)]];
    [junk,I] = sort(keys(:,end));
    keys = keys(I,:);
    keycut = min([10 size(keys,1)]);
    keys = keys(1:keycut,:);
end
