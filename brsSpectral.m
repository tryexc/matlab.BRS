function brs = brsSpectral(ibi,bp,tTotsec,plt)
% 
%
% 
% September 2016
% Implementation by Robert Malinowski  

%try

% ########### Atefact-Filter
%bpFil = artefactFilter(bp,0.1); 
%ibiFil = artefactFilter(ibi,0.1);

bpFil = bp; 
ibiFil = ibi;

Fs = 2;         % Samplingrate in Hz
L = 512;        % Segment-Window-Size (Points)
nfft = 2048;    % FFT Window

% ########### Resampling
bpFilRes = resample(bpFil,tTotsec,Fs,1,1,'spline');
ibiFilRes = resample(ibiFil,tTotsec,Fs,1,1,'spline');

while length(ibiFilRes) < L
    L = 2^nextpow2(L/2);
end

% ########### Detrending Data
bpFilResDetr = detrend(bpFilRes);
ibiFilResDetr = detrend(ibiFilRes);

% ########### PSD
Pxx = pwelch(bpFilResDetr,L,[],nfft,[],'psd');
Pyy = pwelch(ibiFilResDetr,L,[],nfft,[],'psd');

% ########### CSPD 
Pxy = cpsd(bpFilResDetr,ibiFilResDetr,L,[],nfft);

% ########### Coherence
Cxy = mscohere(bpFilResDetr,ibiFilResDetr,128,[],nfft);

% ########### Transfer-Function estimation
Txy = tfestimate(bpFilResDetr,ibiFilResDetr,L,[],nfft);
AngleXY = angle(Txy)/(2*pi)*360;

% ########### Alpha-Index
AlphaXY = zeros(length(Pyy));
for i=1:length(Pyy)
AlphaXY(i) = sqrt(Pyy(i)/Pxx(i));
end


% ########### Frequency-Scale
F = 0:Fs/nfft:0.5; % Max Frequency of interest is 0.5 Hz

% ########### Index-thresholds LF-Band
LFmin = find(abs(F-0.04) < 0.001, 1, 'last' );
LFmax = find(abs(F-0.15) < 0.001, 1 );

alphaLF = 0;
brsLF=0;

nLF = 0;

% ########### LF
% ########### Sum all Values at LF-Band with |K|^2 >= 0.5
areaLF = zeros(LFmax,1);
for i=LFmin:LFmax 
     if Cxy(i)>=0.5 && AngleXY(i)<0
     alphaLF = alphaLF + AlphaXY(i); % for alpha-index
     brsLF = brsLF + abs(Txy(i));      % BRS from Transferfct
     nLF = nLF + 1;
     areaLF(nLF) =  i;
     end
end
 areaLF(areaLF==0)=[];
 
% ####### BRS (LF) ###########
% ############################
alphaLF = alphaLF/nLF      ;%#
brsLF = brsLF/nLF;          %#
%#############################

% ########### Index-thresholds HF-Band
HFmin = find(abs(F-0.15) < 0.001, 1, 'last' );
HFmax = find(abs(F-0.4) < 0.001, 1 );

alphaHF = 0;
brsHF = 0;

nHF = 0;

% ########### HF
% ########### Sum all Values at HF-Band with |K|^2 >= 0.5
areaHF = zeros(HFmax,1);
for i=HFmin:HFmax
     if Cxy(i) >= 0.5 && AngleXY(i)<0
     alphaHF = alphaHF + AlphaXY(i); % for alpha-index
     brsHF = brsHF + abs(Txy(i)); % BRS from Transferfct
     nHF = nHF + 1;
     areaHF(nHF) =  i;
     end
 end
areaHF(areaHF==0)=[];

% ####### BRS (HF) ###########
% ############################
alphaHF = alphaHF/nHF;      %#
brsHF = brsHF/nHF;          %#
%#############################

%######## return #############

brs = struct('alphaLF',alphaLF,'alphaHF',alphaHF,'brsLF',brsLF,'brsHF',brsHF,...
    'F',F,'Pxx',Pxx,'Pyy',Pyy,'Pxy',Pxy,'Cxy',Cxy,'AngleXY',AngleXY,'Txy',Txy,...
    'AlphaXY',AlphaXY,'areaLF',areaLF,'areaHF',areaHF,'LFmin',LFmin,'HFmax',HFmax);


if plt
plotBRSspec(F,Pxx,Pyy,Pxy,Cxy,AngleXY,Txy,AlphaXY,areaLF,areaHF,LFmin,HFmax);
end

% catch
%     brs = struct('alphaLF',NaN,'alphaHF',NaN,'brsLF',NaN,'brsHF',NaN);
% end









