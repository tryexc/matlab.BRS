function brs = brsX(ibi,bp,tTotsec,plt)
%
%
% September 2016
% Implementation by Robert Malinowski

try
% ########### Atefact-Filter
%bpFil = artefactFilter(bp,0.3);
%ibiFil = artefactFilter(ibi,0.3);
bpFil = bp;
ibiFil = ibi;
Fs = 1;   % Samplingrate for resampling in Hz

% ########### Resampling
bpFilRes = resample(bpFil,tTotsec,Fs,1,1,'spline');
ibiFilRes = resample(ibiFil,tTotsec,Fs,1,1,'spline');

k=1;
xbrs = 0;
 for i=1:(length(ibiFilRes)-10-5)
     r=0;
     p=1;
     t=0;
     currMaxCorr = struct('ibi',[],'bp',[],'r',r,'p',p,'tau',t);
     for j=1:6
         x=i+j-1;
     
        bpWin = bpFilRes(i:(i+10));
        rrWin  = ibiFilRes(x:(x+10));
        [r,p]  = corrcoef(bpWin,rrWin);

        if currMaxCorr.r <= r(2,1) && p(2,1) <= 0.01
            currMaxCorr.ibi = rrWin;
            currMaxCorr.bp = bpWin;
            currMaxCorr.r = r(2,1);
            currMaxCorr.p = p(2,1);
            currMaxCorr.tau = j-1;
        end
     
     end 
       if currMaxCorr.r > 0
           
           [x, y, yFit, m, r2] = linReg(currMaxCorr.bp,currMaxCorr.ibi);
    
        if r2 >= 0.6 % R^2 Value 
            regAll(k) =  struct('X',x,'Y', y,'Tau',currMaxCorr.tau,'Fit', yFit,'m',m, 'R2', r2);
            k=k+1;
            xbrs = xbrs + m(2);
        end
             
           
       end
 end
  
 % ######### Return 
 brs= struct('AllVal', regAll ,'xBRS', xbrs/(k-1), 'xBRSN',k-1,'isOK',1);
 
 % ######### Create a figure
 if plt
     plotBRSX(regAll)
 end

catch
     brs= struct('AllVal', NaN ,'xBRS', NaN, 'xBRSN', NaN,'isOK',0);
end
 
end

% ######### Simple linear - Fit 
 function [x, y, yy, m, r2] = linReg(x,y)
    
    X = [ones(length(x),1) x];
    b = X\y;
    
    yy2 = X*b;
    
    Rsq2 = 1 - sum((y - yy2).^2)/sum((y - mean(y)).^2);

        yy=yy2;
        m=b;
        r2=Rsq2;
    
end

 