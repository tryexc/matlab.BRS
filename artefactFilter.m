function y = artefactFilter(x,p)
% For use on non resampled RR data!!
% Calculating reference values RR_ref for each original RR_i
% Check, if difference (abs) between RR_ref and RR_i in % is smaler than p
%
% Y = artefactFilter(X)
% Y = artefactFilter(X,p)
% 
% Idear from Andrea Horn, 
% PhdThesis - "Diagnostik der Herzfrequenz-
%              variabilität in der Sport-medizin - 
%              Rahmenbedingungen und methodische 
%              Grundlagen"
% Implementation by Robert Malinowski  
if nargin < 2
    p = 0.2;
end
isrowx = isrow(x);
if isrowx
    x = x(:);   % If a row, turn into column vector
end

N = size(x,1);
med = zeros(N,1); % array for the median values
xRefVal = zeros(N,1); % array for the reference values

winLength = 5;

% first step: cal. the median
        for k = (winLength+1):(N-winLength)
            dummy = zeros(1,2*winLength);
            for j = 1:winLength
                dummy(j) = x(k-j);
                dummy((2*winLength+1)-j) = x(k+j);
            end
            med(k)= median(dummy);
        end
% second step: cal. the mean of the median
        for k =(2*winLength+1):(N-2*winLength)
            dummy = zeros(1,2*winLength);
            for j = 1:winLength
                dummy(j) = med(k-j);
                dummy((2*winLength+1)-j) = med(k+j);
            end
            xRefVal(k)= mean(dummy);
% Third step: check the differences             
            if(abs(x(k) - xRefVal(k)) / x(k) > p)
                 x(k) =  xRefVal(k);
            end
        end
% return
     y = x; 
 
if isrowx
    y = y.';
end
