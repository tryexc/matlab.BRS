function brs = brsSequence(ibi,bp,lag,p)
% 
%
% 
% September 2016
% Implementation by Robert Malinowski  


% ########### Atefact-Filter
%bpFil = artefactFilter(bp,0.1); 
%ibiFil = artefactFilter(ibi,0.1);
bpFil = bp;
ibiFil = ibi;


% ########### Find up and down Seq.
[Xup,Yup] = findUPSeq(bpFil,ibiFil,lag);
[Xdown,Ydown] = findDOWNSeq(bpFil,ibiFil,lag);

% ########### Calculate lin. regression - parameter
BRSupN = 0;
BRSup = 0;
BpUp = [];
 if not(isempty(Xup))
    BpUp = getRegrVal(Xup,Yup);
    % ########## Calculateing the BRS-Values
    % ########## UP-Seq.
    if isstruct(BpUp)
    BRSupN = length(BpUp);
    BRSup = 0;
        for i=1:BRSupN
            BRSup = BRSup + BpUp(i).m(2);
        end
    BRSup = BRSup / BRSupN;
    end
 end
 
BRSdownN = 0;
BRSdown = 0;
BpDown = [];
 if not(isempty(Xdown))
    BpDown = getRegrVal(Xdown,Ydown);
    % ########## Calculateing the BRS-Values
    % ########## Down-Seq.
    if isstruct(BpDown)
    BRSdownN = length(BpDown);
    BRSdown = 0;
        for i=1:BRSdownN
            BRSdown = BRSdown + BpDown(i).m(2);
           
        end
    BRSdown = BRSdown / BRSdownN;
    end
 end

% ########## Weighted Value.
M = BRSupN + BRSdownN;
wBRS = (BRSupN/M)*BRSup + (BRSdownN/M)*BRSdown;

% ########## Return

brs = struct('BPvalUp',BpUp,'BPvalDown',BpDown,...
             'BRSup',BRSup,'BRSupN',BRSupN,...
             'BRSdown',BRSdown,'BRSdownN',BRSdownN,...
             'wBRS',wBRS,'lag',lag);
if p     
    fig = plotBRSseq(brs);
end





end

% ########## Return all Fits
function RegOut = getRegrVal(X,Y)
 RegOut =  struct('X',[],'Y', [],'Fit', [],'m', [], 'R2', []);
k=1;
    for i=1:size(X,2)
       
        xIn = X(:,i);
        xIn(xIn==0)=[];
        
        yIn = Y(:,i);
        yIn(yIn==0)=[];
        
        [x, y, yFit, m, r2] = linReg(xIn,yIn);
    
        if r2 >= 0.85 % R^2 Value 
        RegOut(k) =  struct('X',x,'Y', y,'Fit', yFit,'m',m, 'R2', r2);
        k=k+1;
        end
    end
    
    if k==1
         RegOut = -1;
    end
    
   
end

% ######### Simple linear - Fit 
 function [x, y, yy, m, r2] = linReg(x,y)

%     b1 = x\y;
%     yy1 = b1*x;
    
    X = [ones(length(x),1) x];
    b = X\y;
    
    yy2 = X*b;
    
%     Rsq1 = 1 - sum((y - yy1).^2)/sum((y - mean(y)).^2);
    Rsq2 = 1 - sum((y - yy2).^2)/sum((y - mean(y)).^2);
    
%     if Rsq1 > Rsq2
%        yy=yy1;
%         m=b1
%        r2=Rsq1;
%     else
%        yy=yy2;
%         m=b;
%        r2=Rsq2;
%     end
    
        yy=yy2;
        m=b;
        r2=Rsq2;
    
end

% ######### Find BP-Up Sequences
function [x,y] = findUPSeq(bp,ibi,lag)
% ######### Define Criteria
dBPmin = 1.00;     % (mmHg)
dIBImin = 6.00;    % (ms)
seqBmin = 3;    % (beats)
    seqX = zeros(10,1000);
    seqY = zeros(10,1000);
    k = 1;
    s = 0; %
    for i=2:length(bp)
            pBP = zeros(10,1);
            pIBI = zeros(10,1);
            l=1;
        
        while i+s+lag<=length(bp) && bp(i+s) > bp(i+s-1) && ...
                ibi(i+s+lag) > ibi(i+s+lag-1) && ...
                abs(bp(i+s)-bp(i+s-1)) >= dBPmin && ...
                abs(ibi(i+s+lag)-ibi(i+s+lag-1)) >= dIBImin
            pBP(l) =  bp(i+s-1);
            pIBI(l) = ibi(i+s+lag-1);
            l=l+1;
            pBP(l) =  bp(i+s);
            pIBI(l) = ibi(i+s+lag);
            s=s+1; % necessary, to avoid generating subsequences
        end

       if l >= seqBmin && l <= 10
       seqX(:,k) = pBP;
       seqY(:,k) = pIBI;
       k=k+1;
       end
    end
% ####### A bit magic ########
% ####### this "thing" deletes all zero-cols and -rows
    seqX(~any(seqX,2),:)=[];
    seqX(:,~any(seqX))=[];
    
    seqY(~any(seqY,2),:)=[];
    seqY(:,~any(seqY))=[];
    
% ####### return   
   x = seqX;
   y = seqY;
end

% ########## Find BP-Down Sequences
% analogous to findUPSeq()
function [x,y] = findDOWNSeq(bp,ibi,lag)
% ######### Define Criteria
dBPmin = 1.00;     % (mmHg)
dIBImin = 6.00;    % (ms)
seqBmin = 3;    % (beats)
    seqX = zeros(10,1000);
    seqY = zeros(10,1000);
    k = 1;
    s = 0;
    for i=2:length(bp)
            pBP = zeros(10,1);
            pIBI = zeros(10,1);
            l=1;
        while i+s+lag<=length(bp) && bp(i+s) < bp(i+s-1) && ...
                ibi(i+s+lag) < ibi(i+s+lag-1) && ...
                abs(bp(i+s)-bp(i+s-1)) >= dBPmin && ...
                abs(ibi(i+s+lag)-ibi(i+s+lag-1)) >= dIBImin
            pBP(l) =  bp(i+s-1);
            pIBI(l) = ibi(i+s+lag-1);
            l=l+1;
            pBP(l) =  bp(i+s);
            pIBI(l) = ibi(i+s+lag);
            s=s+1; % necessary, to avoid generating subsequences
        end
        
       if l >= seqBmin && l <= 10
       seqX(:,k) = pBP;
       seqY(:,k) = pIBI;
       k=k+1;
       end
    end
% ####### A bit magic ########
% ####### this "thing" deletes all zero - cols and rows
    seqX(~any(seqX,2),:)=[];
    seqX(:,~any(seqX))=[];
    
    seqY(~any(seqY,2),:)=[];
    seqY(:,~any(seqY))=[];
    
% ####### return   
   x = seqX;
   y = seqY;
end


