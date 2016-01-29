[num,txt,raw] = xlsread('ctrl_FasGI7_miR217Tg_FasGI7_otros_targets.xlsx');
%[num,txt,raw] = xlsread('ctrl_FasGI7_miR217Tg_FasGI7.xls');
%[num,txt,raw] = xlsread('ctrl_DN_ctrl_FasGI7.xls');
%[num,txt,raw] = xlsread('ctrl_DN_miR217Tg_DN.xls');
[len,wid] = size(txt);
totalExpected = len-1;
totalObserved = 0;
expectedUp = 0;
expectedDn = 0;
observedUp = 0;
observedDn = 0;
indexFC = 1; % COM Fold change
%indexFC = 16; % CUFF Fold change
indexmiR = 36; % mirna scroll
%indexmiR = 33; % mirna algorithms
%DnTh = 0;
%UpTh = 0;
DnTh = -0.01;
UpTh = 0.01;
for i = 2:len,
    if (~strcmp(txt(i,indexmiR), ''))
        totalObserved = totalObserved+1;
        if (num(i-1,indexFC) <= DnTh),
            observedDn = observedDn +1;            
        end;
        if (num(i-1,indexFC) > UpTh),
            observedUp = observedUp +1;            
        end;
    end;
    if (num(i-1,indexFC) <= DnTh),
       expectedDn = expectedDn +1;            
    end;
    if (num(i-1,indexFC) > UpTh),
       expectedUp = expectedUp +1;            
    end;
end

[pValUpReg, pValDnReg, pValDn_Up] = cossgsea(totalExpected, expectedUp,expectedDn, totalObserved, observedUp, observedDn, 'mir217', 1);
disp(['miR217  '  'pValUpReg = ' num2str(pValUpReg) ' pValDnReg = ' num2str(pValDnReg) ' pValDn-Up = ' num2str(pValDn_Up)]);



% CossGsea2

% For the Upregulateds
[realDifference] = cossgsea2(totalExpected, expectedUp, totalObserved, observedUp);
[realDifferenceDn] = cossgsea2(totalExpected, expectedDn, totalObserved, observedDn);
diffExp_Obs = zeros(1,10000);
for rndI = 1:10000,
    % Now a random vector
    myRandIdx = unidrnd(totalExpected,1,totalObserved);

    observedUp = 0;
    observedDn = 0;

    for i = 1:totalObserved,
        if (num(myRandIdx(1,i),indexFC) <= DnTh),
            observedDn = observedDn +1;            
        end;
        if (num(myRandIdx(1,i),indexFC) > UpTh),
            observedUp = observedUp +1;            
        end;
    end
    [diffExp_Obs(rndI)] = cossgsea2(totalExpected, expectedUp, totalObserved, observedUp);
end; 
pvalue = sum(diffExp_Obs < realDifference)/length(diffExp_Obs)

% For the Downregulateds.
diffExp_Obs = zeros(1,10000);
for rndI = 1:10000,
    % Now a random vector
    myRandIdx = unidrnd(totalExpected,1,totalObserved);

    observedUp = 0;
    observedDn = 0;

    for i = 1:totalObserved,
        if (num(myRandIdx(1,i),indexFC) <= DnTh),
            observedDn = observedDn +1;            
        end;
        if (num(myRandIdx(1,i),indexFC) > UpTh),
            observedUp = observedUp +1;            
        end;
    end
    [diffExp_Obs(rndI)] = cossgsea2(totalExpected, expectedDn, totalObserved, observedDn);
end; 
pvalue = sum(diffExp_Obs < realDifferenceDn)/length(diffExp_Obs)
