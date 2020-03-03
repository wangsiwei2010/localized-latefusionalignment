clear
clc
warning off;

path = 'D:\datasets\MultipleKernelDatasets\';
addpath(genpath(path));

DataName = cell(11, 1);

DataName{1} = 'caltech101_nTrain5_48';
DataName{2} = 'caltech101_nTrain10_48';
DataName{3} = 'caltech101_nTrain15_48';
DataName{4} = 'caltech101_nTrain20_48';
DataName{5} = 'caltech101_nTrain25_48';
DataName{6} = 'caltech101_nTrain30_48';
DataName{7} = 'CCV';
DataName{8} = 'flower17';
DataName{9} = 'YALE';
DataName{10} = 'plant';
DataName{11} = 'mfeat';
DataName{12} = 'AR10P';


%%addpath('C:\Program Files\Mosek\8\toolbox\r2014a');
for i = 1:12 
    dataName = DataName{i};
    %%% flower17; flower102; proteinFold,caltech101_mit,UCI_DIGIT,ccv
    %% %% washington; wisconsin; texas; cornell
    %% caltech101_nTrain5_48
    %% proteinFold
    load([path,dataName,'_Kmatrix'],'KH','Y');
    % load([path,'datasets\',dataName,'_Kmatrix'],'KH','Y');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numclass = length(unique(Y));
    numker = size(KH,3);
    num = size(KH,1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    KH = kcenter(KH);
    KH = knorm(KH);
    
    HP = zeros(num,numclass,numker);
    
    qnorm = 2;
    
    opt.disp = 0;
    
    for p=1:numker % m - kernels
        KH(:,:,p) = (KH(:,:,p)+KH(:,:,p)')/2;
        [Hp, ~] = eigs(KH(:,:,p), numclass, 'la', opt);
        HP(:,:,p) = Hp;
    end
    
    gamma0 = ones(numker,1)/numker;
    avgKer  = mycombFun(KH,gamma0);
    
    [H0,~] = eigs(avgKer, numclass, 'la', opt);   
    
    
        %%---The Proposed2 time---%%
    tic
    %%lambdaset9 = 2.^[-15:1:15];
    %%tauset9 = [0.1:0.1:1];
    lambdaset9 = 2.^3;
    tauset9 = 0.3;
    accval9 = zeros(length(lambdaset9),length(tauset9));
    nmival9 = zeros(length(lambdaset9),length(tauset9));
    purval9 = zeros(length(lambdaset9),length(tauset9));
    
    for it =1:length(tauset9)
        neibournum = round(tauset9(it)*num);
        NS = genarateNeighborhood(avgKer,neibournum);
        A = calculateA(NS,neibournum);
        localavg = calculatelocalH(H0,A);
        localH = zeros(num,numclass,numker);
        for i =1:numker
            index = genarateNeighborhood(KH(:,:,p),neibournum);
            s = calculateA(index,neibournum);
            localH(:,:,i) = calculatelocalH(HP(:,:,i),s);
        end
        for ij = 1:length(lambdaset9)
            tic;
            [H_normalized9,WP9,gamma9,obj9] = multikernellocalLatefusionAlignmentclustering(localH,numclass,...
                lambdaset9(ij),Y,localavg);
            res9 = myNMIACC(H_normalized9,Y,numclass);
            accval9(ij,it) = res9(1);
            nmival9(ij,it)= res9(2);
            purval9(ij,it) = res9(3);
            t = toc;
            time = t;
        end
    end
    res(:,9) = [max(max(accval9)); max(max(nmival9));max(max(purval9))];
    
    save(['.\Result\2031_Time' , dataName , '.mat'],'time')
    %%save(['.\Result\2031_Time' , dataName , '.mat'], 'res','accval9','nmival9','purval9')
    mailme(['LMVC-LFA Time ', dataName],'L-MVC-LFA complete');
end


