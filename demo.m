clear
path = '*';
addpath(genpath(path));
dataName = 'caltech101_mit';
for iter = 1:10

    load(['datasets\',dataName,'_Kmatrix'],'KH','Y');

    options.seuildiffsigma=1e-5;        % stopping criterion for weight variation
    options.goldensearch_deltmax=1e-1; % initial precision of golden section search
    options.numericalprecision=1e-16;   % numerical precision weights below this value are set to zero
    options.firstbasevariable='first'; % tie breaking method for choosing the base variable in the reduced gradient method
    options.nbitermax=500;             % maximal number of iteration
    options.seuil=0;                   % forcing to zero weights lower than this
    options.seuilitermax=10;           % value, for iterations lower than this one    
    options.miniter=0;                 % minimal number of iterations
    options.threshold = 1e-4;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    numclass = length(unique(Y));
    Y(Y<1) = numclass;
    numker = size(KH,3);
    num = size(KH,1);
    KH = kcenter(KH);
    KH = knorm(KH);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    epsionset = [0.1:0.1:0.9];
	for ie =1:length(epsionset)
		fprintf('%s: missing_ratio: %d ,iter: %d\n',dataName,epsionset(ie)*100,iter);         

		load([path,'*\generateAbsentMatrix\',dataName,'_missingRatio_',num2str(epsionset(ie)),'_missingIndex_iter_',num2str(iter),'.mat'],'S');
		qnorm = 2;
		%%%%--incomplete-SMKKM  TPAMI%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
		tic
		[H_normalized1,gamma1,obj1,~] = incompleteSimpleMKKM(KH,S,numclass,options);
		[res_mean(:,1),res_std(:,1)] = myNMIACCV2(H_normalized1,Y,numclass);
		timecost(1) = toc;
		clear H_normalized19
		tic;
	end
end
