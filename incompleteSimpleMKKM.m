function [Hstar,Sigma,obj,KH] = incompleteSimpleMKKM(KH,S,numclass,option)

numker = size(KH,3);
num = size(KH,1);

if ~isfield(option,'goldensearch_deltmax')
    option.goldensearch_deltmax=5e-2;
end
if ~isfield(option,'goldensearchmax')
    optiongoldensearchmax=1e-8;
end
if ~isfield(option,'firstbasevariable')
    option.firstbasevariable='first';
end
%%-------------------------
%% Initializing Missing Elememts of KHs
%%-------------------------
KH = initializeKH(KH,S);
flag =1;
iter = 1;
while flag
    %% MKKM with imputed kernel matrix
    [Hstar,Sigma,obj1] = simpleMKKM(KH,numclass,option);
    obj(iter) = obj1(end);
    %% update kernel matrix with Hstar and Sigma
    for p = 1 : numker
        mis_set = S{p}.indx;
        obs_set = setdiff(1:num, mis_set);
        KH(:,:,p) = absentKernelImputationV3(Hstar,KH(obs_set,obs_set,p),mis_set);
    end
    if iter >=2 &&  abs(obj(iter)- obj(iter-1))/obj(iter)<1e-3
        flag =0;
    else
        iter = iter + 1;
    end
end