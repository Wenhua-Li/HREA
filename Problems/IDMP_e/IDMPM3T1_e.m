function varargout = IDMPM3T1_e(Operation,Global,input)
% <problem> <IDMP_e>
% Imbalanced Distance Minimization Problems
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright 2018-2019 Yiping Liu
% This is the code of benchmarks proposed in "Yiping Liu, Hisao Ishibuchi, 
% Gary G. Yen, Yusuke Nojima and Naoki Masuyama, Handling Imbalance Between 
% Convergence and Diversity in the Decision Space in Evolutionary Multi-
% Modal Multi-Objective Optimization, IEEE Transactions on Evolutionary 
% Computation, 2019, Early Access, DOI: 10.1109/TEVC.2019.2938557".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    persistent Points;
    switch Operation
        case 'init'
            Global.M          = 3;
            Global.M          = 3;
            Global.D          = 3;
            Global.D          = 3;
            Global.lower      = [-1,-1,-1];
            Global.upper      = [1,1,1];
            Global.operator   = @EAreal;          
            pgon = nsidedpoly(Global.M);
            psize = [.1,.1,.1,.1];
            center = [-.50,-.50;.50,-.50;.50,.50;-.50,.50]; 
            Points = NaN(Global.M,2,4);
            for i=1:4
                Points(:,:,i) = pgon.Vertices.*psize(i)+center(i,:);
            end
            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1) + repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            N = size(PopDec,1);
            PopObj = NaN(N,Global.M);
            for i=1:Global.M
                temp = pdist2(PopDec(:,1:2),reshape(Points(i,:,:),[2,4])');              
                temp(:,1) = temp(:,1)+1.*(abs(PopDec(:,3)+.6));
                temp(:,2) = temp(:,2)+2.*(abs(PopDec(:,3)+.2));
                temp(:,3) = temp(:,3)+1.*(abs(PopDec(:,3)-.2));
                temp(:,4) = temp(:,4)+2.*(abs(PopDec(:,3)-.6));
                PopObj(:,i) = min(temp,[],2);
            end
            
            index = PopDec(:,2)>0;
            PopObj(index,:) = PopObj(index,:) + 0.03;
            
            PopCon = [];            
            varargout = {input,PopObj,PopCon};
        case 'PF'
            [X,Y]     = ndgrid(linspace(0,1,ceil(sqrt(input))));
            ND        = inpolygon(X(:),Y(:),Points(:,1,3),Points(:,2,3));
            PopObj    = pdist2([X(ND),Y(ND)],Points(:,:,3));
            varargout = {PopObj};
    end
end