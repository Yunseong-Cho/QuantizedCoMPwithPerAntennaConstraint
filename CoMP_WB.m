function [P_Tc, SINR, count_Tc, lambda, D] = CoMP_WB( SP, G, n, method, lambda, noPA)


% First, computation of allocated power
% Second, computation of SINR

% - Input
%     SP: System Parameters
%     H: Channel
%     n: index of target SINR array
%     method: specifies how to calculate the powers ('joint','percell','det')
%     res: Resolution ('infinite', 'low')
%
% - Output
%     P_Tc: Average power
%     SINR: SINR of each user in each cell
%     count_Tc: Average number of iterations to converge

Nc = SP.Nc;
Nu = SP.Nu;
b = SP.b;
TcMax = SP.TcMax ;


initlambda = lambda; %power initialization

switch method
    case 'joint'
        if SP.perantenna == 1
            [lambda, count, SINR, D] = Algo_joint_WB_PA(SP, G, n, initlambda, noPA);
        elseif SP.perantenna == 0
            [lambda, count, SINR] = Algo_joint_WB(SP, G, n, initlambda);
        end
            
    case 'percell'
        [lambda, count, SINR] = Algo_percell_wb(SP, G, n, initlambda);            
end
  

SINR = abs(SINR);
P_Tc = sum(lambda(:));
% P_Tc = sum(lambda(1,1,1))*2*2*64;


count_Tc = 1;%nanmean(count);

end


