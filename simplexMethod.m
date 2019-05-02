function [S,V,rA,rb,rXB] = simplexMethod(C,A,b)
%simplexMethod 单纯形法解线性规划,利用单纯形表，基本思路：
%   step1：找出一初始基及基可行解（一般都是找单位矩阵）
%   step2：判断该可行解是否为最优解，如果是，则停止运算；否就转step3
%   step3：换入换出转换基，转step2
%现在似乎只能解增添的松弛变量构成单位矩阵作为初始基（或者本身存在单位矩阵）的情况
%C--->目标系数向量,价值向量,按行向量输入
%A--->系数矩阵,按普通输入就行
%b--->资源向量，按列向量格式输入
%S--->返回的最优解，在没用写输入参数时会作为ans返回
%V--->返回的最优目标值
%rA--->返回运算结束时的A矩阵，供后续调用
%rb--->返回运算结束时的b向量，供后续调用
%rXB--->返回运算结束时的XB向量，供后续调用
%XB-->定位基变量
%theta--->换出准则向量
%Q--->为检验向量Q=CN-CB * N
lb = length(b); %b向量的长度
lc_list = 1:length(C);
flag = 0;
% 对于单位基矩阵不在矩阵的最后几列，可以人工交换列的次序
N = A;
CN = C;
XN = lc_list;
% 这里的操作以单位阵出现在系数矩阵右端为前提
[row_A, colomn_A] = size(A);
index_A = colomn_A - row_A;
XB = index_A+1:length(C);
XN(XB) = [];
B = A(:,XB);
N(:,XB) = [];
CB = C(XB);
CN(XB) = []; % C中剔除CB即为CN

while 1
    
    Z = CB * b;
    Q = CN - CB * N; % Q为检验向量
    index_Q_Positive = find(Q>0); % Q中可能存在多个正的元素
    
    for i = index_Q_Positive
        if max(N(:,i))<=0
            fprintf('此线性规划问题解无界！\n');
            S=[];V=[];
            flag = 1;
            break;
        end
    end
    if flag
        break; % 如果解无界的话，就跳出大循环，不用再求解了！
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%CTRL + I为批量缩进%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    X = zeros(1, length(C)); % 先把解赋成全0，然后把基向量对应的位置覆盖掉
    n = 1;
    for i = XB
        X(i) = b(n);
        n = n + 1;
    end
    
    if max(Q) < 0 %此时有唯一的最优解，其实最后显示的时候不用显示松弛变量的值，它对函数值也无影响
        S = X;
        V = Z;
        rA = A;
        rb = b;
        rXB = XB;
        disp('最优解为：'),disp(S);
        disp('最优目标值为：'),disp(V);
        break;
    end 
    
    if max(Q) == 0
        %为了后面调用过程中不显示中间过程，把输出注释掉了
        fprintf('此线性规划问题有无穷多最优解\n');
        S = X;
        V = Z;
        rA = A;
        rb = b;
        rXB = XB;
        disp('其中一个最优解为：'),disp(S);
        disp('最优目标值为：'),disp(V);
        break;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 确定换入向量（索引）
    % 要小心这里Q中可能存在多个最大值相同（还不是退化）
    % 当存在多个相同最大值时，就默认取第一个
    index_Q_Max = find(Q == max(Q));
    index_Q_Max = index_Q_Max(1);
    in_A = N(:,index_Q_Max); % in_A为换入向量
    % 确定换出向量（索引）
    theta = zeros(1,lb); %先初始化theta，防止下面有些未赋值的情况
    for i = 1:length(theta) %这里不用计算分母为0或为负数的情况
        if in_A(i) < 1e-8
            continue;
        end
        theta(i) = b(i) / in_A(i);
    end
    % 因为上面theta初始化为零向量，所以这里找的是最小的正theta
    index_theta_Min = find(theta == min(theta(theta~=0))); % 这里至少会有一项正的theta，因为没有的话已经被选中为上面的解无界的情况了
    pivot = N(index_theta_Min, index_Q_Max); % pivot即为系数矩阵中初等行变换的主元
    
    XB(index_theta_Min) = XN(index_Q_Max); % 对XB进行更新，并更新对应的系数
    CB(index_theta_Min) = CN(index_Q_Max);
    XN = lc_list;
    XN(XB) = []; % X中剔除XB即为XN
    % 下面其实也可以这样写,如果先求出XN的话
    %CN = C(XN);
    CN = C;
    CN(XB) = []; % C中剔除CB即为CN
    
    
    %%%%%%%%%接下来要做初等行变换了%%%%%%%%%
    A(index_theta_Min, :) = A(index_theta_Min, :) / pivot;
    b(index_theta_Min) = b(index_theta_Min) / pivot;
  
    for i = 1:row_A
        if i==index_theta_Min
            continue;
        end
        times = N(i,index_Q_Max);
        A(i,:) = A(i,:) - times * A(index_theta_Min, :);
        b(i) = b(i) - times * b(index_theta_Min);
    end
    
    B = A(:,XB);
    % 下面也可以这样写
    % N = A(:,XN);
    N = A;
    N(:,XB) = [];
end
end