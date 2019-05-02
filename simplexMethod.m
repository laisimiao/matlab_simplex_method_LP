function [S,V,rA,rb,rXB] = simplexMethod(C,A,b)
%simplexMethod �����η������Թ滮,���õ����α�����˼·��
%   step1���ҳ�һ��ʼ���������н⣨һ�㶼���ҵ�λ����
%   step2���жϸÿ��н��Ƿ�Ϊ���Ž⣬����ǣ���ֹͣ���㣻���תstep3
%   step3�����뻻��ת������תstep2
%�����ƺ�ֻ�ܽ�������ɳڱ������ɵ�λ������Ϊ��ʼ�������߱�����ڵ�λ���󣩵����
%C--->Ŀ��ϵ������,��ֵ����,������������
%A--->ϵ������,����ͨ�������
%b--->��Դ����������������ʽ����
%S--->���ص����Ž⣬��û��д�������ʱ����Ϊans����
%V--->���ص�����Ŀ��ֵ
%rA--->�����������ʱ��A���󣬹���������
%rb--->�����������ʱ��b����������������
%rXB--->�����������ʱ��XB����������������
%XB-->��λ������
%theta--->����׼������
%Q--->Ϊ��������Q=CN-CB * N
lb = length(b); %b�����ĳ���
lc_list = 1:length(C);
flag = 0;
% ���ڵ�λ�������ھ��������У������˹������еĴ���
N = A;
CN = C;
XN = lc_list;
% ����Ĳ����Ե�λ�������ϵ�������Ҷ�Ϊǰ��
[row_A, colomn_A] = size(A);
index_A = colomn_A - row_A;
XB = index_A+1:length(C);
XN(XB) = [];
B = A(:,XB);
N(:,XB) = [];
CB = C(XB);
CN(XB) = []; % C���޳�CB��ΪCN

while 1
    
    Z = CB * b;
    Q = CN - CB * N; % QΪ��������
    index_Q_Positive = find(Q>0); % Q�п��ܴ��ڶ������Ԫ��
    
    for i = index_Q_Positive
        if max(N(:,i))<=0
            fprintf('�����Թ滮������޽磡\n');
            S=[];V=[];
            flag = 1;
            break;
        end
    end
    if flag
        break; % ������޽�Ļ�����������ѭ��������������ˣ�
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%CTRL + IΪ��������%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    X = zeros(1, length(C)); % �Ȱѽ⸳��ȫ0��Ȼ��ѻ�������Ӧ��λ�ø��ǵ�
    n = 1;
    for i = XB
        X(i) = b(n);
        n = n + 1;
    end
    
    if max(Q) < 0 %��ʱ��Ψһ�����Ž⣬��ʵ�����ʾ��ʱ������ʾ�ɳڱ�����ֵ�����Ժ���ֵҲ��Ӱ��
        S = X;
        V = Z;
        rA = A;
        rb = b;
        rXB = XB;
        disp('���Ž�Ϊ��'),disp(S);
        disp('����Ŀ��ֵΪ��'),disp(V);
        break;
    end 
    
    if max(Q) == 0
        %Ϊ�˺�����ù����в���ʾ�м���̣������ע�͵���
        fprintf('�����Թ滮��������������Ž�\n');
        S = X;
        V = Z;
        rA = A;
        rb = b;
        rXB = XB;
        disp('����һ�����Ž�Ϊ��'),disp(S);
        disp('����Ŀ��ֵΪ��'),disp(V);
        break;
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % ȷ������������������
    % ҪС������Q�п��ܴ��ڶ�����ֵ��ͬ���������˻���
    % �����ڶ����ͬ���ֵʱ����Ĭ��ȡ��һ��
    index_Q_Max = find(Q == max(Q));
    index_Q_Max = index_Q_Max(1);
    in_A = N(:,index_Q_Max); % in_AΪ��������
    % ȷ������������������
    theta = zeros(1,lb); %�ȳ�ʼ��theta����ֹ������Щδ��ֵ�����
    for i = 1:length(theta) %���ﲻ�ü����ĸΪ0��Ϊ���������
        if in_A(i) < 1e-8
            continue;
        end
        theta(i) = b(i) / in_A(i);
    end
    % ��Ϊ����theta��ʼ��Ϊ�����������������ҵ�����С����theta
    index_theta_Min = find(theta == min(theta(theta~=0))); % �������ٻ���һ������theta����Ϊû�еĻ��Ѿ���ѡ��Ϊ����Ľ��޽�������
    pivot = N(index_theta_Min, index_Q_Max); % pivot��Ϊϵ�������г����б任����Ԫ
    
    XB(index_theta_Min) = XN(index_Q_Max); % ��XB���и��£������¶�Ӧ��ϵ��
    CB(index_theta_Min) = CN(index_Q_Max);
    XN = lc_list;
    XN(XB) = []; % X���޳�XB��ΪXN
    % ������ʵҲ��������д,��������XN�Ļ�
    %CN = C(XN);
    CN = C;
    CN(XB) = []; % C���޳�CB��ΪCN
    
    
    %%%%%%%%%������Ҫ�������б任��%%%%%%%%%
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
    % ����Ҳ��������д
    % N = A(:,XN);
    N = A;
    N(:,XB) = [];
end
end