function [ optimal_solution,optimal_value ] = twoPhaseMethod( C,A,b,num_AV )
%twoPhaseMethod �������׶η�������֮ǰ�ĵ����η�����LP��������
%C,A,b���ǻ�Ϊ��׼��֮���ϵ��
%parameter��C--->Ŀ�����ϵ������,�����������룬�����˹�����ϵ��Ϊ-1
%parameter��A--->ϵ������,����ͨ�����������,��������֪A������������С������
%parameter��b--->��Դ����������������ʽ����
%parameter��num_AV--->���ɱ�׼�����˹������ĸ�����0<=num_AV<=row_A��
%return��optimal_solution--->���Ž�
%return��optimal_value--->����ֵ

%[row_A,colomn_A] = size(A);
length_C = length(C);
% flag = 1;
% �����˹�����֮������ĸ���
num_V = length_C - num_AV;
% if num_AV == 0
%     flag = 0;
% end

% ��һ�׶Σ�����֮ǰ��simplexMethod.m������ԭ�����Ƿ��п��н�
% �����ص�����ֵΪ0�������ڶ��׶Σ�����ԭ�����޿��н�
% ����ֻͨ������ɳڱ������ܹ��ɵ�λ����ģ�����Ĭ�Ͻ���ڶ��׶Σ��϶��н⣩
% ��Ϊ��ʵ��һ�׶�Ҳֻ�Ǹ��ڶ��׶���һ����ʼ�⣬����������оͲ���Ҫ��һ�׶���
C11 = zeros(1,num_V);
C12 = -1 * ones(1,num_AV);
C1 = [C11 C12];
% �����֮ǰ�ĵ����η��Ĵ��������ʽ��Ҫ��Ҳ�ǿ��Եģ����ǱȽϴ�ģ�
% ��������÷���֮ǰ����ϵ����������ǵ�λ����ĸ�ʽ
% ���¾��ǰ�ϵ������λ����ƣ����ƶ���Ӧ��CB��XB����������һ�¶�����ֵû��Ӱ�죬
% ֻ�����Ž�ĳЩ˳����һ�£���ʵ�������ƶ��˹�������
index_A_UnitMatrix = seekUnitMatrix(A);
move_index_A = index_A_UnitMatrix(1:length(index_A_UnitMatrix)-num_AV);
length_move_index_A = length(move_index_A);
% ����Ľ���������ְѻ��õ��ֻ���ȥ�����⣬����ͻ������
for i=1:length_move_index_A
    % ����A��C1��C
    A(:,[move_index_A(i),num_V-length_move_index_A+i]) = A(:,[num_V-length_move_index_A+i,move_index_A(i)]);
    C1([move_index_A(i),num_V-length_move_index_A+i]) = C1([num_V-length_move_index_A+i,move_index_A(i)]); %������ʵҪ�õ���C1
    C([move_index_A(i),num_V-length_move_index_A+i]) = C([num_V-length_move_index_A+i,move_index_A(i)]); % C������Ҳ������Ϊ�˱���һ��
end
[rS,rV,rA,rb,rXB] = simplexMethod(C1,A,b);
for i=1:length_move_index_A
    % ����rA
    rA(:,[move_index_A(i),num_V-length_move_index_A+i]) = rA(:,[num_V-length_move_index_A+i,move_index_A(i)]);
end
if rV ~= 0
    disp('ԭ�����޿��н⣡');
    optimal_solution = [];
    optimal_value = [];
    
end

if rV == 0
    disp('ԭ�����п��н⣡');

% �ڶ��׶Σ�ȥ���˹���������ԭĿ�꺯��ϵ�������õ�һ�׶εõ���ʼ�����α�
AA = rA(:,1:num_V);
[row_AA,colomn_AA] = size(AA);
C = C(1:num_V);
% ��ʱ��Ϊ���ʺ�simplexMethod.m�������ʽ���ֵý�����
index_AA_UnitMatrix = seekUnitMatrix(AA);
W = 1:num_V;
W(index_AA_UnitMatrix) = [];
not_index_AA_UnitMatrix = W;
length_move_index_AA = length(index_AA_UnitMatrix);
AAA=zeros(row_AA,colomn_AA);
CCC=zeros(1,length(C));
for i=1:length_move_index_AA
    % �ڴյ����η��������ʽ
    AAA(:,colomn_AA-length_move_index_AA+i) = AA(:,index_AA_UnitMatrix(i));
    CCC(colomn_AA-length_move_index_AA+i) = C(index_AA_UnitMatrix(i));
end
for j = 1:length(not_index_AA_UnitMatrix)
    AAA(:,j) = AA(:,not_index_AA_UnitMatrix(j));
    CCC(j) = C(not_index_AA_UnitMatrix(j));
end
[finalS,finalV,~,~,~] = simplexMethod(CCC,AAA,rb);
finalSS = zeros(1,length(finalS));
% ���ڵõ��Ľ⻹�ý���˳�򣬺�����֮ǰ�Ľ��˳���Ӧ
% ����1
for j = 1:length(not_index_AA_UnitMatrix)
    finalSS(j) = finalS(not_index_AA_UnitMatrix(j));
end
% ����2
for i=1:length_move_index_AA
    finalSS(index_AA_UnitMatrix(i)) = finalS(colomn_AA-length_move_index_AA+i);
end
% ����3
disp('ԭ�������Ž�Ϊ��'),disp(finalSS);
disp('ԭ��������ֵΪ��'),disp(finalV);
end

end

