function [ optimal_solution,optimal_value ] = twoPhaseMethod( C,A,b,num_AV )
%twoPhaseMethod 运用两阶段法，调用之前的单纯形法进行LP问题的求解
%C,A,b都是化为标准型之后的系数
%parameter：C--->目标变量系数向量,按行向量输入，其中人工变量系数为-1
%parameter：A--->系数矩阵,按普通逐行输入就行,经分析可知A几乎都是行数小于列数
%parameter：b--->资源向量，按列向量格式输入
%parameter：num_AV--->化成标准型中人工变量的个数（0<=num_AV<=row_A）
%return：optimal_solution--->最优解
%return：optimal_value--->最优值

%[row_A,colomn_A] = size(A);
length_C = length(C);
% flag = 1;
% 计算人工变量之外变量的个数
num_V = length_C - num_AV;
% if num_AV == 0
%     flag = 0;
% end

% 第一阶段：调用之前的simplexMethod.m来检验原问题是否有可行解
% 若返回的最优值为0，则进入第二阶段；否则，原问题无可行解
% 若是只通过添加松弛变量就能构成单位矩阵的，那样默认进入第二阶段（肯定有解）
% 因为其实第一阶段也只是给第二阶段找一个初始解，如果本来就有就不需要第一阶段了
C11 = zeros(1,num_V);
C12 = -1 * ones(1,num_AV);
C1 = [C11 C12];
% 因不想改之前的单纯形法的代码输入格式（要改也是可以的，就是比较大改）
% 所以这里得符合之前输入系数矩阵最后是单位矩阵的格式
% 以下就是把系数矩阵单位阵后移，并移动相应的CB和XB，这样交换一下对最优值没用影响，
% 只是最优解某些顺序换了一下（其实并不会移动人工变量）
index_A_UnitMatrix = seekUnitMatrix(A);
move_index_A = index_A_UnitMatrix(1:length(index_A_UnitMatrix)-num_AV);
length_move_index_A = length(move_index_A);
% 这里的交换不会出现把换好的又换出去的问题，下面就会出现了
for i=1:length_move_index_A
    % 交换A，C1，C
    A(:,[move_index_A(i),num_V-length_move_index_A+i]) = A(:,[num_V-length_move_index_A+i,move_index_A(i)]);
    C1([move_index_A(i),num_V-length_move_index_A+i]) = C1([num_V-length_move_index_A+i,move_index_A(i)]); %下面其实要用的是C1
    C([move_index_A(i),num_V-length_move_index_A+i]) = C([num_V-length_move_index_A+i,move_index_A(i)]); % C在这里也交换了为了保持一致
end
[rS,rV,rA,rb,rXB] = simplexMethod(C1,A,b);
for i=1:length_move_index_A
    % 交换rA
    rA(:,[move_index_A(i),num_V-length_move_index_A+i]) = rA(:,[num_V-length_move_index_A+i,move_index_A(i)]);
end
if rV ~= 0
    disp('原问题无可行解！');
    optimal_solution = [];
    optimal_value = [];
    
end

if rV == 0
    disp('原问题有可行解！');

% 第二阶段，去掉人工变量，还原目标函数系数，利用第一阶段得到初始单纯形表
AA = rA(:,1:num_V);
[row_AA,colomn_AA] = size(AA);
C = C(1:num_V);
% 此时，为了适合simplexMethod.m的输入格式，又得交换列
index_AA_UnitMatrix = seekUnitMatrix(AA);
W = 1:num_V;
W(index_AA_UnitMatrix) = [];
not_index_AA_UnitMatrix = W;
length_move_index_AA = length(index_AA_UnitMatrix);
AAA=zeros(row_AA,colomn_AA);
CCC=zeros(1,length(C));
for i=1:length_move_index_AA
    % 在凑单纯形法的输入格式
    AAA(:,colomn_AA-length_move_index_AA+i) = AA(:,index_AA_UnitMatrix(i));
    CCC(colomn_AA-length_move_index_AA+i) = C(index_AA_UnitMatrix(i));
end
for j = 1:length(not_index_AA_UnitMatrix)
    AAA(:,j) = AA(:,not_index_AA_UnitMatrix(j));
    CCC(j) = C(not_index_AA_UnitMatrix(j));
end
[finalS,finalV,~,~,~] = simplexMethod(CCC,AAA,rb);
finalSS = zeros(1,length(finalS));
% 现在得到的解还得交换顺序，好与最之前的解的顺序对应
% 交换1
for j = 1:length(not_index_AA_UnitMatrix)
    finalSS(j) = finalS(not_index_AA_UnitMatrix(j));
end
% 交换2
for i=1:length_move_index_AA
    finalSS(index_AA_UnitMatrix(i)) = finalS(colomn_AA-length_move_index_AA+i);
end
% 交换3
disp('原问题最优解为：'),disp(finalSS);
disp('原问题最优值为：'),disp(finalV);
end

end

