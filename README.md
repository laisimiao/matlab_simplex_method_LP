# matlab_simplex_method_LP
% author：赖思邈
% date：2019/5/1
说明：
基于第一版本的改进版
优点：
	1、能解基本的LP问题
	2、输入更加自由，只要按标准化的C,A,b输入即可（C中松弛/剩余/人工变量的系数为0）
	3、加入了两阶段法，对解的情况进行全覆盖
	
缺点：
	1、输入仍然受限，多了个人工变量个数参数num，需自己输入；在没有人工变量的时候也需要加入人工变量
	为了输入格式的要求
	2、对于退化问题，运用Bland法则还没加进去，后一版本考虑加进；
	3、对于最后工作空间显示结果的说明：
		1、最上面的“此线性规划问题有无穷多最优解、其中一个最优解为：、最优目标值为：”不予理会
		但中间若还有这句话，表明该LP问题真的有无穷多解
		2、中间若出现“最优解为：、最优目标值为：”字眼，则为有唯一解
		3、最后若出现“此线性规划问题解无界！”字眼，并有少量报错信息，则表示该LP问题解无界