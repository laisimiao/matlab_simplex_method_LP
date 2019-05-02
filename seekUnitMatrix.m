function [ indexUnitMatrix ] = seekUnitMatrix( inputMatrix )
%seekUnitMatrix �Ӿ��󣨿��Բ��Ƿ������ҵ����Թ��ɵ�λ�������
%parameter:inputMatrix---->�������
%return:indexUnitMatrix--->����˳�򹹳ɵ�λ������У��������ɵ�λ�����򷵻�0

indexUnitMatrix = [];
[row_inputMatrix,colomn_inputMatrix] = size(inputMatrix);
if row_inputMatrix <= colomn_inputMatrix
    for i = 1:row_inputMatrix
        E = zeros(row_inputMatrix,1);
        E(i) = 1;
        for j = 1:colomn_inputMatrix
            if all(inputMatrix(:,j)==E)
                indexUnitMatrix = [indexUnitMatrix,j];
                break; % ���ܴ��ڶ�������ͬ�ģ�ֻҪһ�е���������
            end
        end
    end
    
elseif row_inputMatrix > colomn_inputMatrix
    inputMatrix = inputMatrix';
    indexUnitMatrix = seekUnitMatrix( inputMatrix ); % �ݹ����

end
if isempty(indexUnitMatrix)
    indexUnitMatrix = 0;
end
end

