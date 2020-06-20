function [rwg1,gradrwg] = rwg(point,edge,trianglep,trianglem,pointtable)
%RWG ���������
%   �˴���ʾ��ϸ˵��
%����Ϊ�������κ͸������ε������㣬point�Ǽ�������ĵ㣬edge�ǹ����ߡ�
%���еĵ��point�ⶼ�ǵ㼯������
vertexp = setdiff(trianglep,edge);
vertexm = setdiff(trianglem,edge);%�����������εö��㣬setdiff ���󲹼�����
%% �жϵ��Ƿ������������� ͨ������ж�
vp = point - pointtable(vertexp,:);
ve1 = pointtable(edge(1),:) - pointtable(vertexp,:);
ve2 =  pointtable(edge(2),:) - pointtable(vertexp,:);
tag1 = cross(vp,ve2) * cross(ve1,ve2)';
e1v = -ve1;
e1p =  point - pointtable(edge(1),:);
e1e2 = pointtable(edge(2),:) - pointtable(edge(1),:);
tag2 = cross(e1p,e1e2)*cross(e1v,e1e2)';

%% �жϵ��Ƿ��ڸ���������
mvp = point - pointtable(vertexp,:);
mve1 = pointtable(edge(1),:) - pointtable(vertexm,:);
mve2 =  pointtable(edge(2),:) - pointtable(vertexm,:);
tag3 = cross(mvp,mve2) * cross(mve1,mve2)';
me1v = -mve1;
me1p =  point - pointtable(edge(1),:);
me1e2 = pointtable(edge(2),:) - pointtable(edge(1),:);
tag4 = cross(me1p,me1e2)*cross(me1v,me1e2)';
%% ������һ����ϵ�����ں����⣬�����rwg�������Ѿ�����ԭʼ������
if tag1>=0 && tag2>=0
    A = abs(cross(e1v,e1e2))/2;
    pv = -vp;
    current = 0.5 * pv ;
    gradrwg = -1;
elseif tag3>=0 && tag4>=0
    A = abs(cross(melv,me1e2))/2;
    current = 0.5 * mvp ;
    gradrwg = 1;
else
    current = [0,0,0];
    gradrwg = 0;
end
rwg1 = current;
end