function [rwg1,gradrwg] = rwg(point,edge,trianglep,trianglem,pointtable)
%RWG 计算基函数
%   此处显示详细说明
%输入为正三角形和负三角形的三个点，point是计算电流的点，edge是公共边。
%所有的点除point外都是点集的索引
vertexp = setdiff(trianglep,edge);
vertexm = setdiff(trianglem,edge);%求两个三角形得顶点，setdiff 是求补集函数
%% 判断点是否在正三角形内 通过叉乘判断
vp = point - pointtable(vertexp,:);
ve1 = pointtable(edge(1),:) - pointtable(vertexp,:);
ve2 =  pointtable(edge(2),:) - pointtable(vertexp,:);
tag1 = cross(vp,ve2) * cross(ve1,ve2)';
e1v = -ve1;
e1p =  point - pointtable(edge(1),:);
e1e2 = pointtable(edge(2),:) - pointtable(edge(1),:);
tag2 = cross(e1p,e1e2)*cross(e1v,e1e2)';

%% 判断点是否在负三角形内
mvp = point - pointtable(vertexp,:);
mve1 = pointtable(edge(1),:) - pointtable(vertexm,:);
mve2 =  pointtable(edge(2),:) - pointtable(vertexm,:);
tag3 = cross(mvp,mve2) * cross(mve1,mve2)';
me1v = -mve1;
me1p =  point - pointtable(edge(1),:);
me1e2 = pointtable(edge(2),:) - pointtable(edge(1),:);
tag4 = cross(me1p,me1e2)*cross(me1v,me1e2)';
%% 计算结果一部分系数放在函数外，这里的rwg基函数已经不是原始函数。
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