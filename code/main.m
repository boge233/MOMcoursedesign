clear all;
%读取stl 三维文件，内部已有剖好的网格，
%使用matlab自带的剖分也可以，但需要编程。
init();
model = stlread('dazuoye.stl');
[m,n] = size(model.ConnectivityList);
edgetable = [];
%% 提取边，一条边由两个点确定
for i=1:m
    triangleedge = combntns(model.ConnectivityList(i,:),2);
    edgetable = [edgetable;triangleedge];
end
edgetable = (sort(edgetable'))';
edgetable = unique(edgetable,'row');
%% 由公共边确定三角形对
trp = [];
[m,n] = size(edgetable);
for i =1:m
    [ri1,ci] = find(model.ConnectivityList==edgetable(i,1));
    [ri2,ci] = find(model.ConnectivityList==edgetable(i,2));
    ci = intersect(ri1,ri2);
    trp = [trp;ci'];%三角形对数组,两行n列，两行分别为正负三角形在连接表中的索引，并且与公共边顺次对应。
end
global triangletable
global pointtable
triangletable = model.ConnectivityList;
pointtable = model.Points;
global test;
test = 0;
Z = zeros(m,m);
V = zeros(m,1);
for i = 1:m
    V(i) = excitation(0,edgetable(i,:),trp(i,:));
    
end
for i = 1:m
    for j = 1:m
        Z(i,j) = integral(edgetable(i,:),edgetable(j,:),trp(i,:),trp(j,:));
    end
end
K = Z\V;