clear all;
%��ȡstl ��ά�ļ����ڲ������ʺõ�����
%ʹ��matlab�Դ����ʷ�Ҳ���ԣ�����Ҫ��̡�
init();
model = stlread('dazuoye.stl');
[m,n] = size(model.ConnectivityList);
edgetable = [];
%% ��ȡ�ߣ�һ������������ȷ��
for i=1:m
    triangleedge = combntns(model.ConnectivityList(i,:),2);
    edgetable = [edgetable;triangleedge];
end
edgetable = (sort(edgetable'))';
edgetable = unique(edgetable,'row');
%% �ɹ�����ȷ�������ζ�
trp = [];
[m,n] = size(edgetable);
for i =1:m
    [ri1,ci] = find(model.ConnectivityList==edgetable(i,1));
    [ri2,ci] = find(model.ConnectivityList==edgetable(i,2));
    ci = intersect(ri1,ri2);
    trp = [trp;ci'];%�����ζ�����,����n�У����зֱ�Ϊ���������������ӱ��е������������빫����˳�ζ�Ӧ��
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