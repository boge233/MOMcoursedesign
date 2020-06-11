clear all;
%��ȡstl ��ά�ļ����ڲ������ʺõ�����
%ʹ��matlab�Դ����ʷ�Ҳ���ԣ�����Ҫ��̡�
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
%% �ɱ�ȷ�������ζ�
trp = []
[m,n] = size(edgetable);
for i =1:m
    [ri1,ci] = find(model.ConnectivityList==edgetable(i,1));
    [ri2,ci] = find(model.ConnectivityList==edgetable(i,2));
    ci = intersect(ri1,ri2);
    trp = [trp;ci'];%�����ζ�����
end