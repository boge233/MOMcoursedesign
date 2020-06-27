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
parfor i =1:m
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
    V(i) = excitation(pi/2,edgetable(i,:),trp(i,:));    
end
save('trit.mat','triangletable');
save('pointt.mat','pointtable');
parfor i = 1:m
    disp(i);
    for j = 1:m
        Z(i,j) = integral(edgetable(j,:),edgetable(i,:),trp(j,:),trp(i,:));        
    end
end
K = Z\V;
[l,n] = size(triangletable);
for h = 1:4
    for i = 1:l
        alpha = [0.3333 0.6 0.2 0.2];
        beta = [0.3333 0.2 0.6 0.2];
        gamma = [0.3332 0.2 0.2 0.6];
        tri = triangletable(i,:);%第i个三角形点索引
        obp = (alpha(h)*pointtable(tri(1),:)+beta(h)*pointtable(tri(2),:)+gamma(h)*pointtable(tri(3),:));%观察点坐标
        [ki,ci] = find(trp == i);%得到包含第i个三角形的三角形对
        v1 = pointtable(tri(2),:)-pointtable(tri(1),:);
        v2 = pointtable(tri(3),:)-pointtable(tri(1),:);
        A = norm(cross(v1,v2))/2;
        Jn = [0,0,0];
        Vm = 0;
        for ei = 1:3
            Lm = norm(pointtable(edgetable(ki(ei),1),:)-pointtable(edgetable(ki(ei),2),:));
            trv = trp(ki(ei),:);
            [J,cg] = rwg(obp,edgetable(ki(ei),:),triangletable(trv(1),:),triangletable(trv(2),:),pointtable);
            Jn = Jn + J*Lm/(A*2)*K(ki(ei));
        end
        Jnc = norm(Jn);
        scatter3(obp(1),obp(2),obp(3),10,abs(Vm),'filled');
        hold on;        
    end
end



