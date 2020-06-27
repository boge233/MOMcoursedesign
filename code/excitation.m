function output = excitation(initphase,edge,trip)
%point是观察点
%输出为观察点处的电场
%输出矩量法右边的向量Vm
global pointtable;
global triangletable;
global w;
global k;
mag = 1;
j = complex(0,1);
ri = [0,0,1];
alpha = [0.3333 0.6 0.2 0.2];
beta = [0.3333 0.2 0.6 0.2];
gamma = [0.3332 0.2 0.2 0.6];
weight = [-0.5625 0.5208333 0.5208333 0.5208333];
%基函数积分点确定
a1 = triangletable(trip(1),1);
b1 = triangletable(trip(1),2);
c1 = triangletable(trip(1),3);%正三角形顶点索引
a2 = triangletable(trip(2),1);
b2 = triangletable(trip(2),2);
c2 = triangletable(trip(2),3);%负三角形顶点索引
point1 = zeros(4,3);
point2 = zeros(4,3);
point = pointtable;
for i= 1:4
    point1(i,:) =  gamma(i)*point(a1,:)+alpha(i)*point(b1,:)+beta(i)*point(c1,:);
    point2(i,:) =  gamma(i)*point(a2,:)+alpha(i)*point(b2,:)+beta(i)*point(c2,:);
end
I = 0;
%A1 = norm(cross((point(a1,:)-point(b1,:)),(point(c1,:)-point(b1,:))))/2;
%A2 = norm(cross((point(a2,:)-point(b2,:)),(point(c2,:)-point(b2,:))))/2;
for i = 1:4
    [rou1,s1] = rwg(point1(i,:),edge,triangletable(trip(1),:),triangletable(trip(2),:),pointtable);
    [rou2,s2] = rwg(point2(i,:),edge,triangletable(trip(1),:),triangletable(trip(2),:),pointtable);
    Lm = pointtable(edge(1),:) - pointtable(edge(2),:);
    Lm = norm(Lm);
    I = I + Lm/2*[1,0,0]*rou1'*weight(i)*exp(j*k*[0,0,1]*point1(i,:)'+initphase*j);
    I = I + Lm/2*[1,0,0]*rou2'*weight(i)*exp(j*k*[0,0,1]*point2(i,:)'+initphase*j);%极化为x方向
    if(I == 0)
    disp('**********');
    disp(rou1);
    disp(rou2);
    disp('**********');
    end
end
output = I;
end