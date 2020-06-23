function [output] = integral(edge1,edge2,trib,trit)
%本函数用来计算填充矩阵的元素,利用EFIE方程
%trib : 基函数三角形对索引
%trit : 检验函数三角形对索引
%triangletable  : 三角形表，由点的索引构成
%pointtable : 点的表
%edge1 : 基函数三角形公共边在pointtable中的索引，由两个点的索引构成
%edge2  :  势函数三角形公共边在pointtable中的索引，有两个点的索引构成
global pointtable;
global triangletable;
point = pointtable;
jc = complex(0,1);
global epsilon
global mu
global k
global w
%% 确定积分点 采用四点积分
alpha = [0.3333 0.6 0.2 0.2];
beta = [0.3333 0.2 0.6 0.2];
gamma = [0.3332 0.2 0.2 0.6];
weight = [-0.5625 0.5208333 0.5208333 0.5208333];
%基函数积分点确定
a1 = triangletable(trib(1),1);
b1 = triangletable(trib(1),2);
c1 = triangletable(trib(1),3);%正三角形顶点索引
a2 = triangletable(trib(2),1);
b2 = triangletable(trib(2),2);
c2 = triangletable(trib(2),3);%负三角形顶点索引
point1 = zeros(4,3);
point2 = zeros(4,3);
for i= 1:4
    point1(i,:) =  gamma(i)*point(a1,:)+alpha(i)*point(b1,:)+beta(i)*point(c1,:);
    point2(i,:) =  gamma(i)*point(a2,:)+alpha(i)*point(b2,:)+beta(i)*point(c2,:);
end
%势函数积分点确定
a1 = triangletable(trit(1),1);
b1 = triangletable(trit(1),2);
c1 = triangletable(trit(1),3);%正三角形顶点索引
a2 = triangletable(trit(2),1);
b2 = triangletable(trit(2),2);
c2 = triangletable(trit(2),3);%负三角形顶点索引
point3 = zeros(4,3);
point4 = zeros(4,3);
for i = 1:4
    point3(i,:) =  gamma(i)*point(a1,:)+alpha(i)*point(b1,:)+beta(i)*point(c1,:);
    point4(i,:) =  gamma(i)*point(a2,:)+alpha(i)*point(b2,:)+beta(i)*point(c2,:);
end
%% 利用RWG函数进行积分  
overlap = intersect(trib,trit);%重合的三角形
I = 0;
edgep1 = pointtable(edge1(1),:) - pointtable(edge1(2),:);
edgep2 = pointtable(edge2(1),:) - pointtable(edge2(2),:);
Ln = norm(edgep1);
Lm = norm(edgep2);
global test
if length(overlap) == 0
    Ippmm = 0 ;%两个正三角形之间的积分
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point1(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point3(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point1(i,:)-point3(z,:)))/norm(point1(i,:)-point3(z,:));
        end
    end
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point1(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point4(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point1(i,:)-point4(z,:)))/norm(point1(i,:)-point4(z,:));
        end
    end
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point2(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point3(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point2(i,:)-point3(z,:)))/norm(point2(i,:)-point3(z,:));
        end
    end
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point2(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point4(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point2(i,:)-point4(z,:)))/norm(point2(i,:)-point4(z,:));
        end
    end    
    I = Ippmm * Lm*Ln/(4*pi);
elseif length(overlap) == 1
    Ippmm = 0;
    v1 = pointtable(triangletable(overlap,1),:);
    v2 = pointtable(triangletable(overlap,2),:);
    v3 = pointtable(triangletable(overlap,3),:);
    index1 = find(trib==overlap);
    index2 = find(trit==overlap);
    sign = (index1-1.5)*(index2-1.5)*4;%当两个三角形不同种时为负
    vn = pointtable(setdiff(triangletable(overlap,:),edge1),:);
    vm = pointtable(setdiff(triangletable(overlap,:),edge2),:);
    %奇异性处理
    [I1,I2] = analyI12(v1,v2,v3,vn,vm,Lm,Ln,sign);
    pointb = zeros(4,3,2);
    pointt = zeros(4,3,2);
    pointb(:,:,1) = point1 ; pointb(:,:,2) = point2;
    pointt(:,:,1) = point3 ; pointt(:,:,2) = point4;
    for page1 = 1:2
        for  page2 = 1:2
            for i = 1:4
                for z = 1:4
                    [rou1,s1] = rwg(pointb(i,:,page1),edge1,triangletable(trib(page1),:),triangletable(trib(page2),:),pointtable);
                    [rou2,s2] = rwg(pointt(z,:,page2),edge2,triangletable(trit(page1),:),triangletable(trit(page2),:),pointtable);
                    edgep1 = pointtable(edge1(1),:) - pointtable(edge1(2),:);
                    edgep2 = pointtable(edge2(1),:) - pointtable(edge2(2),:);
                    if page1 == index1&&page2 == index2                        
                        Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*jc*k*(-1)*norm(edgep1)*norm(edgep2)/(4*pi);
                    else
                        Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*rou2'-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(pointb(i,:,page1)-pointt(z,:,page2)))/norm(pointb(i,:,page1)-pointt(z,:,page2))*norm(edgep1)*norm(edgep2)/(4*pi);
                    end
                end
            end
        end
    end
    I = Ippmm+I1+I2;
else
    Ippmm = 0;
    edgep1 = pointtable(edge1(1),:) - pointtable(edge1(2),:);
    edgep2 = pointtable(edge2(1),:) - pointtable(edge2(2),:);
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point1(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point3(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * rou1*(rou2') -s1*s2*jc/(w*epsilon))*(-jc*k)*norm(edgep1)*norm(edgep2)/(4*pi);
        end
    end
    vn = pointtable(setdiff(triangletable(trit(1),:),edge1),:);
    vm = pointtable(setdiff(triangletable(trib(1),:),edge2),:);
    v1 = pointtable(a1,:);v2 = pointtable(b1,:);v3 = pointtable(c1,:);
    sign = 1;
    [I1,I2] = analyI12(v1,v2,v3,vn,vm,Lm,Ln,sign);
    Ippmm = Ippmm + I1+I2;
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point1(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point4(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * (rou1*rou2')-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point1(i,:)-point4(z,:)))/norm(point1(i,:)-point4(z,:))*norm(edgep1)*norm(edgep2)/(4*pi);
        end
    end
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point2(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point3(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * (rou1*rou2')-s1*s2*jc/(w*epsilon))*exp(-jc*k*norm(point2(i,:)-point3(z,:)))/norm(point2(i,:)-point3(z,:))*norm(edgep1)*norm(edgep2)/(4*pi);
            
        end
    end
    for i = 1:4
        for z = 1:4
            [rou1,s1] = rwg(point2(i,:),edge1,triangletable(trib(1),:),triangletable(trib(2),:),pointtable);
            [rou2,s2] = rwg(point4(z,:),edge2,triangletable(trit(1),:),triangletable(trit(2),:),pointtable);
            Ippmm = Ippmm + weight(i)*weight(z)* (jc*w*mu/4 * (rou1*rou2')-s1*s2*jc/(w*epsilon))*(-jc*k)*norm(edgep1)*norm(edgep2)/(4*pi);
        end
    end
    vn = pointtable(setdiff(triangletable(trit(2),:),edge1),:);
    vm = pointtable(setdiff(triangletable(trib(2),:),edge2),:);
    v1 = pointtable(a2,:);v2 = pointtable(b2,:);v3 =pointtable(c2,:);
    sign = 1;
    [I1,I2] = analyI12(v1,v2,v3,vn,vm,Lm,Ln,sign);
    Ippmm = Ippmm + I1+I2;
    I = Ippmm ;
end
    output = I;
end
function [I1,I2] = analyI12(v1,v2,v3,vn,vm,Lm,Ln,sign)
%输入为点本身 sign是I1的符号
global w;
global mu;
global epsilon;
syms a b c d ;
ja = complex(0,1);
 %% 初始计算公式表达式
 l1_l1p = log((b + a^0.5*c^0.5)/(b-c+c^0.5*(a-2*b+c)^0.5))/(40*c^0.5) + log((-b + c+c^0.5*(a-2*b+c)^0.5)/(-b + a^0.5*c^0.5))/(40*c^0.5) +...
 (a^0.5*(a-2*b+c)^0.5 - c^0.5*(a-2*b+c)^0.5)/(120*(a-2*b+c)^1.5)+(2*a-5*b+3*c)*log(((a-b+a^0.5*(a-2*b+c)^0.5)*(c - b+c^0.5*(a-2*b+c)^0.5))/((b-a+(a^0.5*(a-2*b+c)^0.5))*(b-c+c^0.5*(a-2*b+c)^0.5)))/(120*(a-2*b+c)^1.5)+...
 +(2*a+b)*log(((b+a^0.5*c^0.5)*(a-b+a^0.5*d^0.5))/((-b+a^0.5*c^0.5)*(-a+b+a^0.5*d^0.5)))/(120*a^1.5)+(-a^0.5*c^0.5+a^0.5*d^0.5)/(120*a^1.5);
 l1_l2p = log((b+a^0.5*c^0.5)/(b-c+c^0.5*d^0.5))/(120*c^0.5)+ log((a-b+a^0.5*d^0.5)/(-b+a^0.5*c^0.5))/(120*a^0.5) + (-a^0.5*d^0.5+c^0.5*d^0.5)/(120*d^1.5)+...
     (2*a+2*b)*log((b+a^0.5*c^0.5)/(-a+b+a^0.5*d^0.5))/(120*a^1.5)+(a^0.5*d^0.5-c^0.5*d^0.5)/(120*d^1.5)+(a-3*b+2*c)*log((-b+c+c^0.5*d^0.5)/(-a+b+a^0.5*d^0.5))/(120*d^1.5)+...
     (-3*a^0.5*c^0.5+3*c^0.5*d^0.5)/(120*c^1.5)+(3*b+2*c)*log((-b+a+c^0.5*d^0.5)/(-b+a^0.5*c^0.5))/(120*c^1.5)+(-3*a^0.5*c^0.5+3*a^0.5*d^0.5)/(120*a^1.5)+...
     (2*a - 3*b +c)*log((a-b+a^0.5*d^0.5)/(b-c+c^0.5*d^0.5))/(120*d^1.5);
 l1_p = -log((-b+a^0.5*c^0.5)/(a-b+a^0.5*d^0.5))/(24*a^0.5)+log((b+a^0.5*c^0.5)/(b-c+c^0.5*d^0.5))/(24*c^0.5)+(-a^0.5*c^0.5+a^0.5*d^0.5)/(24*a^1.5)+...
     (a+b)*log((b+a^0.5*c^0.5)/(-a+b+a^0.5*d^0.5))/(24*a^1.5)+log((a-b+a^0.5*d^0.5)/(b-c+c^0.5*d^0.5))/(24*d^0.5)-log((-b + a^0.5*c^0.5)/(-b+c+c^0.5*d^0.5))/(12*c^0.5)+...
     (a^0.5*d^0.5-c^0.5*d^0.5)/(24*d^1.5)+(a-3*b+2*c)*log((-b+c+c^0.5*d^0.5)/(-a+b+a^0.5*d^0.5))/(24*d^1.5);
 I2_f = log((a-b+a^0.5*d^0.5)*(b+a^0.5*c^0.5)/((-b+a^0.5*c^0.5)*(-a+b+a^0.5*d^0.5)))/(6*a^0.5)+log((b+a^0.5*c^0.5)*(-b+c+c^0.5*d^0.5)/((b-c+c^0.5*d^0.5)*(-b+a^0.5*c^0.5)))/(6*c^0.5)+...
     log((a-b+a^0.5*d^0.5)*(-b+c+c^0.5*d^0.5)/((b-c+c^0.5*d^0.5)*(-a+b+a^0.5*d^0.5)))/(6*d^0.5);

 %% 具体计算
 a = (v1 - v3)*(v1 - v3)';
 b = (v1 - v3)*(v1 - v2)';
 c = (v1 - v2)*(v1 - v2)';
 d = a-2*b+c;
 l1_l1pn = eval(l1_l1p);
 l1_l2pn = eval(l1_l2p);
 l1_pn = eval(l1_p);
 I2 = eval(I2_f);
 %为了计算其他的项需要将v2 v3调换位置，表现在这里就是交换a，c的值，
 I2 = ja/(w*epsilon)*Lm*Ln/(4*pi)*I2*sign; 
 temp = a;
 a = c;
 c = temp;
 l2_l2pn = eval(l1_l1p);
 l1p_l2n = eval(l1_l2p);
 l2_pn = eval(l1_p);
 a11 = v1*v1' ;a12 = v1*v2';a13 = v1*v3' ;a22 = v2*v2';a23 = v2*v3';
 a33 = v3*v3' ;a1n = v1*vn';a1m = v1*vm' ;a2n = v2*vn';a2m = v2*vm';
 a3n = v3*vn' ;a3m = v3*vm';amn = vm*vn' ;
 I1 = l1_l1pn*(a11-2*a12+a22) + l1_l2pn*(a11 -  a13 -a12 +a23) +l2_l2pn*(a11 - 2*a13 +a33) + l1p_l2n*(a11-a12-a13+a23)+...
     l1_pn*(-a11+a1n+a12-a2n)+l2_pn*(-a11+a1n+a13-a3n)+l1_pn*(-a11+a1m+a12-a2m)+l2_pn*(-a11+a1m+a13-a3m)+a11-a1n-a1m+amn ;
 I1 = ja*w*mu*Lm*Ln/(4*pi)*I1;
 I1 = I1*sign;
end