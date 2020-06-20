function  init()
%UNTITLED 此处显示有关此函数的摘要
%   此处显示详细说明
global epsilon
global mu
global c
global k
global w
epsilon = 8.854187817*1e-12;
mu = 4*pi*1e-7;
c = 3*10*1e8;
f = 2e9;
w = f*2*pi;
lamda = c/f;
k = 2*pi/lamda;
end

