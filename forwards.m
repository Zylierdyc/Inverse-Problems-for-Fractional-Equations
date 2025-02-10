%% 线性分数阶波方程2D

function [u_end] = forwards(fractional_order,a_1)
%% 网格定义
beta = fractional_order; %导数阶数
alpha = beta-1; %降阶
sigma = 1-alpha/2; %偏移点
L = pi; %空间区间长度

M = 80; %空间区间数目
h = L/M; %空间步长
x = h:h:L-h; %x方向空间节点
y = h:h:L-h; %y方向空间节点

e = ones(M-1,1); %向量
C = spdiags([-1/3*e 8/3*e -1/3*e],[-1 0 1],M-1,M-1);
D = spdiags([-1/3*e -1/3*e -1/3*e],[-1 0 1],M-1,M-1);
I = eye(M-1);
ee = ones(M-1,1);
S = spdiags([ee ee],[-1 1],M-1,M-1);
A = kron(I,C)+kron(S,D); %刚度大矩阵

E = spdiags(h^2*[1/9*e 4/9*e 1/9*e],[-1 0 1],M-1,M-1); 
F = spdiags(h^2*[1/36*e 1/9*e 1/36*e],[-1 0 1],M-1,M-1);
B = kron(I,E)+kron(S,F); %质量大矩阵

z = [1000]; %时间区间数目
for i = 1:length(z)
    N = z(i);
    T = 1;
    tau = T/N; %时间步长
    t = 0:tau:T; %时间节点
       
    %源项 
    tt = (1-sigma)*t(1:end-1) + sigma*t(2:end);
    f = zeros(M-1,M-1,N);    
                   
    %数值解存储空间
    u_numerical = zeros(M-1,M-1,N+1);
    u_numerical(:,:,1) = zeros(M-1,M-1);
    v_numerical = zeros(M-1,M-1,N+1);
    v_numerical(:,:,1) = a_1;

    %第一层的数值解
       
    g = coe_a(tau,alpha,sigma,1);
    
    left = 2/tau*g(1)*B + sigma*A; %B质量矩阵 A刚度矩阵
       
    right1 = reshape(f(:,:,1),(M-1)*(M-1),1); %源项
    
    right2 = 2/tau*g(1)*u_numerical(:,:,1) + 2*g(1)*v_numerical(:,:,1); %u v历史项
    right2 = reshape(right2,(M-1)*(M-1),1);
              
    right3 = reshape(u_numerical(:,:,1),(M-1)*(M-1),1);
    right3 = -(1-sigma)*A*right3; %空间前项
       
    u_numerical_temp = left\( B*( right1 + right2 ) + right3 );    
    
    u_numerical(:,:,2) = reshape(u_numerical_temp,(M-1),(M-1));
    
    v_numerical(:,:,2) = 2/tau*(u_numerical(:,:,2) - u_numerical(:,:,1)) - v_numerical(:,:,1); %u v差分替换
       
    for j = 3:N+1
        
       g = coe_a(tau,alpha,sigma,j-1);
       
       left = (2*sigma+1)/2/tau/sigma*g(j-1)*B + sigma*A; %B质量矩阵 A刚度矩阵
       
       right1 = reshape(f(:,:,j-1),(M-1)*(M-1),1); %源项
              
       right2 = zeros(M-1,M-1);
       for jj = 2:j-1
           right2 =  right2 + v_numerical(:,:,jj)*(g(jj)-g(jj-1));%积分里求和的项
       end
       right2 = reshape(right2,(M-1)*(M-1),1);
       
       right3 = v_numerical(:,:,1)*g(1) + 2/tau*g(j-1)*u_numerical(:,:,j-1) ...
           - (2*sigma-1)/2/tau/sigma*g(j-1)*u_numerical(:,:,j-2) + (1-sigma)/sigma*g(j-1)*v_numerical(:,:,j-1);
       right3 = reshape(right3,(M-1)*(M-1),1); %降阶变量替换剩余项
       
       right4 = reshape(u_numerical(:,:,j-1),(M-1)*(M-1),1);
       right4 = -(1-sigma)*A*right4; %空间前项
                                   
       u_numerical_temp = left\( B*( right1 + right2 + right3 ) + right4 );
             
       u_numerical(:,:,j) = reshape(u_numerical_temp,M-1,M-1);
       
       v_numerical(:,:,j) = 1/2/tau/sigma*((2*sigma+1)*u_numerical(:,:,j) - 4*sigma*u_numerical(:,:,j-1) ...
           + (2*sigma-1)*u_numerical(:,:,j-2)) - (1-sigma)/sigma*v_numerical(:,:,j-1); %u v差分替换
    end   
           
end

u_end = u_numerical(:,:,end);

end

