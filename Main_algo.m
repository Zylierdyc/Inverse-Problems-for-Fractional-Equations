
clear,clc

load('a1_true_6_12rd40_1.mat')
L = 1;
M = 80; % 空间网格数目
h = L/M; % 空间步长
x = h:h:L-h; % x方向空间节点
y = h:h:L-h; % y方向空间节点

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

labmda = (79^2)^(-4*(1+0)/(4*(1+0)+2))

fractional_order = 1.2;
a_1 = zeros(M-1,M-1);

% M1 = 80;
% h1 = L/M1;
% x1 = h1:h1:L-h1; 
% y1 = h1:h1:L-h1; 
% for i = 1:M1-1
%     for j = 1:M1-1
%         x_r = x(find(x>x1(i),1));
%         x_l = x(find(x>x1(i),1)-1);
%         y_u = y(find(y>y1(j),1));
%         y_l = y(find(y>y1(j),1)-1);
%         temp1 = (x1(i)-x_l)/h*u_end(find(x>x1(i),1),find(y>y1(j),1))...
%             +(x_r-x1(i))/h*u_end(find(x>x1(i),1)-1,find(y>y1(j),1));
%         temp2 = (x1(i)-x_l)/h*u_end(find(x>x1(i),1),find(y>y1(j),1)-1)...
%             +(x_r-x1(i))/h*u_end(find(x>x1(i),1)-1,find(y>y1(j),1)-1);        
%         a1_true_1(i,j) = (y_u-y1(j))/h*temp2+(y1(j)-y_l)/h*temp1;
%     end
% end

count = 0;
tic
beta = 10^3; 

while 1
    u_end = forwards(fractional_order,a_1);
    temp_1 = reshape(u_end - a1_true_6_12rd40_1,(M-1)*(M-1),1);
    temp_2 = 1/((M-1)*(M-1))*inv(B)*temp_1;
    u_end_u_delta = reshape(temp_2,M-1,M-1);
    w_ini = backwards(fractional_order,u_end_u_delta);

    while 1
        a_1_update = a_1 - beta*(labmda*a_1-w_ini); 
        u_end_test = forwards(fractional_order,a_1_update);
        temp_1_1 = reshape(a_1_update,(M-1)*(M-1),1);
        temp_2_2 = reshape(a_1,(M-1)*(M-1),1);
        
        if (norm(u_end_test-a1_true_6_12rd40_1,'fro')/(M-1))^2 + labmda*(temp_1_1'*B*temp_1_1) < (norm(u_end-a1_true_6_12rd40_1,'fro')/(M-1))^2 + labmda*(temp_2_2'*B*temp_2_2)
            break
        else
            beta = beta/2;
        end
    end
    
    a_1_update_lie = reshape(a_1_update,(M-1)*(M-1),1);
    labmda_update = ( (79^2)^(-0.5)*norm(u_end_test-a1_true_6_12rd40_1,'fro')/(M-1)/sqrt(a_1_update_lie'*B*a_1_update_lie) )^( 8*(1+0)/(4*(1+0)+2) );
    
    if abs(labmda-labmda_update) < 10^(-6)       
        break
    else
        labmda = labmda_update
        a_1 = a_1_update;
    end
end


toc
count

figure
a_1_full = zeros(M+1,M+1);
a_1_full(2:M,2:M) = a_1;
surf([0;x';1],[0,y,1],a_1_full)
% legend('数值解')
shading interp

% figure
% surf(x',y,u_end)
% shading interp
% 
% figure
% surf(x1',y1,a1_true_1)
% shading interp


% figure
% b_1 = (exp(4*x).*sin(pi*x))'*sin(2*pi*y); % 精确初值
% surf(x',y,b_1)
% shading interp
% 
% norm(a_1-b_1)/norm(b_1)

% c = zeros(M-1,M-1);
% for i = 1:M-1
%     for j = 1:M-1
%         if x(i) >= 0.25 && x(i) <= 0.75 && y(j) >= 0.25 && y(j) <= 0.75
%             c(i,j) = 1;
%         end
%     end
% end
% 
% norm(a_1-c)/norm(c)
% 
% figure
% surf(x',y,c)
% shading interp
% b_2 = 8*(x.^1.01.*(x-1))'*sin(pi*y); % 初始猜测

b_2 = 7*(x.^0.75.*(x-1))'*sin(2*pi*y); % 初始猜测

temp_3 = reshape(a_1-b_2,(M-1)*(M-1),1);
temp_4 = reshape(b_2,(M-1)*(M-1),1);
sqrt(temp_3'*inv(A)*temp_3)/sqrt(temp_4'*inv(A)*temp_4)

figure
b_2_full = zeros(M+1,M+1);
b_2_full(2:M,2:M) = b_2;
surf([0;x';1],[0,y,1],b_2_full)
% legend('精确解')
shading interp


figure
surf(x',y,a_1-b_2)
legend('误差')
shading interp

% b_3 = (exp(4*x).*sin(pi*x))'*(1-abs(2*y-1)); % 初始猜测
% norm(a_1-b_3)/norm(b_3)
% figure
% surf(x',y,b_3)
% shading interp

% figure
% surf(x',y,c)
% shading flat
% 
% figure
% surf(x',y,u_end)
% shading interp
% % 
% figure
% surf(x1',y1,a1_true_1)
% shading interp
% % 
% figure
% surf(x1',y1,a1_true_2_2_5)
% shading interp

