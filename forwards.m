%% ���Է����ײ�����2D

function [u_end] = forwards(fractional_order,a_1)
%% ������
beta = fractional_order; %��������
alpha = beta-1; %����
sigma = 1-alpha/2; %ƫ�Ƶ�
L = pi; %�ռ����䳤��

M = 80; %�ռ�������Ŀ
h = L/M; %�ռ䲽��
x = h:h:L-h; %x����ռ�ڵ�
y = h:h:L-h; %y����ռ�ڵ�

e = ones(M-1,1); %����
C = spdiags([-1/3*e 8/3*e -1/3*e],[-1 0 1],M-1,M-1);
D = spdiags([-1/3*e -1/3*e -1/3*e],[-1 0 1],M-1,M-1);
I = eye(M-1);
ee = ones(M-1,1);
S = spdiags([ee ee],[-1 1],M-1,M-1);
A = kron(I,C)+kron(S,D); %�նȴ����

E = spdiags(h^2*[1/9*e 4/9*e 1/9*e],[-1 0 1],M-1,M-1); 
F = spdiags(h^2*[1/36*e 1/9*e 1/36*e],[-1 0 1],M-1,M-1);
B = kron(I,E)+kron(S,F); %���������

z = [1000]; %ʱ��������Ŀ
for i = 1:length(z)
    N = z(i);
    T = 1;
    tau = T/N; %ʱ�䲽��
    t = 0:tau:T; %ʱ��ڵ�
       
    %Դ�� 
    tt = (1-sigma)*t(1:end-1) + sigma*t(2:end);
    f = zeros(M-1,M-1,N);    
                   
    %��ֵ��洢�ռ�
    u_numerical = zeros(M-1,M-1,N+1);
    u_numerical(:,:,1) = zeros(M-1,M-1);
    v_numerical = zeros(M-1,M-1,N+1);
    v_numerical(:,:,1) = a_1;

    %��һ�����ֵ��
       
    g = coe_a(tau,alpha,sigma,1);
    
    left = 2/tau*g(1)*B + sigma*A; %B�������� A�նȾ���
       
    right1 = reshape(f(:,:,1),(M-1)*(M-1),1); %Դ��
    
    right2 = 2/tau*g(1)*u_numerical(:,:,1) + 2*g(1)*v_numerical(:,:,1); %u v��ʷ��
    right2 = reshape(right2,(M-1)*(M-1),1);
              
    right3 = reshape(u_numerical(:,:,1),(M-1)*(M-1),1);
    right3 = -(1-sigma)*A*right3; %�ռ�ǰ��
       
    u_numerical_temp = left\( B*( right1 + right2 ) + right3 );    
    
    u_numerical(:,:,2) = reshape(u_numerical_temp,(M-1),(M-1));
    
    v_numerical(:,:,2) = 2/tau*(u_numerical(:,:,2) - u_numerical(:,:,1)) - v_numerical(:,:,1); %u v����滻
       
    for j = 3:N+1
        
       g = coe_a(tau,alpha,sigma,j-1);
       
       left = (2*sigma+1)/2/tau/sigma*g(j-1)*B + sigma*A; %B�������� A�նȾ���
       
       right1 = reshape(f(:,:,j-1),(M-1)*(M-1),1); %Դ��
              
       right2 = zeros(M-1,M-1);
       for jj = 2:j-1
           right2 =  right2 + v_numerical(:,:,jj)*(g(jj)-g(jj-1));%��������͵���
       end
       right2 = reshape(right2,(M-1)*(M-1),1);
       
       right3 = v_numerical(:,:,1)*g(1) + 2/tau*g(j-1)*u_numerical(:,:,j-1) ...
           - (2*sigma-1)/2/tau/sigma*g(j-1)*u_numerical(:,:,j-2) + (1-sigma)/sigma*g(j-1)*v_numerical(:,:,j-1);
       right3 = reshape(right3,(M-1)*(M-1),1); %���ױ����滻ʣ����
       
       right4 = reshape(u_numerical(:,:,j-1),(M-1)*(M-1),1);
       right4 = -(1-sigma)*A*right4; %�ռ�ǰ��
                                   
       u_numerical_temp = left\( B*( right1 + right2 + right3 ) + right4 );
             
       u_numerical(:,:,j) = reshape(u_numerical_temp,M-1,M-1);
       
       v_numerical(:,:,j) = 1/2/tau/sigma*((2*sigma+1)*u_numerical(:,:,j) - 4*sigma*u_numerical(:,:,j-1) ...
           + (2*sigma-1)*u_numerical(:,:,j-2)) - (1-sigma)/sigma*v_numerical(:,:,j-1); %u v����滻
    end   
           
end

u_end = u_numerical(:,:,end);

end

