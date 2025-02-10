% �����׺��򲨷���

function [w_ini] = backwards(fractional_order,u_end_u_delta)
alpha = fractional_order;
T = 1;
L = pi;
M = 80; % �ռ�������Ŀ
h = L/M; % �ռ䲽��
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

z = [1000]; % ʱ������
for i = 1:length(z)
   N = z(i); % ʱ��������Ŀ
   tau = T/N; % ʱ�䲽�� 
   t = 0:tau:T; % ʱ������
  
   f = zeros(M-1,M-1,N+1);

   
   W = zeros(M-1,M-1,N+1); % ��ֵ��
   W(:,:,N+1) = zeros(M-1,M-1,1); 
   W_der_end = u_end_u_delta;
   
   U = zeros(M-1,M-1,N+1); % ��ֵ��
   U(:,:,N+1) = zeros(M-1,M-1,1); % ��ֵ��
   
   l = 0:N-1;
   a = tau^(2-alpha)/gamma(3-alpha)*( (l+1).^(2-alpha) - l.^(2-alpha) );
   
   W(:,:,N) = W(:,:,N+1) - W_der_end*tau;
   U(:,:,N) = 1/a(1)*W(:,:,N);
   
   for j = 1:N-1
       left = a(1)*B+tau^2/4*A;
       right1_pre = zeros(M-1,M-1);
       for jj = 2:j+1
           right1_pre =  right1_pre + U(:,:,N-j+jj-1)*a(jj); %��������͵���
       end
       right1 = 2*W(:,:,N-j+1) - W(:,:,N-j+2) - right1_pre + tau^2*f(:,:,N-j+1);
       right1 = reshape(right1,(M-1)*(M-1),1);       
       right2 = tau^2/2*U(:,:,N-j+1) + tau^2/4*U(:,:,N-j+2);
       right2 = reshape(right2,(M-1)*(M-1),1);
       u_numerical_temp = left\( B*right1 - A*right2 );
       U(:,:,N-j) = reshape(u_numerical_temp,(M-1),(M-1));
       right3 = zeros(M-1,M-1);
       for jj = 1:j+1
           right3 =  right3 + U(:,:,N-j+jj-1)*a(jj); %��������͵���
       end
       W(:,:,N-j) = right3;
   end
   
   
end

w_ini = W(:,:,1);

end






