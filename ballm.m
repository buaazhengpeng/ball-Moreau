clc;
m=1;
R1=1;
R2=0.2;
g=9.8;
u=0.3;
Jc=1/2*m*R2^2;
M=[m,0,0;0,m,0;0,0,Jc];
R=[0;m*g;0];
h=0.5e-3;
TT=3;
r=1;
tolerance=10^(-10);
n=floor(TT/h)+1;
t=zeros(1,n);
for i=1:n
    t(1,i)=(i-1)*h;
end
qd=zeros(3,n);
qd(1,1)=R2-R1;
qv=zeros(3,n);
qa=zeros(3,n);
g_T=zeros(1,n);%相对切向加速度
g_N=zeros(1,n);%法向切向速度
qa(:,1)=inv(M)*R;
qn=zeros(1,n);
qf=zeros(1,n);
qh=[0;m*g;0];%广义力
for i=2:n%i-1步到i步
    %计算中间步
    tt=t(1,i-1)+0.5*h;
    qdd=qd(:,i-1)+0.5*h.*qv(:,i-1);
    iM=inv(M);
    %迭代变量的
    qn(:,i)=0.0;%qn(:,i-1);
    qf(:,i)=0.0;%qf(:,i-1);
    qvv=zeros(3,1);
    %迭代过程中不变的量W_T,W_N
    W_N=[-qdd(1,1)/sqrt(qdd(1,1)^2+qdd(2,1)^2);...
         -qdd(2,1)/sqrt(qdd(1,1)^2+qdd(2,1)^2);0];
    W_T=[qdd(2,1)/sqrt(qdd(1,1)^2+qdd(2,1)^2);...
         -qdd(1,1)/sqrt(qdd(1,1)^2+qdd(2,1)^2);R2];
    error=1.0;
    while (error>tolerance)
        %更新速度
        qvv=qv(:,i-1)+iM*(qh.*h+W_N*qn(:,i)+W_T*qf(:,i));
        g_N(:,i)=W_N'*qvv;
        g_T(:,i)=W_T'*qvv;
        Projc=projc(qf(1,i)-r.*g_T(1,i),qn(1,i),u);
        dProjc=dprojc(qf(1,i)-r.*g_T(1,i),qn(1,i),u);
        Projn=projn(qn(1,i)-r.*g_N(1,i));
        dProjn=dprojn(qn(1,i)-r.*g_N(1,i));
        H1=qn(:,i)-Projn;
        H2=qf(:,i)-Projc;
        H=[H1;H2];
        D=[1-dProjn*(1-r.*W_N'*iM*W_N),r.*dProjn*W_N'*iM*W_T;...
            r.*dProjc*W_T'*iM*W_N,1-dProjc*(1-r.*W_T'*iM*W_T)];
        as=zeros(2,1);
        as=[qn(:,i);qf(:,i)]-inv(D)*H;
        qn(:,i)=as(1);
        qf(:,i)=as(2);
        error=max(abs([H1;H2]));
    end
    qv(:,i)=qv(:,i-1)+iM*(qh.*h+W_N*qn(:,i)+W_T*qf(:,i));
    qd(:,i)=qdd+0.5*h.*qv(:,i);
    disp(t(1,i));
end

%计算误差
% errd1=0;errd2=0;
% errv1=0;errv2=0;
% erra1=0;erra2=0;mm=floor(h/1.0e-6);
% for i=1:n
%     ii=mm*(i-1)+1;
%     errd1=errd1+abs(qd(1,i)-rqd(1,ii));
%     errd2=errd2+abs(rqd(1,ii));
%     errv1=errv1+abs(qv(1,i)-rqv(1,ii));
%     errv2=errv2+abs(rqv(1,ii));
% %     erra1=erra1+abs(qa(i,1)-rqa(1,mm*(i-1)+1));
% %     erra2=erra2+abs(rqa(1,mm*(i-1)+1));
% end
% errd=errd1/errd2
% errv=errv1/errv2
% erra=erra1/erra2
%----------------------画图----------------------
figure;
grid on
hold on 
plot(t,qv(1,:),'r-','linewidth',2)   
plot(t,qf(1,:)./h,'b-','linewidth',2)
plot(t,qv(1,:),'r-','linewidth',2,'MarkerIndices',round(linspace(1,length(t),100))); % 接触点1相对切向速度      
xlabel('{\itt} [s]'); ylabel('{\itg_T} [m/s]');
% figure;
% grid on
% hold on 
% plot(t,qd(1,:),'k-','linewidth',2,'MarkerIndices',round(linspace(1,length(t),100))); % 接触点1位移     
% xlabel('{\itt} [s]'); ylabel('{\itx} [m/s]');
% figure;
% grid on
% hold on 
% plot(t,qf(1,:)./h,'b-','linewidth',2,'MarkerIndices',round(linspace(1,length(t),100))); % 接触点1摩擦力      
% xlabel('{\itt} [s]'); ylabel('{\itF_f} [m/s]');
% figure;
% grid on
% hold on 
% plot(t,qa(1,:),'k-','linewidth',2,'MarkerIndices',round(linspace(1,length(t),100))); % 接触点1加速度     
% xlabel('{\itt} [s]'); ylabel('{\ita} [m/s]');