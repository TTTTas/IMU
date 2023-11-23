% %% 读取示例文件
% f1=fopen("IMU.bin");
% [f,N]=fread(f1,'double');
% f=reshape(f,7,N/7)';
% fclose(f1);
% 
% Len=N/7;
% l=0;
% for i=1:Len
%     if f(i,1)>91620
%         l=i; 
%         break;
%     end
% end
% stat=f(l:end,:);
% time=stat(:,1)-91620;
% 
% %初值设置
% rou=180/pi;
% we=7.292115e-5;
% B0=23.1373950708/rou;           %初始纬度
% L0=113.3713651222/rou;          %初始经度
% H0=2.175;                       %初始高度
% phi0=0.0107951084511778/rou;    %横滚
% theta0=-2.14251290749072/rou;   %俯仰
% pusi0=-75.7498049314083/rou;    %航向
% dt=0.005;                       %时间历元间隔
% Cnb0=Eu2Dcm(phi0,theta0,pusi0); %初始的方向余弦矩阵
% v0=[0;0;0];                     %初始速度矩阵
% 
% B1=B0;L1=L0;H1=H0;v1=v0;Cnbk1=Cnb0;

%% 读取小推车采集数据
stat=importdata("IMU.txt");
time=stat(:,1)-274119.9986884708;

% 初值设置
rou=180/pi;
we=7.292115e-5;
B0=30.5278036089/rou;           %初始纬度
L0=114.3557909607/rou;          %初始经度
H0=20.9960000000000;                      %初始高度
%[phi0,theta0,pusi0]=ssi_alignment(stat,we,B0,H0);  %初始对准计算姿态角static-state-initial-alignment
phi0=0.33539381/rou;
theta0=0.68297002/rou;
pusi0=358.04142691/rou;

dt=0.005;                       %时间历元间隔
Cnb0=Eu2Dcm(phi0,theta0,pusi0); %初始的方向余弦矩阵
v0=[0;0;0];

B1=B0;L1=L0;H1=H0;v1=v0;Cnbk1=Cnb0;

%% 逐历元解算
fout=fopen("Result1.txt",'wt');
for i=60715:301000
    
    %该历元角增量、速度增量
    tk=[stat(i,2);stat(i,3);stat(i,4)];
    vk=[stat(i,5);stat(i,6);stat(i,7)];

    %前一历元角增量、速度增量
    if(i==60715)
        tk1=[0;0;0];vk1=[0;0;0];
    else
        tk1=[stat(i-1,2);stat(i-1,3);stat(i-1,4)];
        vk1=[stat(i-1,5);stat(i-1,6);stat(i-1,7)];
    end

   
    %姿态更新
    phik=tk+cross(tk1,tk)/12;
    [RM,RN]=CalRadium(B1);
    wie1=[we*cos(B1);0;-we*sin(B1)];
    wen1=[v1(2)/(RN+H1);-v1(1)/(RM+H1);-v1(2)*tan(B1)/(RN+H1)];
    ksik=(wie1+wen1)*dt;
    phik1=MCross(phik);
    ksik1=MCross(ksik);
    if(norm(phik)==0)
        Cbb=eye(3);
    else
        Cbb=eye(3)+sin(norm(phik))/norm(phik)*phik1+(1-cos(norm(phik)))/norm(phik)^2*phik1^2;
    end
    if(norm(ksik)==0)
        Cnn=eye(3);
    else
        Cnn=eye(3)-sin(norm(ksik))/norm(ksik)*ksik1+(1-cos(norm(ksik)))/norm(ksik)^2*ksik1^2;
    end
    Cnbk=Cnn*Cnbk1*Cbb;
    [phitk,thetak,pusik]=Dcm2Eu(Cnbk);

    %线性外推
    if(i==60715)
        vk12=v1;
        B12=B1;H12=H1;
        wie12=wie1;wen12=wen1;
    else
        vk12=1.5*v1-0.5*v2;
        B12=1.5*B1-0.5*B2;
        H12=1.5*H1-0.5*H2;
        wie12=1.5*wie1-0.5*wie2;
        wen12=1.5*wen1-0.5*wen2;
    end

    if((109878<=i&&i<=123458)||(173379<=i&&i<=188172)||(237250<=i&&i<=250177)||i>299308)
            v=[0;0;0];
    else
    %速度更新
    g=gCal(B12,H12);
    ksinn=(wie12+wen12)*dt;
    dvfb=vk+0.5*cross(tk,vk)+(cross(tk1,vk)+cross(vk1,tk))/12;
    dvgk=(g-cross((2*wie12+wen12),vk12))*dt;
    dvfk=(eye(3)-0.5*MCross(ksinn))*Cnbk1*dvfb;
    v=v1+dvfk+dvgk;
    end
    
    
    %位置更新
    H=H1-0.5*(v1(3)+v(3))*dt;
    B=B1+(v(1)+v1(1))/(2*RM+H+H1)*dt;
    [RM12,RN12]=CalRadium(0.5*(B+B1));
    L=L1+(v(2)+v1(2))/(2*RN12+H+H1)/cos(0.5*(B+B1))*dt;
    

    %输出到文件
    output=[time(i),B*rou,L*rou,H,v',phitk,thetak,pusik];
    fprintf(fout,'%4f %9.8f %9.8f %5.3f %8.7f %8.7f %8.7f %8.7f %8.7f %8.7f\n',output(1),output(2),output(3),output(4),output(5),output(6),output(7),output(8),output(9),output(10));

    %更新迭代值
    v2=v1;v1=v;
    B2=B1;B1=B;
    L1=L;
    H2=H1;H1=H;
    Cnbk1=Cnbk;
    wie2=wie1;wen2=wen1;

end
fclose(fout);

% %% 示例数据参考值、计算值读取
% 
% %读入示例数据真值
% f2=fopen("PureINS.bin");
% [Ref,NN]=fread(f2,'double');
% Ref=reshape(Ref,10,NN/10)';
% fclose(f2);
% 
% %读入计算结果
% cal=importdata('Result.txt');
% 
% % 结果作图
% %轨迹图
% figure(1);
% geoplot(cal(:,2),cal(:,3));
% title('Path');
% figure(2);
% subplot(1,3,1);
% plot(cal(:,1),cal(:,4));
% xlabel('time/s');ylabel('H/m');title('H');
% subplot(1,3,2);
% plot(cal(:,1),cal(:,5),cal(:,1),cal(:,6),cal(:,1),cal(:,7));
% xlabel('time/s');ylabel('v m/s');legend('vn','ve','vd');title('Velocity');
% subplot(1,3,3);
% plot(cal(:,1),cal(:,8),cal(:,1),cal(:,9),cal(:,1),cal(:,10));
% xlabel('time/s');ylabel('φ-θ-ψ /deg');legend('roll','pitch','yaw');title('Posture');
% figure(3);
% res=cal-Ref;
% subplot(2,2,1)
% plot(time,resP,'lineWidth',1);legend('B','L');
% xlabel('time/s');ylabel('Lon-Lat residual/deg');title('B-L Bias');
% subplot(2,2,2)
% plot(time,resH,'lineWidth',1);legend('H');title('Height Bias');
% xlabel('time/s');ylabel('Height residual/m');title('H Bias');
% subplot(2,2,3)
% plot(time,resV,'lineWidth',1);legend('vn','ve','vd');
% xlabel('time/s');ylabel('Velocity residual/(m/s)');title('N-E-D·Velocity Bias');
% subplot(2,2,4)
% plot(time,resC,'lineWidth',1);legend('roll','pitch','yaw');
% xlabel('time/s');ylabel('Posture residual/deg');title('φ-θ-ψ Bias');

%% 小推车数据参考值、计算值读取

Ref=importdata("refpose.txt");
cal=importdata("Result1.txt");
Ref(:,1)=Ref(:,1)-274119.999;
ref=Ref(61016:302502,:);

% 处理参考数据
for i=2:240286
    if ref(i,1)-ref(i-1,1) <0.002
        ref(i,:)=[];
    end
end

%% 处理航向角取值
for i=1:240286
    if(ref(i,10)>180)
        ref(i,10)=ref(i,10)-360.0;
    end
end
%%  将经纬度转换到站心ENU坐标系
for  i=1:length(cal)
    [Rm,Rn]=CalRadium(cal(i,2)/rou);
    cal(i,11)=(cal(i,2)/rou-B0)*(Rm+cal(i,4));
    cal(i,12)=(cal(i,3)/rou-L0)*(Rn+cal(i,4))*cos(cal(i,2)/rou);
end

for  i=1:length(ref)
    [Rm,Rn]=CalRadium(ref(i,5)/rou);
    ref(i,41)=(ref(i,5)/rou-B0)*(Rm+ref(i,7));
    ref(i,42)=(ref(i,6)/rou-L0)*(Rn+ref(i,7))*cos(ref(i,5)/rou);
end

%%
plot(cal(:,3),cal(:,2),ref(:,6),ref(:,5));
%%  小车数据绘图
%轨迹\高程图
figure(1);
plot(cal(:,12),cal(:,11),'.-b',ref(:,42),ref(:,41),'.-r');legend('Cal','Ref');
xlabel('E(m)');ylabel('N(m)');title('Trace');
figure(2);
subplot(1,3,1);
plot(cal(:,1),cal(:,4));
xlabel('time(s)');ylabel('H(m)');title('H');
subplot(1,3,2);
plot(cal(:,1),cal(:,5),cal(:,1),cal(:,6),cal(:,1),cal(:,7));
xlabel('time(s)');ylabel('v(m/s)');legend('vn','ve','vd');title('Velocity');
subplot(1,3,3);
plot(cal(:,1),cal(:,8),cal(:,1),cal(:,9),cal(:,1),cal(:,10));
xlabel('time(s)');ylabel('φ-θ-ψ(deg)');legend('roll','pitch','yaw');title('Posture');
figure(3);
plot(cal(:,1),cal(:,4),'.-b',ref(:,1),ref(:,7),'.-r');legend('Cal','Ref');
xlabel('time(s)');ylabel('H(m)');title('H-change');
%误差序列图
figure(4);
resP=cal(:,2:3)-ref(:,5:6);resH=cal(:,4)-ref(:,7);resV=cal(:,5:7)-ref(:,18:20);resC=cal(:,8:10)-ref(:,8:10);%经纬度、高程、速度、姿态角
subplot(2,2,1)
plot(cal(:,1),resP,'lineWidth',1);legend('B','L');
xlabel('time(s)');ylabel('Lon-Lat residual(deg)');title('B-L Bias');
subplot(2,2,2)
plot(cal(:,1),resH,'lineWidth',1);legend('H');title('Height Bias');
xlabel('time(s)');ylabel('Height residual(m)');title('H Bias');
subplot(2,2,3)
plot(cal(:,1),resV,'lineWidth',1);legend('vn','ve','vd');
xlabel('time(s)');ylabel('Velocity residual/(m/s)');title('N-E-D·Velocity Bias');
subplot(2,2,4)
plot(cal(:,1),resC,'lineWidth',1);legend('roll','pitch','yaw');
xlabel('time(s)');ylabel('Posture residual(deg)');title('φ-θ-ψ Bias');
%% function

%计算子午圈和卯酉圈半径函数
function [RM,RN]=CalRadium(B0)
    a=6378137.0;
    e=0.08181919104;
    RM=a*(1-e*e)/power(1-e*e*sin(B0)^2,3/2);
    RN=a/sqrt(1-e*e*sin(B0)^2);
end

%欧拉角转姿态矩阵
function Cnb=Eu2Dcm(phi,theta,pusi)
    ct=cos(theta);cu=cos(pusi);cp=cos(phi);
    st=sin(theta);su=sin(pusi);sp=sin(phi);
    Cnb=[ct*cu,-cp*su+sp*st*cu,sp*su+cp*st*cu;
        ct*su,cp*cu+sp*st*su,-sp*cu+cp*st*su;
        -st,sp*ct,cp*ct];
end

%姿态矩阵转欧拉角
function [phi,theta,pusi]=Dcm2Eu(Cnb)
    phi=atan2d(Cnb(3,2),Cnb(3,3));
    theta=atand(-Cnb(3,1)/sqrt(Cnb(3,2)^2+Cnb(3,3)^2));
    pusi=atan2d(Cnb(2,1),Cnb(1,1));
end

%反对称阵
function Mc=MCross(m)
    Mc=[0,-m(3),m(2);
        m(3),0,-m(1);
        -m(2),m(1),0];
end

%计算正常重力
function g=gCal(B,H)
    g0=9.7803267715*(1+0.0052790414*sin(B)^2+0.0000232718*sin(B)^4);
    gbh=g0-(3.087691089e-6-4.397731e-9*sin(B)^2)*H+0.721e-12*H^2;
    g=[0;0;gbh];
end

%计算初始姿态角static-state-initial-alignment
function [phi,theta,pusi]=ssi_alignment(stat,Ome,lat,H)
    gn=gCal(lat,H);
    ome_ie_n=[Ome*cos(lat);0;-Ome*sin(lat)];
    vg=gn/norm(gn);
    vw=cross(gn,ome_ie_n)/norm(cross(gn,ome_ie_n));
    vgw=(cross(cross(gn,ome_ie_n),gn))/norm(cross(cross(gn,ome_ie_n),gn));

    SumS=sum(stat(1:60714,:))/60714;
    gb0=[-SumS(1,5);-SumS(1,6);-SumS(1,7)];
    ome_ie_b0=[SumS(1,2);SumS(1,3);SumS(1,4)];
    wg=gb0/norm(gb0);
    ww=(cross(gb0,ome_ie_b0))/norm(cross(gb0,ome_ie_b0));
    wgw=(cross(cross(gb0,ome_ie_b0),gb0))/norm(cross(cross(gb0,ome_ie_b0),gb0));

    V=[vg,vw,vgw];
    W=[wg.';ww.';wgw.'];
    Cb_n=V*W;
    phi=atan(-Cb_n(3,1)/sqrt(Cb_n(3,2)^2+Cb_n(3,3)^2));
    theta=atan2(Cb_n(3,2),Cb_n(3,3));
    pusi=atan2(Cb_n(2,1),Cb_n(1,1));
end