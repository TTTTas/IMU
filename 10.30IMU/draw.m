figure('WindowState','maximized')
subplot(3,1,1)
hold on
plot(IMU_result(:,1),IMU_result(:,2),'-','LineWidth',2,'DisplayName','Pure-IMU_{B}')
plot(real(:,1),real(:,3),'-.','LineWidth',2,'DisplayName','Real_{B}')
title("纬度差异")
xlabel("Time/s",'FontName','Times New Roman')
ylabel("B/deg",'FontName','Times New Roman')
legend
grid on
subplot(3,1,2)
hold on
plot(IMU_result(:,1),IMU_result(:,3),'-','LineWidth',2,'DisplayName','Pure-IMU_{L}')
plot(real(:,1),real(:,2),'-.','LineWidth',2,'DisplayName','Real_{L}')
title("经度差异")
xlabel("Time/s",'FontName','Times New Roman')
ylabel("L/deg",'FontName','Times New Roman')
legend
grid on
subplot(3,1,3)
hold on
plot(IMU_result(:,1),IMU_result(:,4),'-','LineWidth',2,'DisplayName','Pure-IMU_{H}')
plot(real(:,1),real(:,4),'-.','LineWidth',2,'DisplayName','Real_{H}')
title("高程差异")
xlabel("Time/s",'FontName','Times New Roman')
ylabel("H/m",'FontName','Times New Roman')
legend
grid on