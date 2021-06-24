clc
clear 
% Ut=Uxx+Uyy
%Explicit Method
%Ut=(Un+1-Un)/dt, Uxx=(Ui+1-2*Ui+Ui-1)/dx^2, Uyy=Uj+1-2*Uj+Uj-1)/dy^2
%dx=dy, h=dt/dx^2<=0.25
n=10;
y=linspace(0,1,n);
x=linspace(0,1,n);

dx=length(x)/n;

dt=0.1;
h=dt/dx^2;
T=200;
u=zeros(length(x),length(y),T/dt);
u(:,:,1)=0;
u(1,:,:)=0; u(length(x),:,:)=25; u(:,1,:)=10;

 

for k=1:T/dt-1
%     while Error>1e-3
%     z=z+1;
%     uold=u;

    for j=2:length(y)-1
        for i=2:length(x)-1
            u(i,j,k+1) =(1-4*h)*u(i,j,k)+h*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k));
            u(:,length(y),:)=u(:,length(y)-1,:);
        end
    end
%        Error= max(abs(uold-u)); 
end

% end 


Tx200 = u(:, : ,200/dt );
Tx50 = u(:, : ,50/dt);
Tx5 = u(:, : ,5/dt);
Tx2 = u(:, : ,2/dt);
Tx1 = u(:, : ,1/dt);


figure (1)
subplot(2, 1, 1);
contourf(Tx200)
title('Explicit Method t= 200s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (2)
subplot(2, 1, 1);
contourf(Tx50)
title('Explicit Method t= 50s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (3)
subplot(2, 1, 1);
contourf(Tx5)
title('Explicit Method t= 5s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (4)
subplot(2, 1, 1);
contourf(Tx2)
title('Explicit Method t= 2s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (5)
subplot(2, 1, 1);
contourf(Tx1)
title('Explicit Method t= 1s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

clear 
% Ut=Uxx+Uyy
%implicit Method
%Ut=(Un-Un-1)/dt, Uxx=(Ui+1-2*Ui+Ui-1)/dx^2, Uyy=Uj+1-2*Uj+Uj-1)/dy^2
%dx=dy, h=dt/dx^2<=0.25
n=10;
y=linspace(0,1,n);
x=linspace(0,1,n);
dx=length(x)/n;

disp('********************************************************')

dt=0.1;
h=dt/dx^2
T=200;
u=zeros(length(x),length(y),T/dt);
u(:,:,1)=0;
u(1,:,:)=0; u(length(x),:,:)=25; u(:,1,:)=10;

z=0;

for k=2:T/dt
    Error =1;
    while Error>1e-3
        z=z+1;
        uold(:, :, k)=u(:, :, k);

        for j=2:length(y)-1 
            for i=2:length(x)-1
                u(i,j,k) =u(i,j,k-1)/(1+4*h)+((h/(1+4*h))*(u(i+1,j,k)+u(i-1,j,k)+u(i,j+1,k)+u(i,j-1,k)));
                u(:,length(y),:)=u(:,length(y)-1,:);
            end
        end
        Error= max(max(abs(uold(:, :, k)- u(:, :, k))));
    end 
end
Tx200 = u(:, : ,200/dt );
Tx50 = u(:, : ,50/dt);
Tx5 = u(:, : ,5/dt);
Tx2 = u(:, : ,2/dt);
Tx1 = u(:, : ,1/dt);

figure (1)
subplot(2, 1, 2);
contourf(Tx200)
title('Implicit Method t= 200s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (2)
subplot(2, 1, 2);
contourf(Tx50)
title('Implicit Method t= 50s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (3)
subplot(2, 1, 2);
contourf(Tx5)
title('Implicit Method t= 5s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (4)
subplot(2, 1, 2);
contourf(Tx2)
title('Implicit Method t= 2s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])

figure (5)
subplot(2, 1, 2);
contourf(Tx1)
title('Implicit Method t= 1s')
set(gca,'Xticklabel',[])
set(gca,'Yticklabel',[])
