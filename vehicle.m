%
% vehicle.m
%
% created on: 2023
%     author: Daniel A.
%
% see readme file for more information on the vehicle example
%
% you need to run ./vehicle binary first 
%
% so that the files: vehicle_ss.bdd 
%                    vehicle_obst.bdd
%                    vehicle_target.bdd
%                    vehicle_controller.bdd 
% are created
%

% function vehicle
clear set
close all



%% simulation

% target set
lb=[9 0];
ub=lb+0.5;
v=[9 0; 9.51  0; 9 0.51; 9.51 .51];
% initial state
x0=[0, 0.2, 0];

%data driven controller printed in black
controller=SymbolicSet('vehicle_controller.bdd','projection',[1 2 3]);

%model based controller printed in red
controller2=SymbolicSet('vehicleModel_controller.bdd','projection',[1 2 3]);
%target region
target=SymbolicSet('vehicle_target.bdd');
target2=SymbolicSet('vehicleModel_target.bdd');

y=x0;
v=[];
y2=x0;
v2=[];
uu2=[];
%plot path for data driven control vehicle
  while(1)        
      
      if (target.isElement(y(end,:))) 
        break;
      end 
    
      u=controller.getInputs(y(end,:));
      v=[v; u(1,:)];
      [t x]=ode45(@unicycle_ode,[0 .25], y(end,:), [],u(1,:));
    
      y=[y; x(end,:)];
   end
   disp('DD complete')
%plot path for model based control vehicle
disp('start')
while(1)
  
  if (target2.isElement(y2(end,:)))
    break;
  end 
  
  u2= controller2.getInputs(y2(end,:));
%   disp(u2);
  v2=[v2; u2(1,:)];
  [t2 x2]=ode45(@unicycle_ode,[0 .25], y2(end,:), [],u2(1,:));
  
  y2=[y2; x2(end,:)];
  uu2=[uu2; u2(1,:)];
end
% profile on
disp('MB complete')
% save('control.mat', uu2)


%% plot the vehicle domain
% colors
colors=get(groot,'DefaultAxesColorOrder');


% load the symbolic set containig the abstract state space
set=SymbolicSet('vehicle_ss.bdd','projection',[1 2]);
plotCells(set,'facecolor','none','edgec',[0.8 0.8 0.8],'linew',.1)
hold on

% load the symbolic set containig obstacles
set=SymbolicSet('vehicle_obst.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(1,:)*0.5+0.1,'edgec',colors(1,:),'linew',.1)

% plot the real obstacles and target set
plot_domain

% load the symbolic set containig target set
set=SymbolicSet('vehicle_target.bdd','projection',[1 2]);
plotCells(set,'facecolor',colors(2,:)*0.5+0.1,'edgec',colors(2,:),'linew',.1)

% plot initial state  and trajectory
plot(y(:,1),y(:,2),'k.-','LineWidth',1.0) %plot data driven controller
plot(y2(:,1),y2(:,2),'r.-','LineWidth',1.0) %plot model based controller
plot(y(1,1),y(1,2),'.','color',colors(5,:),'markersize',20)

box on
axis([-.5 10.5 -.5 10.5])
axis tight
ylabel('$y(m)$','Interpreter','latex', 'FontSize',18)
xlabel('$x(m)$','Interpreter','latex', 'FontSize',18)

% profile off
% profile viewer

% end

function dxdt = unicycle_ode(t,x,u)
  w = [0.0 0.0 0.0];
  dxdt = zeros(3,1);
  c=atan(tan(u(2))/2);
  %c=0;

  dxdt(1)=u(1)*cos(c+x(3))/cos(c)+w(1);
  dxdt(2)=u(1)*sin(c+x(3))/cos(c)+w(2);
  dxdt(3)=u(1)*tan(u(2))+w(3);


end


function plot_domain

colors=get(groot,'DefaultAxesColorOrder');

v=[9 0; 9.5  0; 9 0.5; 9.5 .5];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(2,:),'edgec',colors(2,:));


v=[1     0  ;1.2  0   ; 1     9    ; 1.2 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   0  ;2.4  0   ; 2.2   5    ; 2.4 5   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[2.2   6  ;2.4  6   ; 2.2   10   ; 2.4 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[3.4   0  ;3.6  0   ; 3.4   9    ; 3.6 9   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[4.6   1  ;4.8  1   ; 4.6   10   ; 4.8 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   0  ;6    0   ; 5.8   6    ; 6   6   ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[5.8   7  ;6    7   ; 5.8   10   ; 6   10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[7     1  ;7.2  1   ; 7     10   ; 7.2 10  ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.2   0  ;8.4  0   ; 8.2   8.5  ; 8.4 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   8.3;9.3  8.3 ; 8.4   8.5  ; 9.3 8.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   7.1;10   7.1 ; 9.3   7.3  ; 10  7.3 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   5.9;9.3  5.9 ; 8.4   6.1  ; 9.3 6.1 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   4.7;10   4.7 ; 9.3   4.9  ; 10  4.9 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[8.4   3.5;9.3  3.5 ; 8.4   3.7  ; 9.3 3.7 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));
v=[9.3   2.3;10   2.3 ; 9.3   2.5  ; 10  2.5 ];
patch('vertices',v,'faces',[1 2 4 3],'facec',colors(1,:),'edgec',colors(1,:));


end
%%}
