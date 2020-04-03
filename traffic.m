clear all
close all

global dmin dmax vmax;
dmin = 15; % distance within which car stops (meters)
dmax = 100; % distance above which cars goes at max speed (meters)
vmax = 30; % maximum speed (meters/second)
L = 800; % length of one block (meters)
nb= 10; % number of blocks
radius = nb*L/(2*pi); % radius of ring (meters)
R = .075; % rate of arrival of cars at each entry point (1/second)
tmax = 10000; %time of simulation (seconds)
dt = 0.5; %time step (seconds)
clockmax=ceil(tmax/dt);

if(dt*vmax > dmin)
  error('dt*vmax > 0.1*dmin')
end
if(dt*R > 0.2)
  error('dt*R>0.1')
end

ncmax=ceil(2*tmax*R*nb); % max number of cars allowed (2x expected number)
ncwmax=100; % max length of queue allowed at each entry point

firstcar=zeros(1,nb);
lastcar=zeros(1,nb);
nextcar=zeros(1,ncmax);
onroad=zeros(1,ncmax);
x=zeros(1,ncmax);
xx=zeros(1,ncmax);
yy=zeros(1,ncmax);
tzero=zeros(1,ncmax);
tstart=zeros(1,ncmax);
texit=zeros(1,ncmax);
ncw=zeros(nb,1);
cw=zeros(nb,ncwmax);
xxw=zeros(1,nb*ncwmax);
yyw=zeros(1,nb*ncwmax);
tsave=zeros(1,clockmax);
vsave=zeros(1,clockmax);
nc_onroad=zeros(1,clockmax);
nc_waiting=zeros(1,clockmax);

hold off
figure(1)
set(gcf,'double','on')
h=plot(0,0,'r.');
axis(1.5*[-radius,radius,-radius,radius])
axis equal
axis manual
hold on
hw=plot(0,0,'b.');

nskip=5;


nc=0;
for clock = 1:clockmax
  t=clock*dt;

  %loop over blocks for arrival and entry of cars into ring:
  for b=1:nb
    bb=b-nb*(b==nb);ba=bb+1; %ba is index of block ahead of b
    
    %random arrival of cars:
    if(rand<R*dt)
      nc=nc+1;ncw(b)=ncw(b)+1; %one more car, and one more car waiting
      if(nc>ncmax)
        error('too many cars')
      end
      if(ncw(b)>ncwmax)
        error('too many cars waiting')
      end
      cw(b,ncw(b))=nc;         %put new car on queue for this block
      x(nc)=bb*L;bd(nc)=ceil(nb*rand); %locate car, assign exit destination
      tzero(nc)=t;nextcar(nc)=0;onroad(nc)=0; %note arrival time
    end
    
    %entry to ring when there is room to go:
    if(ncw(b)>0)&&((lastcar(ba)==0)||(L*bb+dmin<x(lastcar(ba))))
      c=cw(b,1); % index of car that is first in line
      cw(b,1:ncw(b))=[cw(b,2:ncw(b)),0]; %move up the whole line
      ncw(b)=ncw(b)-1; %one fewer car waiting
      tstart(c)=t;onroad(c)=1; %note time that car starts moving
      %insert car c into data structure of the block ba:
      if(lastcar(ba)==0) %if the block ba is empty, then
        firstcar(ba)=c; %car c is the first car on that block;
      else              %otherwise,
        nextcar(lastcar(ba))=c; %car c is behind the last car on that block
      end
      lastcar(ba)=c;nextcar(c)=0; %car c is the last car on block ba
    end
  end

  %loop over blocks for moving the cars:
  %note use of mod( ,nb*L), since road is a closed ring of length nb*L
  vsum=0; % initialize to find average velocity
  for b=1:nb
    if(firstcar(b)>0) %if the block is not empty:
      bb=b-nb*(b==nb);ba=bb+1; %ba is index of the block ahead
      %the first car is special, since it has to look at the block ahead
      c=firstcar(b);
      if(lastcar(ba)==0) %if the block ahead is empty, then
        d=dmax;          %set d to go at maximum speed;
      else               %otherwise set d to be
        d=mod(x(lastcar(ba))-x(c),nb*L); %distance to last car on next block
      end
      x(c)=mod(x(c)+dt*v(d),nb*L); %move car c at speed v(d)
      vsum=vsum+v(d);
      
      %now that first car has moved, move all other cars on block:
      ca=c; %index of car ahead of the car that will be moved next
      c=nextcar(c); %index of car that will be moved next
      while(c>0) %loop through all cars on block
        d=mod(x(ca)-x(c),nb*L); %distance to car ahead
        x(c)=mod(x(c)+dt*v(d),nb*L); %move car c at speed v(d)
        vsum=vsum+v(d);
        ca=c; %index of car ahead of the car that will be moved next
        c=nextcar(c); %index of car that will be moved next
      end
    end
  end
  vsave(clock)=vsum/sum(onroad);
  
  %loop over blocks to check whether first car has left the block:
  %if, so, it must be removed from the data structure of its former block,
  %and, if it has reached its destination, it should exit;
  %but if it has not reached its destination,
  %it should be written into the data structure of its new block.
  for b=1:nb
    if(firstcar(b)>0)
      bb=b-nb*(b==nb);ba=bb+1;
      c=firstcar(b);
      if((b<nb)&&(x(c)>=L*b))||((b==nb)&&(x(c)<L))
      %that is, if car c has reached or passed the end of block b, 
        firstcar(b)=nextcar(c); %then car behind it is now the first car,
        if(firstcar(b)==0) %but if there was no such car behind car c,
          lastcar(b)=0; %then block b is now empty
        end
        nextcar(c)=0; %there is now no car behind car c on its new block

        if(b==bd(c)) %if car c just passed is destination,
          texit(c)=t; %note exit time
          onroad(c)=0; %and remove car c from the road
        else %otherwise, write car c into data structure of block ba:
          if(lastcar(ba)==0) %if block ba was empty,
            firstcar(ba)=c; %then car c is now the first car on block ba
          else %otherwise, car c is now behind the previous last car:
            nextcar(lastcar(ba))=c; 
          end
          %either way, car c is now the last car on block ba:
          lastcar(ba)=c;
        end
      end
    end
  end
  tsave(clock)=t;
  nc_onroad(clock)=sum(onroad);
  nc_waiting(clock)=sum(ncw);
  %plot cars that are on the road, if any
  if (nc_onroad(clock)>0) && (mod(clock,nskip)==0)
    c_onroad=find(onroad);
    xx=radius*cos(x(c_onroad)/radius);
    yy=radius*sin(x(c_onroad)/radius);
    set(h,'xdata',xx,'ydata',yy)
    count=0;
    for b =1:nb
      theta = (b/nb)*2*pi;
      xxw(count+(1:ncw(b)))=(radius+dmin*(1:ncw(b)))*cos(theta);
      yyw(count+(1:ncw(b)))=(radius+dmin*(1:ncw(b)))*sin(theta);
      count=count+ncw(b);
    end
    if(count>0)
      set(hw,'xdata',xxw(1:count),'ydata',yyw(1:count))
    else
      set(hw,'xdata',10*radius,'ydata',10*radius) %plot offscreen
    end
    drawnow
  end
end

figure(2);plot(tsave,vsave)
figure(3);plot(tsave,nc_onroad)
figure(4);plot(tsave,nc_waiting)
c_exited = find(texit>0);
figure(5);plot(tzero(c_exited),texit(c_exited)-tzero(c_exited))


%filename v.m
function vv=v(d)
    global dmin dmax vmax;
    if(d<dmin)
      vv=0;
    elseif(d<dmax)
      vv=vmax*log(d/dmin)/log(dmax/dmin);
    else
      vv=vmax;
    end
end
