clear
n=32;
x=0:2*pi/n:2*pi;

mtg

figure(1)
title('the solution 0n 8 x 8 x 8 grid points')
subplot(4,2,1)
surf(x,x,p1)
text(6,8,.3,'x = 2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,2)
surf(x,x,p2)
text(6,8,.3,'x = 2*2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,3)
surf(x,x,p3)
text(6,8,.3,'x = 3*2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,4)
surf(x,x,p4)
text(6,8,.0,'x = 4*2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,5)
surf(x,x,p5)
text(6,8,.3,'x = 5*2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,6)
surf(x,x,p6)
text(6,8,.3,'x = 6*2\pi/8')
xlabel('y')
ylabel('z')

subplot(4,2,7)
surf(x,x,p7)
text(6,8,.3,'x = 7*2\pi/8')
xlabel('y')
ylabel('z')



x1 = x(3)
for j = 1:n
  for m = 1:n
    b(j,m) = 0;
    for k=1:4
       b(j,m)=b(j,m)+sin(k*x1)*sin(k*x(j))*sin(k*x(m));
    end
  end
end
size(b)
