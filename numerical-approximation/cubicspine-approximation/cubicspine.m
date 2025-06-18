%code by jerry
function cubicspline(x,fx,M0,Mn)
n=length(x);M(1)=M0;M(n)=Mn;
for i=1:n-1
    h(i)=x(i+1)-x(i);
end
for i=1:n-2
    g(i)=h(i)/(h(i)+h(i+1));
end
A=zeros(n,n);A(1:n,1)=fx;
for j=1:2
    for i=1:n-j
        A(i,j+1)=(A(i+1,j)-A(i,j))/(x(i+j)-x(i));
    end
end
b(1)=6*A(1,3)-g(1)*M0;b(n-2)=6*A(n-2,3)-(1-g(n-2))*Mn;
b(2:n-3)=A(2:n-3);
M(2:n-1)=chase(g(2:n-2),2*ones(n-2,1),1-g(1:n-3),b);
syms y
for i=1:n-1
    al=A(i+1,1)/h(i)-h(i)*M(i+1)/6;be=A(i,1)/h(i)-h(i)*M(i)/6;
    f=M(i+1)*(y-x(i)).^3./(6*h(i))-M(i)*(y-x(i+1)).^3./(6*h(i))+al*(y-x(i))+be*(x(i+1)-y);
    fplot(f,[x(i),x(i+1)]);
    hold on;
end
plot(x,fx,'ro');
grid on;
hold off;
end
