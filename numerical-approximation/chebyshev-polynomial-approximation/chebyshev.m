function result=chebyshev(func,a,b,n)
if ischar(func)||isstring(func)
    %首先，判断输入的函数是字符串还是句柄。如果是字符串就转化成句柄；如果是句柄就保持不变。
    f1=str2func(['@(x)' func]);
elseif isa(func,'function_handle')
    f1=func;
else
    error('输入错误');
end
%将函数从[a,b]区间转移到[-1,1]区间
f=@(t) f1((b-a)/2*t+(a+b)/2);
%现在，可以计算多项式系数，并将它们存储在向量c里。然后，再赋值给函数p。
c=zeros(n+1,1);
for k=1:n+1
    c(k)=2/pi*integral(@(x) f(x).*chebyshevT(k-1,x)./sqrt(1-x.^2),-1,1);
end
p=@(x) c(1)/2;
for i=1:n
    p=@(x) p(x)+c(i+1)*chebyshevT(i,x);
end
%将函数转换回原区间并画图
figure;
fplot(@(x) f((x-(b+a)/2)/(b-a)*2),[a,b]);
hold on;
fplot(@(x) p((x-(b+a)/2)/(b-a)*2),[a,b]);
legend('f(x)',['P',num2str(n),'(x)']);
grid on;
%计算加权均方误差
result=['加权均方误差为' num2str(sqrt(integral(@(x) (f(x)-p(x)).^2,-1,1)))];
end
