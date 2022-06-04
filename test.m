xi=(0:0.5:100);

ro=sqrt(2);

X=set_x(xi,ro); % poly4aem iz XI

disp(X)

function res=set_x(xi,ro)
    res=zeros(length(xi),1);
    for i=1:length(xi)
        res(i)=abs(sqrt(xi(i)-ro^2));
    end
end