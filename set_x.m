function res=set_x(xi,ro)
    res=zeros(length(xi),1);
    for i=1:length(xi)
        res(i)=sqrt(xi(i)-ro^2);
    end
end
