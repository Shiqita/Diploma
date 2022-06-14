
disp(df());

function res=v1(arg)
    LAM=3827.582;
    res=(w^2-(LAM)/(x^2-1)-(3)/(x^2-1)^2+(4+((LAM)/(2*w))^2)/(x^2-1)^2)^(1/2);
end
function res=df(arg)
    arg
    LAM=3827.582;
    res=((12*x)/(x^2-1)^3-(4*x*(4+LAM^2/(4*w^2)))/(x^2-1)^3+(2*LAM*x)/(x^2-1)^2)/(2*(w^2-LAM/(x^2-1)-(3)/(x^2-1)^2+(4+LAM^2/(4*w^2))/(x^2-1)^2)^(1/2));
end