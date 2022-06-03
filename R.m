%!! start-1 <-  !!!  start+1 -> !!
function res=R(Omega,start,endl)
    res(start)=1;
    if start>endl
        k=1; %start=n-1 endl=1
        for i=start-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    else
        k=-1; %start=2 endl=n
        for i=start-k:endl
            res(i)=Omega(i+k)*res(i+k);
        end
    end
end

