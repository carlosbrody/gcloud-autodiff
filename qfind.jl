"""
function y=qfind(x,ts)
    does a binary search: assumes that x is sorted low to high and unique.

x is a vector, ts is the target (can be one or many targets),
y is same length as ts  and returns the indices of the elements of x closest to each element of ts

COPIED FROM JEFF'S QFIND FUNCTION -- SEEMS TO HAVE A BUG IN THAT IT CAN RETURN -1?
"""

function qfind(x,ts)

    ys=Array(Int, 1, length(ts));
    for i=1:length(ts)
        t=ts[i];
        if isnan(t)
            y=NaN;
        else
            high = length(x);
            low = -1;
            if t>=x[end]
                y=length(x);
            else
                try
                    while (high - low > 1) 
                        probe = convert(Int, ceil((high + low) / 2));
                        if (x[probe] > t)
                            high = probe;
                        else
                            low = probe;
                        end
                    end
                    
                    y=low;
                catch 
                    y=low;
                end
            end
        end

        ys[i]=convert(Int, y);
    end

    return ys
end
