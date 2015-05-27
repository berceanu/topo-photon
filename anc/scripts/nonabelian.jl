import HH

#reload("HH")

const q = 11

H = zeros(Complex{Float64}, (q,q))
HH.momsphhmat!(H, 1., 1.)

ishermitian(H)
#H == H'
#H == ctranspose(H)

function chop(A, eps)
    for ix in eachindex(A)
        if A[ix] < eps; A[ix] = 0
        end
    end
end
