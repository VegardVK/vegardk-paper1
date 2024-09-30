using Combinatorics
using Printf
using PrettyTables

function get_elevation_matrix(n, m)
    M_n = get_M(n)
    M_m = get_M(m)

    M_m_inv = inv(M_m)
    j_array = j_nm(n, m)

    elevation_matrix = round.(M_n * j_array * M_m_inv, digits=2)
    return elevation_matrix
end

function get_M(n)
    m_array = zeros(n+1, n+1)
    for i in 0:n
        for j in 0:(n-i)
            m_array[i+1, i+j+1] = m_ij(n, i, j)
        end
    end
    return m_array
    # println(m_array)
end

function m_ij(n, i, j)
    (-1)^j * binomial(n, i) * binomial(n-i, j)
end

function j_nm(n,m)
    j_array = zeros(n+1, m+1)
    for i in 0:n
        j_array[i+1, i+1] = 1
    end
    return j_array
end


function get_k_matrix(n)
    k_array = zeros(n+1, n)
    for i in 1:n
        k_array[i, i] = -n
        k_array[i+1, i]  = n
    end
    return k_array
end

# n = 2
# # m = 3
# # matrix = get_elevation_matrix(n, m)
# matrix = get_k_matrix(n)
# pretty_table(matrix)
