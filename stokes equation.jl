module ApproxOperator

include("elements.jl")  # 確保有這行
include("stokes.jl")    # 確保有這行

end  # module ApproxOperator

using ApproxOperator.AbstractElement

#===== 粘性项算子：μ∫∇u:∇v dΩ → 对应矩阵 A =====#
function (op::Operator{:∫∫μ∇u∇vdxdy})(aᵤ::T; k::AbstractMatrix{Float64}) where T<:AbstractElement
    𝓒 = aᵤ.𝓒; 𝓖 = aᵤ.𝓖
    μ = op.μ
    for ξ in 𝓖
        B₁ = ξ[:∂𝝭∂x]  # 速度形函数 x 导数
        B₂ = ξ[:∂𝝭∂y]  # 速度形函数 y 导数
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼   # 速度自由度全局索引（每个节点2自由度）
            for (j,xⱼ) in enumerate(𝓒)
                J = xⱼ.𝐼
                # 粘性项贡献：μ ∫ (∇u_x ⋅ ∇v_x + ∇u_y ⋅ ∇v_y) dΩ
                k[I-1,2J-1] += μ * (B₁[i]*B₁[j] + B₂[i]*B₂[j]) * 𝑤
                k[I,2J]     += μ * (B₁[i]*B₁[j] + B₂[i]*B₂[j]) * 𝑤
            end
        end
    end
end

#===== 压力-速度耦合项：∫p div(v) dΩ → 矩阵 B =====#
function (op::Operator{:∫pdivvdxdy})(aᵤ::T, aₚ::S; k::AbstractMatrix{Float64}) where {T<:AbstractElement, S<:AbstractElement}
    𝓒ᵤ = aᵤ.𝓒; 𝓒ₚ = aₚ.𝓒
    𝓖 = aᵤ.𝓖
    for ξ in 𝓖
        Nₚ = ξ[:𝝭]     # 压力形函数
        B₁ = ξ[:∂𝝭∂x]  # 速度形函数 x 导数
        B₂ = ξ[:∂𝝭∂y]  # 速度形函数 y 导数
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒ₚ)
            I = xᵢ.𝐼   # 压力自由度索引
            for (j,xⱼ) in enumerate(𝓒ᵤ)
                J = xⱼ.𝐼
                # B 矩阵：-∫ p (∂v_x/∂x + ∂v_y/∂y) dΩ
                k[I, 2J-1] -= Nₚ[i] * B₁[j] * 𝑤  # x 方向
                k[I, 2J]   -= Nₚ[i] * B₂[j] * 𝑤  # y 方向
                # B^T 矩阵：-∫ q (∂u_x/∂x + ∂u_y/∂y) dΩ
                k[2J-1, I] -= Nₚ[i] * B₁[j] * 𝑤  # 转置对称性
                k[2J, I]   -= Nₚ[i] * B₂[j] * 𝑤
            end
        end
    end
end

#===== 体力项：∫b⋅v dΩ → 向量 f =====#
function (op::Operator{:∫bvdxdy})(aᵤ::T; f::AbstractVector{Float64}) where T<:AbstractElement
    𝓒 = aᵤ.𝓒; 𝓖 = aᵤ.𝓖
    b = op.b  # 体力向量 [b_x, b_y]
    for ξ in 𝓖
        N = ξ[:𝝭]  # 速度形函数
        𝑤 = ξ.𝑤
        for (i,xᵢ) in enumerate(𝓒)
            I = xᵢ.𝐼
            # 体力贡献到 x 和 y 方向
            f[2I-1] += N[i] * b[1] * 𝑤
            f[2I]   += N[i] * b[2] * 𝑤
        end
    end
end