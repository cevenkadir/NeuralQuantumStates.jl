abstract type AbstractOperator end
abstract type DiscreteOperator{T} <: AbstractOperator end

#abstract type IsingOperator{T} <: DiscreteOperator{T} end

abstract type LocalOperator{T} <: DiscreteOperator{T} end

function build(symbol::Symbol, args...; kwargs...)
    return build(Val(symbol), args...; kwargs...)
end
