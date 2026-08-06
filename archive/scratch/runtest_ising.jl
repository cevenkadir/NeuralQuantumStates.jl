# %% import required packages
using NeuralQuantumStates: Lattices, Hilberts, Operators
using NeuralQuantumStates.Extras
using Random
using MetaGraphsNext
using Graphs
using Test
using StaticArrays
using LinearAlgebra

using BenchmarkTools

# %% Cell 1
lat = Lattices.build(:Hypercube, [8], 1.0; periodic=[true])
hs = Hilberts.build(:Spin, 1 // 2, Lattices.nv(lat); array_type=Array)
ha = Operators.build(:TransverseFieldIsing, hs, lat; J=1.0, h_x=1.0, h_z=1.0)

a = Operators.connected_basis_configs(
    ha,
    [1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, -1 // 2, 1 // 2, -1 // 2];
    #force_unique=true
)

s_prime, mels = Operators.connected_basis_configs(
    ha,
    [1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, -1 // 2, 1 // 2, -1 // 2]
)

a = connected_basis_configs(
    ha,
    [-1 // 2, 1 // 2, 1 // 2, 1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2];
    #force_unique=true
)

a = connected_basis_configs(
    ha,
    [-1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2];
    #force_unique=true
)

aa = connected_basis_configs(
    ha,
    [
        [1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, -1 // 2, 1 // 2, -1 // 2],
        [-1 // 2, 1 // 2, 1 // 2, 1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2],
        [-1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2],
        [1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2, -1 // 2, 1 // 2, 1 // 2],
        [1 // 2, -1 // 2, 1 // 2, 1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2]
    ];
)

@benchmark aa = connected_basis_configs(
    ha,
    $([
        [1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, -1 // 2, 1 // 2, -1 // 2],
        [-1 // 2, 1 // 2, 1 // 2, 1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2],
        [-1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2, -1 // 2, 1 // 2],
        [1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2, -1 // 2, 1 // 2, 1 // 2],
        [1 // 2, -1 // 2, 1 // 2, 1 // 2, 1 // 2, -1 // 2, -1 // 2, -1 // 2]
    ]);
)

aa = Operators.connected_basis_configs(
    ha,
    [
        1//2 1//2 -1//2 1//2 -1//2 -1//2 1//2 -1//2;
        -1//2 1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2;
        -1//2 1//2 -1//2 1//2 -1//2 1//2 -1//2 1//2;
        1//2 1//2 -1//2 -1//2 -1//2 -1//2 1//2 1//2;
        1//2 -1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2
    ];
)
s_prime, mels = Operators.connected_basis_configs(
    ha,
    [
        1//2 1//2 -1//2 1//2 -1//2 -1//2 1//2 -1//2;
        -1//2 1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2
    ];
)

@benchmark aa = connected_basis_configs(
    ha,
    $([
        1//2 1//2 -1//2 1//2 -1//2 -1//2 1//2 -1//2;
        -1//2 1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2;
        -1//2 1//2 -1//2 1//2 -1//2 1//2 -1//2 1//2;
        1//2 1//2 -1//2 -1//2 -1//2 -1//2 1//2 1//2;
        1//2 -1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2
    ]);
)

aa = connected_basis_configs(
    ha,
    reshape([
            1//2 1//2 -1//2 1//2 -1//2 -1//2 1//2 -1//2;
            -1//2 1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2;
            -1//2 1//2 -1//2 1//2 -1//2 1//2 -1//2 1//2;
            1//2 1//2 -1//2 -1//2 -1//2 -1//2 1//2 1//2;
            1//2 -1//2 1//2 1//2 1//2 -1//2 -1//2 -1//2
        ], 5, :, 8)
)
