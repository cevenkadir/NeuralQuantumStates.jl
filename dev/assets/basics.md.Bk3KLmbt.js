import{_ as a,c as n,o as i,ai as p}from"./chunks/framework.B5eQE0wt.js";const g=JSON.parse('{"title":"","description":"","frontmatter":{},"headers":[],"relativePath":"basics.md","filePath":"basics.md","lastUpdated":null}'),l={name:"basics.md"};function t(e,s,h,k,c,r){return i(),n("div",null,s[0]||(s[0]=[p(`<h2 id="installation" tabindex="-1">Installation <a class="header-anchor" href="#installation" aria-label="Permalink to &quot;Installation&quot;">​</a></h2><p>If you still want to try it out, you can install it from the Julia REPL by running:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">julia</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">&gt;</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> import</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Pkg; Pkg</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">add</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#032F62;--shiki-dark:#9ECBFF;">&quot;https://github.com/cevenkadir/NeuralQuantumStates.jl&quot;</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><h2 id="basics" tabindex="-1">Basics <a class="header-anchor" href="#basics" aria-label="Permalink to &quot;Basics&quot;">​</a></h2><p>You can start to use the package by runnung:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">using</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> NeuralQuantumStates</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">:</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Lattices, Hilberts, Operators</span></span></code></pre></div><p>As an example, let&#39;s construct the Hamiltonian of a transverse-field Ising chain. For that, first initialize a chain lattice of length 8:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">lat </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Lattices</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">build</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Hypercube</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, [</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">8</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">], </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">; periodic</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">[</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">true</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">])</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>NeuralQuantumStates.Lattices.Lattice{Int64, Float64, 1, 1}</span></span>
<span class="line"><span>  metagraph: MetaGraphsNext.MetaGraph{Int64, Graphs.SimpleGraphs.SimpleGraph{Int64}, Tuple{Int64}, StaticArraysCore.SVector{1, Float64}, Int64, Nothing, MetaGraphsNext.var&quot;#11#13&quot;, Float64}</span></span>
<span class="line"><span>  shape: StaticArraysCore.SVector{1, Int64}</span></span>
<span class="line"><span>  basis: NeuralQuantumStates.Lattices.LatticeBasis{Float64, 1, 1}</span></span>
<span class="line"><span>  periodic: StaticArraysCore.SVector{1, Bool}</span></span></code></pre></div><p>Then, choose a spin-1/2 as the basis for the Hilbert space:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">hil </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Hilberts</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">build</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:Spin</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, </span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> //</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, Lattices</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">nv</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(lat))</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>NeuralQuantumStates.Hilberts.FiniteUniformHilbert{Rational{Int64}, Array, 8, 2}</span></span>
<span class="line"><span>  lDoF: StaticArraysCore.SVector{2, Rational{Int64}}</span></span>
<span class="line"><span>  constraint: NeuralQuantumStates.Hilberts.NoDiscreteHilbertConstraint NeuralQuantumStates.Hilberts.NoDiscreteHilbertConstraint()</span></span>
<span class="line"><span>  type: Symbol Spin</span></span></code></pre></div><p>With these two, construct the Hamiltonian operator of the model:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">ham </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Operators</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">build</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">:TransverseFieldIsing</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, hil, lat; J</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, h_x</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">, h_z</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1.0</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">)</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>NeuralQuantumStates.Operators.TransverseFieldIsingOperator{Float64, Rational{Int64}, Array, Int64, Float64, Float64, 8, 2, 1, 1}</span></span>
<span class="line"><span>  hilbert: NeuralQuantumStates.Hilberts.FiniteUniformHilbert{Rational{Int64}, Array, 8, 2}</span></span>
<span class="line"><span>  lattice: NeuralQuantumStates.Lattices.Lattice{Int64, Float64, 1, 1}</span></span>
<span class="line"><span>  J: Float64 1.0</span></span>
<span class="line"><span>  h_x: Float64 1.0</span></span>
<span class="line"><span>  h_z: Float64 1.0</span></span></code></pre></div><p>To test it, let&#39;s find out the connected basis configurations, <code>s_prime</code>, to the Hamiltonian for the following 2 basis configurations with their matrix elements, <code>mels</code>:</p><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">s_prime, mels </span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">=</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;"> Operators</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">.</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">connected_basis_configs</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">(</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ham,</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    [</span></span>
<span class="line"><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">        1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">        -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;"> 1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;"> -</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">1</span><span style="--shiki-light:#D73A49;--shiki-dark:#F97583;">//</span><span style="--shiki-light:#005CC5;--shiki-dark:#79B8FF;">2</span><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">;</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">    ]</span></span>
<span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">);</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">s_prime</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>9×2×8 Array{Union{Missing, Rational{Int64}}, 3}:</span></span>
<span class="line"><span>[:, :, 1] =</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 2] =</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 3] =</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 4] =</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 5] =</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 6] =</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 7] =</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span> -1//2   1//2</span></span>
<span class="line"><span>  1//2  -1//2</span></span>
<span class="line"><span></span></span>
<span class="line"><span>[:, :, 8] =</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span> -1//2  -1//2</span></span>
<span class="line"><span>  1//2   1//2</span></span></code></pre></div><div class="language-julia vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang">julia</span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span style="--shiki-light:#24292E;--shiki-dark:#E1E4E8;">mels</span></span></code></pre></div><div class="language- vp-adaptive-theme"><button title="Copy Code" class="copy"></button><span class="lang"></span><pre class="shiki shiki-themes github-light github-dark vp-code" tabindex="0"><code><span class="line"><span>9×2 Matrix{Union{Missing, Float64}}:</span></span>
<span class="line"><span> -4.0  4.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span>
<span class="line"><span>  1.0  1.0</span></span></code></pre></div>`,21)]))}const o=a(l,[["render",t]]);export{g as __pageData,o as default};
