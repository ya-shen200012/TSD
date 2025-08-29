using Random
using LinearAlgebra
using Printf
using Dates
using Pkg
Pkg.add("CUTEst")
Pkg.add("DataFrames")
Pkg.add("CSV")
using CUTEst
using DataFrames, CSV
using NLPModels

# TSD method

function zoom(phi, dphi, alpha_lo, alpha_hi, phi0, dphi0, c1, c2; max_iter=35)
    for _ in 1:max_iter
        alpha_j = (alpha_lo + alpha_hi) / 2
        phi_j = phi(alpha_j)

        if (phi_j > phi0 + c1 * alpha_j * dphi0) || (phi_j >= phi(alpha_lo))
            alpha_hi = alpha_j
        else
            dphi_j = dphi(alpha_j)

            if abs(dphi_j) <= -c2 * dphi0
                return alpha_j,0
            end

            if dphi_j * (alpha_hi - alpha_lo) >= 0
                alpha_hi = alpha_lo
            end

            alpha_lo = alpha_j
        end
    end

    #@warn "zoom: max_iter reached without satisfying strong Wolfe conditions"
    return (alpha_lo + alpha_hi) / 2,0
end

function line_search_strong_wolfe(f, grad_f, x, p; c1=1e-4, c2=0.1, alpha_init=1.0, alpha_max=1e4, max_iter=20)
    phi(alpha) = f(x + alpha * p)
    dphi(alpha) = dot(grad_f(x + alpha * p), p)

    phi0 = phi(0.0)
    dphi0 = dphi(0.0)

    alpha_prev = 0.0
    alpha_curr = alpha_init
    phi_prev = phi0

    for i in 1:max_iter
        phi_curr = phi(alpha_curr)

        if (phi_curr > phi0 + c1 * alpha_curr * dphi0) || (i > 1 && phi_curr >= phi_prev)
            return zoom(phi, dphi, alpha_prev, alpha_curr, phi0, dphi0, c1, c2; max_iter=max_iter)
        end

        dphi_curr = dphi(alpha_curr)

        if abs(dphi_curr) <= -c2 * dphi0
            return alpha_curr,0
        end

        if dphi_curr >= 0
            return zoom(phi, dphi, alpha_curr, alpha_prev, phi0, dphi0, c1, c2; max_iter=max_iter)
        end

        alpha_prev = alpha_curr
        phi_prev = phi_curr
        alpha_curr = min(2 * alpha_curr, alpha_max)
    end

    @warn "line search: max_iter reached without satisfying strong Wolfe conditions"
    flag=1
    return alpha_curr,flag
end

function sd_2_algorithm(x0::Vector{Float64}, fobj::Function, grad::Function, SD_2N::Int, tol::Float64, MaxIter::Int)
    x = copy(x0)
    iter = 0
    g = grad(x)
    gnorm = maximum(abs.(g))
    p = copy(g)
    g00 = zeros(length(g))
    g01 = zeros(length(g))
    flag=0
    # println("begin while")
    while ((gnorm > tol) && (iter < MaxIter))
        # println("A new cycle")
        # println(gnorm)
        alpha00 = 0.0
        alpha01 = 0.0
        for i_sd_2 in 1:SD_2N
            # about TSD_algorithm 
            # recall function Strong_Wolfe_Condition(phi::Function,grad_swc::Function,gp::Float64,p::Vector{Float64})
            if i_sd_2 == 1
                # println(iter)
                # phi = alpha -> fobj(x - alpha * p) # phi(alpha)=fobj(x-alpha*p)
                # grad_swc = alpha -> grad(x - alpha * p) # grad_swc(alpha)=grad(x-alpha*p)
                # gp=dot(p, g) # positive number
                # alpha=Strong_Wolfe_Condition(phi,grad_swc,gp,p,20) # positive number, with - descent direction, format: x_new=x-alpha*p
 
                # function line_search_strong_wolfe(f, grad_f, x, p;
                                #   c1=1e-4, c2=0.1,
                                #   alpha_init=1.0,
                                #   alpha_max=100.0,
                                #   max_iter=20)
                alpha,flag0=line_search_strong_wolfe(fobj,grad,x,-p)
                flag=max(flag,flag0)
                # alpha = gp / dot(p, A .* p) # original stepsize, positive number, with - descent direction, format: x_new=x-alpha*p
                x_new = x - alpha * p
            else
                # phi = alpha -> fobj(x - alpha * g) # phi(alpha)=fobj(x-alpha*g)
                # grad_swc = alpha -> grad(x - alpha * g) # grad_swc(alpha)=grad(x-alpha*g)
                # gp=dot(g, g)
                # p=copy(g)
                # alpha=Strong_Wolfe_Condition(phi,grad_swc,gp,p,20) # positive number, format: x_new=x+alpha*p
                p=copy(g)                
                alpha,flag0=line_search_strong_wolfe(fobj,grad,x,-p)
                if isnan(alpha)||isinf(alpha)
                    error("alpha is NaN or ±Inf in the $(iter) th iteration！")
                end
                flag=max(flag,flag0)
                x_new = x - alpha * p
                bad_indices = findall(xi -> isnan(xi) || isinf(xi), x)
                if !isempty(bad_indices)
                    error("Detect NaN or ±Inf in x in $(iter) th iteration, location: $bad_indices")
                end
                bad_indices_2 = findall(xi -> isnan(xi) || isinf(xi), p)
                if !isempty(bad_indices_2)
                    # println(x)
                    # println(p)
                    error("Detect NaN or ±Inf in p in $(iter) th iteration, location: $bad_indices_2")
                end
                # alpha = dot(g, g) / dot(g, A * g) # original stepsize, positive number, with - descent direction, format: x_new=x-alpha*p
                # x_new = x - alpha * g
            end

            if i_sd_2 == SD_2N - 1
                alpha00 = alpha
                g00 = copy(g)
            end

            if i_sd_2 == SD_2N
                alpha01 = alpha
                g01 = copy(g)
            end
            # update data
            x = x_new
            iter += 1
            g = grad(x)
            gnorm = maximum(abs.(g))

            if gnorm < tol
                return x, iter, gnorm, flag    
            end
            # println(gnorm)
        end # for i_sd_2 in 1:SD_2N

        p = alpha00 * g00 + alpha01 * g01

        if iter == MaxIter
            @warn "iter == MaxIter"
        end
    end # while ((gnorm > tol) && (iter < MaxIter))

    return x, iter, gnorm,flag
end

function sd_2(x0::Vector{Float64}, A::AbstractMatrix, SD_2N::Int, grad::Function, tol::Float64, MaxIter::Int)
    x = copy(x0)
    iter = 0
    g = grad(x)
    gnorm = maximum(abs.(g))
    p = copy(g)
    n = size(A, 2)
    g00 = zeros(length(g))
    g01 = zeros(length(g))

    while ((gnorm > tol) && (iter < MaxIter))
        alpha00 = 0.0
        alpha01 = 0.0

        for i_sd_2 in 1:SD_2N
            if i_sd_2 == 1
                if n == 1
                    alpha = dot(p, g) / dot(p, A .* p)
                else
                    alpha = dot(p, g) / dot(p, A * p)
                end
                x_new = x - alpha * p
            else
                if n == 1
                    alpha = dot(g, g) / dot(g, A .* g)
                else
                    alpha = dot(g, g) / dot(g, A * g)
                end
                x_new = x - alpha * g
            end

            if i_sd_2 == SD_2N - 1
                alpha00 = alpha
                g00 = copy(g)
            end

            if i_sd_2 == SD_2N
                alpha01 = alpha
                g01 = copy(g)
            end
            # update data
            x = x_new
            iter += 1
            g = grad(x)
            gnorm = maximum(abs.(g))

            if gnorm < tol
                return x, iter, gnorm
            end
        end # for i_sd_2 in 1:SD_2N

        p = alpha00 * g00 + alpha01 * g01

        if iter == MaxIter
            @warn "iter == MaxIter"
        end
    end # while ((gnorm > tol) && (iter < MaxIter))

    return x, iter, gnorm
end

# BBQ Method for Unconstrained Optimization in Julia

# You need to provide two functions:
#   fobj(x)   - the objective function
#   grad(x)   - the gradient function

function BBQ_unc(x0, fobj, grad; opts=Dict(), LinesearchParms=Dict())
    # Set default options
    opts = merge(Dict(
        "PrintLevel" => 1,
        "MaxIter" => 200000,
        "MaxFunc" => 999999,
        "tol" => 1e-6,
        "bbmin" => 1e-10,
        "bbmax" => 1e6,
        "tau" => 0.2,
        "gama" => 1.02,
        "record" => 0
    ), opts)

    LinesearchParms = merge(Dict(
        "delta" => 1e-4,
        "sigma" => 0.5,
        "M" => 10,
        "maxlsf" => 100
    ), LinesearchParms)

    # Initialize
    t0=time()
    x = copy(x0)
    fk = fobj(x)
    g = grad(x)
    nfunc = 1
    ngrad = 1
    iter = 0
    gnorm = maximum(abs.(g))
    lambda = 1.0
    alp = 1.0 / gnorm
    xbest = copy(x)
    fbest = fk
    gnorm_best = gnorm
    LastFvalues = fill(-1e99, LinesearchParms["M"])
    LastFvalues[1] = fk
    mem_sty = zeros(2)
    record = opts["record"] == 1
    xtemp=similar(x)
    ftemp=0
    code=0
    bb1=0
    bb2=0
    # History
    hist = Dict("bbs" => Float64[], "step" => Float64[], "obj" => Float64[], "gnorm" => Float64[])
    
    if opts["PrintLevel"] >= 1
        println(" iter        obj        gnorm       step")
    end
    # test optimality condition for the initial point
    CheckStop = 0;
    if ( gnorm <= opts["tol"] )
        CheckStop = 1;
    end
    if ( CheckStop == 1 )
        if opts["PrintLevel"] == 2
            @printf("%5d  %12.4e  %10.4e  %8.4f\n", iter, fk, gnorm, lambda)
        end
        out = Dict(
        "code" => 0,
        "iter" => iter,
        "nfunc" => nfunc,
        "ngrad" => ngrad,
        "x" => xbest,
        "obj" => fbest,
        "gnorm" => gnorm_best,
        "cputime" => time()-t0# you may replace with @elapsed
        )
        return out, hist
    else
        xnorm = maximum(abs.(x))
        if ( xnorm > 0.0 )
            alp = xnorm/gnorm;
        else
            alp = 1.0/gnorm;
        end
    end
    # History
    hist = Dict("bbs" => Float64[], "step" => Float64[], "obj" => Float64[], "gnorm" => Float64[])

    while iter < opts["MaxIter"]
        if record
            push!(hist["bbs"], alp)
            push!(hist["step"], lambda)
            push!(hist["obj"], fk)
            push!(hist["gnorm"], gnorm)
        end

        if gnorm <= opts["tol"]
            fbest = fk;
            xbest = copy(x);
            gnorm_best = gnorm;
            code = 0;
            break
        end

        iter += 1
        if nfunc >= opts["MaxFunc"]
            code = 1
            break
        end

        # GLL line search
        px = x .- alp .* g
        d = px .- x
        gtd = dot(g, d)
        fr = maximum(LastFvalues)
        alpha = 1.0
        nfunco = nfunc

        while true
            xtemp = x .+ alpha .* d
            ftemp = fobj(xtemp)
            nfunc += 1
            if ftemp <= fr + LinesearchParms["delta"] * alpha * gtd
                break
            end

            if alpha <= 0.1
                alpha *= LinesearchParms["sigma"]
            else
                atemp = -gtd * alpha^2 / (2 * (ftemp - fk - alpha * gtd))
                if (atemp < 0.1 || atemp > (0.9 * alpha))
                    atemp = LinesearchParms["sigma"] * alpha
                end
                alpha = atemp
            end

            if (nfunc >= (nfunco + LinesearchParms["maxlsf"]))
                code=3
                error("Line search failed")
            end
        end

        lambda = alpha
        sx = xtemp .- x
        gold = g
        # update x and gradient
        x = xtemp
        fk = ftemp
        g = grad(x)
        ngrad += 1
        gnorm = maximum(abs.(g))

        if opts["PrintLevel"] == 2
            @printf("%5d  %12.4e  %10.4e  %8.4f\n", iter, fk, gnorm, lambda)
        end

        # record recent M objective values
        LastFvalues[mod(iter+1, LinesearchParms["M"]) + 1] = fk
        #LastFvalues[mod(iter, LinesearchParms["M"]) + 1] = fk

        # Update best
        if fk < fbest
            xbest = copy(x)
            fbest = fk
            gnorm_best = gnorm
        end

        #compute new stepsize
        y = g - gold
        sty = dot(sx, y)
        if ( sty > 0 )
            mem_sty[mod(iter,2)+1] = 1
            if ( iter>0 && mem_sty[mod(iter-1,2)+1]==1 )
                bb1_old = bb1
                bb2_old = bb2
            end
        
            bb1 = dot(sx,sx)/sty
            bb2 = (sty)/dot(y,y)     
            if ( iter > 1 )
                if ( bb2/bb1 < opts["tau"] && sum(mem_sty)==2 )
                    phi13 = (bb2_old - bb2)/(bb2_old * bb2 * (bb1_old - bb1))
                    phi23 = phi13 * bb1 + 1/bb2
                    bb_new = 2/( phi23 + sqrt(phi23^2 - 4*phi13) )
                    if ( bb_new <= 0 )
                        alp = min(bb2_old,bb2)
                    else
                        alp = min(bb_new,bb2,bb2_old)
                    end
                    opts["tau"] = opts["tau"]/opts["gama"]
                else
                    alp = bb1
                    opts["tau"] = opts["tau"]*opts["gama"]
                end
            else
                alp = bb1
            end
        else
            mem_sty = zeros(2,1)
            xnorm = maximum(abs.(x))
            alp = min(1.0/gnorm, xnorm/gnorm)
        end

        alp = clamp(alp, opts["bbmin"], opts["bbmax"])
    end

    if ( iter == opts["MaxIter"] )
        code= 2
        if ( opts["PrintLevel"] > 1 )
            println("Maximum number of function evaluations reached.\n")
        end
    end

    if opts["PrintLevel"] >= 1
        @printf("%5d  %12.4e  %10.4e  %8.4f\n", iter, fk, gnorm, lambda)
    end

    if record
        push!(hist["bbs"], alp)
        push!(hist["step"], lambda)
        push!(hist["obj"], fk)
        push!(hist["gnorm"], gnorm)
    end

    out = Dict(
        "code" => code,
        "iter" => iter,
        "nfunc" => nfunc,
        "ngrad" => ngrad,
        "x" => xbest,
        "obj" => fbest,
        "gnorm" => gnorm_best,
        "cputime" => time()-t0 # you may replace with @elapsed
    )

    return out, hist
end

# test begins 
t_initial_start = now()
problem_inx=2 # 1 for quadratic problems. 2 for non-quadratic problems.
const kappa = 1e2
const tola = [1e-6, 1e-9, 1e-12]
const Maxiter = Int(1e4)
if problem_inx==1
    const n = Int(1e4)
    const ITER = 10
    s_iter = zeros(16, ITER+1)
    s_gnorm = zeros(16, ITER+1)
    s_time2 = zeros(16, ITER+1)
    s1_iter=zeros(3)
    s1_time2=zeros(3)
    s2_iter=zeros(3)
    s2_time2=zeros(3)
else
    const ITER = 1
    s_i_g_t=zeros(16,3)
end
# initial values, global
# 提前声明为函数类型（或更宽泛）Declare it as a function type (or more general) in advance
grad_my = nothing
fobj = nothing
# 初始化空字典，键是算法名，值是二维数组 Initialize an empty dictionary, the key is the algorithm name, and the value is a two-dimensional array
results = Dict{String, DataFrame}()
failed = String[]                   # 存储失败问题名称 Storage failure problem name
tol=0

# find the problem 
filtered_problems = select_sif_problems(contype="unc",max_var=10000)
for filtered_problems_i in eachindex(filtered_problems[1:95])
    println(filtered_problems_i)
    postinx=0
    # println(length(filtered_problems))

    nlp = CUTEstModel{Float64}(filtered_problems[filtered_problems_i+postinx])
    n=size(nlp.meta.x0)[1]

    fx = obj(nlp, nlp.meta.x0)
    gx = grad(nlp, nlp.meta.x0)
    Hx = hess(nlp, nlp.meta.x0)
    cx = cons(nlp, nlp.meta.x0)
    Jx = jac(nlp, nlp.meta.x0)
    print(nlp.meta)

    # A = zeros(2, 3, 4)     # 全部是 0.0 的三维数组, 引用格式 A[i,j,k]. A three-dimensional array of all 0.0 values, referenced in the format A[i,j,k].
    # 设置随机种子 Setting the random seed
    Random.seed!(20230629)

    if problem_inx==1
        # 构造 d（对应不同的条件数设置）Construct d (corresponding to different condition number settings)
        d = zeros(n)
        d[1:div(n,2)] .= 1 .+ (kappa - 1) .* (0.8 .+ 0.2 .* rand(div(n,2)))
        d[(div(n,2)+1):end] .= 1 .+ (kappa - 1) .* (0.2 .* rand(div(n,2)))
        d[1] = 1
        d[end] = kappa

        # A = d
        A=reshape(d,n,1)
        xs = -10 .+ 20 .* rand(n)
    end

    # 梯度和目标函数 Gradients and Objective Functions
    global grad_my
    global fobj
    if problem_inx==1
        if size(A,2) == 1
            grad_my = x->A .* (x .- xs)
            fobj =x-> dot(x .- xs, A .* (x .- xs))
        else
            grad_my = x->A * (x .- xs)
            fobj = x->dot(x .- xs, A * (x .- xs))
        end
    else
        fobj = x->obj(nlp, x)
        grad_my = x->grad(nlp, x)
    end

    # tau_list = 0.1:0.1:0.9 # for BBQ method
    tau_list = [0.2] # for BBQ method
    TSD_list = [10,50]
    # TSD_list = [10]

    for j in 1:1
        try
            println("------ tol index $j: tol = ", tola[j], " ------")
            for i in 1:ITER
                println("Iteration $i")
                if problem_inx==1
                    x0 = -10 .+ 20 .* rand(n)
                else
                    x0=nlp.meta.x0
                end
                global tol
                if problem_inx==1
                    tol = tola[j] * maximum(abs.(grad_my(x0)))
                else
                    tol = tola[j]
                end
                
                for (tau_idx, tau) in enumerate(tau_list)
                    println("Running BBQ with τ = $tau")
                    opts = Dict(
                        "PrintLevel" => 0,
                        "tol" => tol,
                        "gama" => 1.0,
                        "tau" => tau
                    )

                    t_start = now()
                    out, hist = BBQ_unc(x0, fobj, grad_my, opts=opts)
                    t_end = now()
                    if problem_inx==1
                        s_iter[tau_idx, i] = out[:"iter"]
                        s_gnorm[tau_idx, i] = out[:"gnorm"]
                        s_time2[tau_idx, i] = Millisecond(t_end - t_start).value / 1000.0  # seconds
                    else
                        s_i_g_t[tau_idx,1]=out[:"iter"]
                        s_i_g_t[tau_idx,2]=out[:"gnorm"]
                        s_i_g_t[tau_idx,3]=Millisecond(t_end - t_start).value / 1000.0  # seconds
                    end
                end
                if problem_inx==1
                    for (TSD_idx, TSD) in enumerate(TSD_list)
                        println("Running TSD with parameter = $TSD")

                        t_start = now()
                        x1,iter,gnorm = sd_2(x0, A, TSD, grad_my, tol, Maxiter)
                        t_end = now()
                        s_iter[TSD_idx+9, i] = iter
                        s_gnorm[TSD_idx+9, i] = gnorm
                        s_time2[TSD_idx+9, i] = Millisecond(t_end - t_start).value / 1000.0  # seconds
                    end
                else 
                    for (TSD_idx, TSD) in enumerate(TSD_list)
                        println("Running TSD with parameter = $TSD")

                        t_start = now()
                        # function sd_2_algorithm(x0::Vector{Float64}, fobj::Function, grad::Function, SD_2N::Int, tol::Float64, MaxIter::Int)
                        # println("sd_2_algorithm begins")
                        x1,iter,gnorm,flag = sd_2_algorithm(x0, fobj, grad_my, TSD, tol, Maxiter)
                        t_end = now()
                        s_i_g_t[TSD_idx+9, 1] = iter
                        s_i_g_t[TSD_idx+9, 2] = gnorm
                        s_i_g_t[TSD_idx+9, 3] = Millisecond(t_end - t_start).value / 1000.0  # seconds
                        if flag==1
                            timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
                            msg = "[$timestamp] TSD: Strong Wolfe line search failed on 问题 $(string(filtered_problems_i+postinx,"  ",nlp.meta.name))\n"
                            println(msg)

                            # 写入 log 文件（以追加模式）Write to log file (in append mode)
                            open("error_log_TSD.txt", "a") do io
                                write(io, msg)
                            end
                        end
                    end
                end
            end
            # 存入字典 Store in dictionary
            colnames = [:Iterations, :GradientNorm, :Time]
            s_i_g_t_df = DataFrame(s_i_g_t, colnames)
            results[string(filtered_problems_i+postinx,"  ",nlp.meta.name)] = s_i_g_t_df

            # 初始化一个空 DataFrame 用于合并 Initialize an empty DataFrame for merging
            all_results = DataFrame()

            # 遍历每个算法及其结果 Iterate over each algorithm and its results
            for (algoname, df) in results
                # 给每个 DataFrame 添加一列 Algorithm，值为算法名 Add a column Algorithm to each DataFrame, with the value being the algorithm name
                df_alg = copy(df)
                df_alg.Algorithm = fill(algoname, nrow(df))  # 添加列 Add columns
                append!(all_results, df_alg)
            end

            # 把列顺序调整一下（可选）Adjust the column order (optional)
            select!(all_results, :Algorithm, Not(:Algorithm))

            # 输出到 CSV 文件 Output to CSV file
            CSV.write("combined_results_20250819.csv", all_results)
        catch e
            timestamp = Dates.format(now(), "yyyy-mm-dd HH:MM:SS")
            msg = "[$timestamp] skip problem $(string(filtered_problems_i+postinx,"  ",nlp.meta.name)), error type:$(typeof(e)), error information:$(e)\n"
            println(msg)

            # 写入 log 文件（以追加模式）Write to log file (in append mode)
            open("error_log_20250819.txt", "a") do io
                write(io, msg)
            end

            push!(failed, string(filtered_problems_i+postinx,"  ",nlp.meta.name))
        end
        if problem_inx==1
            for i in 1:length(tau_list)
                s_iter[i,ITER+1]=sum(s_iter[i,1:ITER])/ITER
                s_gnorm[i,ITER+1]=sum(s_gnorm[i,1:ITER])/ITER
                s_time2[i,ITER+1]=sum(s_time2[i,1:ITER])/ITER
            end
            (s1_iter[j], s2_iter[j])=findmin(s_iter[1:9,ITER+1])
            (s1_time2[j], s2_time2[j])=findmin(s_time2[1:9,ITER+1])
        end
    end
    if problem_inx==1
        println(s1_iter)
        println(s1_time2)
        for i in 1:16
            println(s_iter[i,:])
        end
        for i in 1:16
            println(s_time2[i,:])
        end
    end

    finalize(nlp)
end #for
t_final_end=now()
println("The code has ran for $(Millisecond(t_final_end - t_initial_start).value / 1000.0) seconds.")
