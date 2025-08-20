
function obtain_ensol(L0, H0, delta, T, alpha, Ls, params, tspan, e_ops, tlist)
        Ls_par, H_par, He_par = mo.obtain_parametric_unraveling_operators(L0, H0, T, alpha, NLEVELS)
        sys = cs.System(H_par(delta), Ls, NLEVELS, nchannels)
        return mc.get_sol_jumps(sys, params, tspan, e_ops, tlist;
                save_on=false, save_start=false, save_end=false)
end

function get_lindbladsol(delta, tspan, u0; saveat)
        function rf_ode!(dr, r, p, t)
                dr[1] = -delta * r[2] - 0.5 * GAMMA * r[1]
                dr[2] = -OMEGA * r[3] + delta * r[1] - 0.5 * GAMMA * r[2]
                dr[3] = OMEGA * r[2] - GAMMA * (r[3] + 1)
        end
        prob = ODEProblem(rf_ode!, u0, tspan; saveat=tlist)
        return solve(prob)
end

function check_convergencetolindblad(delta, tspan, u0, sol_lindblad, sim, tlist, tolerance)
        r_mean = mo.average_expvals(sim)
        belowtolerance = true
        ntimes = length(tlist)
        # r_lindblad = sol_lindblad.u
        for k in 1:ntimes
                # difference_x = 0.5 * (r_lindblad[1] - r_mean[1, k])
                # difference_y = 0.5 * (r_lindblad[2] - r_mean[2, k])
                # difference_z = 0.5 * (r_lindblad[3] - r_mean[3, k])
                difference = 0.5 * ((sol_lindblad.u[k][1] - r_mean[1, k]) * e_ops[1] +
                                    (sol_lindblad.u[k][2] - r_mean[2, k]) * e_ops[2] +
                                    (sol_lindblad.u[k][3] - r_mean[3, k]) * e_ops[3])

                trace_distance = abs(0.5 * tr(sqrt(adjoint(difference) * difference)))
                if trace_distance > tolerance
                        return !belowtolerance
                end
        end
        return belowtolerance
end

function extractfisample!(fisample, sim)
        ntraj = length(fisample)
        for k in 1:ntraj
                normalize!(sim[k].u[end])
                fisample[k] = (2 * real(dot(sim[k].u[end],
                        sim[k].prob.kwargs[:callback].continuous_callbacks[1].affect!.cache_phi)))^2
        end
end


function extracttimessample!(tsample, sim)
        ntraj = size(tsample)
        for k in 1:ntraj
                normalize!(sim[k].u[end])
                fisample[k] = (2 * real(dot(sim[k].u[end],
                        sim[k].prob.kwargs[:callback].continuous_callbacks[1].affect!.cache_phi)))^2
        end
end



function generate_checkedsamples!(samplefi, sampleexpvals,
        u0, L0, H0, delta, ddelta, T, alpha, Ls, dLs, params, tspan, e_ops, tlist, tolerance)
        sim = obtain_ensol(L0, H0, delta, ddelta, T, alpha, Ls, dLs, params, tspan, e_ops, tlist)
        # sol_lindblad = get_lindbladsol(delta, tspan, u0; saveat=tlist)
        # convergence_flag = check_convergencetolindblad(delta, tspan, u0, sol_lindblad, sim, tlist, tolerance)
        # if !convergence_flag
        #         @warn "Tolerance surpassed with delta=$delta, alpha=$alpha"
        # end
        extractexpvalsample!(sampleexpvals, sim)
        extractfisample!(samplefi, sim)
end
# (L0, H0, delta, ddelta, T, alpha, Ls, dLs, params, tspan, e_ops, tlist)

function generate_samples!(samplefi, sampleexpvals,
        L0, H0, delta, ddelta, T, alpha, Ls, dLs, params, tspan, e_ops, tlist)
        sim = obtain_ensol(L0, H0, delta, ddelta, T, alpha, Ls, dLs, params, tspan, e_ops, tlist)
        extractexpvalsample!(sampleexpvals, sim)
        extractfisample!(samplefi, sim)
end


function extractexpvalsample!(sampleexpvals, sim)
        nops, ntimes, ntraj = size(sampleexpvals)
        for k in 1:ntraj
                for j in 1:ntimes
                        for i in 1:nops
                                sampleexpvals[i, j, k] = sim[k].prob.kwargs[:callback].discrete_callbacks[1].affect!.func.expvals[i, j]
                        end
                end
        end
end

function obtain_ensol(L0, H0, delta, T, alpha, params, tspan, e_ops, tlist)
        return mc.get_sol_jumps(cs.System(H0(delta), L0(delta), T, alpha),
                params, tspan, e_ops, tlist;
                save_on=false, save_start=false, save_end=false)
end

macro timeerr(ex)
        quote
                local stats = Base.gc_num()
                local t0 = time_ns()
                local v = $(esc(ex))
                local t1 = time_ns()
                local diff = Base.GC_Diff(Base.gc_num(), stats)

                local elapsed = (t1 - t0) / 1e9
                local allocd = diff.allocd
                local gc_time = diff.total_time / 1e9

                println(stderr,
                        "  Time: ", round(elapsed, digits=6), " s, ",
                        "GC time: ", round(gc_time, digits=6), " s, ",
                        "alloc: ", Base.format_bytes(allocd))

                v
        end
end
