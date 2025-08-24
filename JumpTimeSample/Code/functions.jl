function get_lindbladsol(delta, tspan, u0; saveat)
        function rf_ode!(dr, r, p, t)
                dr[1] = -delta * r[2] - 0.5 * GAMMA * r[1]
                dr[2] = -OMEGA * r[3] + delta * r[1] - 0.5 * GAMMA * r[2]
                dr[3] = OMEGA * r[2] - GAMMA * (r[3] + 1)
        end
        prob = ODEProblem(rf_ode!, u0, tspan; saveat=tlist)
        return solve(prob)
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

function get_a(theta)
        Htilde = H0(theta) - 1.0im * conj(cfield) * L0(theta) - 0.5im * adjoint(L0(theta)) * L0(theta)
        a = [tr(pauli_basis[k] * Htilde) / tr(pauli_basis[k]^2) for k in 1:4]
        return a[1], a[2], a[3], a[4]
end


function get_Heffexponential(theta, tau)
        a0, a1, a2, a3 = get_a(theta)
        xi = sqrt(a1^2 + a2^2 + a3^2)
        dotproduct = a1 * util.sigma_x + a2 * util.sigma_y + a3 * util.sigma_z
        return exp(-1.0im * a0 * tau) * (cos(tau * xi) * I(2) - 1.0im * dotproduct / xi * sin(tau * xi)) #the alpha factor may be ommited by normalization
end

function get_Heffexponential!(cache, theta, tau)
        a0, a1, a2, a3 = get_a(theta)
        xi = sqrt(a1^2 + a2^2 + a3^2)
        cache .= exp(-1.0im * a0 * tau) * (cos(tau * xi) * I(2) - 1.0im * (a1 * util.sigma_x +
                                                                           a2 * util.sigma_y + a3 * util.sigma_z) / xi * sin(tau * xi)) #the alpha factor may be ommited by normalization
end

function monitoringstep!(cache_exp, cache_dexp, L, dL, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
        mul!(psi, cache_exp, cache_state)  # This is exp(-i\tau H_e)\psi_n
        ####### PHI UPDATE without rescaling
        # Obtain  \partial_\theta exp(-i\tau*H_eff(\theta))*\psi , store in aux1
        mul!(cache_aux1, cache_dexp, cache_state)
        # Obtain \partial_\theta exp(-i\tau*H_eff(\theta))*\psi + exp(-i\tau*H_eff(theta))*\phi, store where the derivative was
        mul!(cache_aux1, cache_exp, cache_phi, 1.0, 1.0)
        # Multiply by the jump operator and store in phi_cache
        mul!(cache_phi, L, cache_aux1)
        # Prepare the last term
        mul!(cache_aux2, dL, psi)
        # Now put everything together and store in cache_phi
        cache_phi .+= cache_aux2
        ###### STATE UPDATE without normalization
        mul!(cache_state, L, psi)
        ##### NORMALIZATION
        # Normalize phi
        normalization = norm(cache_state)
        lmul!(1 / normalization, cache_phi)
        # Normalize the after jump state
        lmul!(1 / normalization, cache_state)
end

function finalmonitoringstep!(cache_exp, cache_dexp, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
        mul!(psi, cache_exp, cache_state)
        mul!(cache_aux1, cache_dexp, cache_state)
        mul!(cache_aux2, cache_exp, cache_phi)
        cache_phi .= cache_aux1 + cache_aux2
        lmul!(1 / norm(psi), cache_phi)
        normalize!(psi)
end


function psiphi_finaltime!(sol, L, dL, cache_exp, cache_dexp, cache_aux1, cache_aux2, psi, cache_state, cache_phi, delta)
        jumptimes = sol.prob.kwargs[:callback].continuous_callbacks[1].affect!.jump_times
        njumps = sol.prob.kwargs[:callback].continuous_callbacks[1].affect!.jump_counter[] - 1
        if njumps == 0
                tau = tf
                get_Heffexponential!(cache_exp, delta, tau)
                ForwardDiff.derivative!(cache_dexp, theta -> get_Heffexponential(theta, tau), delta)
                finalmonitoringstep!(cache_exp, cache_dexp, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
        else
                tau = jumptimes[1]
                get_Heffexponential!(cache_exp, delta, tau)
                ForwardDiff.derivative!(cache_dexp, theta -> get_Heffexponential(theta, tau), delta)
                monitoringstep!(cache_exp, cache_dexp, L, dL, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
                if njumps > 1
                        for k in 2:njumps
                                tau = jumptimes[k] - jumptimes[k-1]
                                get_Heffexponential!(cache_exp, delta, tau)
                                ForwardDiff.derivative!(cache_dexp, theta -> get_Heffexponential(theta, tau), delta)
                                monitoringstep!(cache_exp, cache_dexp, L, dL, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
                        end
                end

                tau = tf - jumptimes[end]
                get_Heffexponential!(cache_exp, delta, tau)
                ForwardDiff.derivative!(cache_dexp, theta -> get_Heffexponential(theta, tau), delta)
                finalmonitoringstep!(cache_exp, cache_dexp, cache_state, psi, cache_aux1, cache_aux2, cache_phi)
        end
end
