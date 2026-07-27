# ── regressions/ — out-of-library textbook end-to-end cases ──────────────
#
# See README.md for the rules. In short: each case is ONE flat main() that
# consumes the PUBLIC library only, checks a cited LITERATURE value, runs in
# <= 60 s, and exits non-zero on failure.
#
# These are deliberately NOT registered with add_test(): CLAUDE.md's rule is
# that tests/ is unit-only and "if it is registered, it runs". Regressions are
# end-to-end, so they get their own target instead:
#
#   cmake -B build -DATV_BUILD_REGRESSIONS=ON
#   cmake --build build --target regressions      # builds AND runs every case

add_custom_target(regressions
    COMMENT "Running all end-to-end textbook regression cases")

# atv_add_regression(<domain> <case>)
#   expects regressions/<domain>/<case>/reg_<case>.cpp
function(atv_add_regression domain case)
    set(target reg_${case})
    add_executable(${target} ${domain}/${case}/${target}.cpp)
    target_link_libraries(${target} PRIVATE atomiverse)
    set_target_properties(${target} PROPERTIES CXX_STANDARD 20)

    add_custom_target(run_${target}
        COMMAND ${target}
        DEPENDS ${target}
        COMMENT "regression [${domain}/${case}]")
    add_dependencies(regressions run_${target})
endfunction()

# Layout guarantee (always on): out-of-library consumer of the core contract
# headers — backend-mirror-01-core moved these from facade/ into core/.
atv_add_regression(layout core_contract_includes)
atv_add_regression(layout backend_owned_longrange)

if(ATV_BUILD_MD)
    # Layout: full-word backend include paths (backend-mirror-03-paths).
    atv_add_regression(layout path_words)
    atv_add_regression(layout driver_names)
    atv_add_regression(layout no_backend)
    atv_add_regression(layout architecture_gates)
    atv_add_regression(md lj_argon)
    atv_add_regression(md const-pot-bo-parallel-plate)
    atv_add_regression(md per-electrode-eta)
    atv_add_regression(md const-charge-electrode)
    atv_add_regression(md const-pot-cg-large-electrode)
    atv_add_regression(md electrode-charge-predictor)
    atv_add_regression(md pr_npt_lj)
    atv_add_regression(md nph_lj)
    atv_add_regression(md surface-tension-coupling)
    atv_add_regression(md pe_bonded_nve)
    atv_add_regression(md eam_cu)
    atv_add_regression(md sw_si)
    atv_add_regression(md rigid_water)
    atv_add_regression(md dihedral_charmm)
    atv_add_regression(md improper_harmonic)
    atv_add_regression(md gaff_protocol_smoke)
    atv_add_regression(md charging-work-energy)
    atv_add_regression(md gaussian-electrode-charge)
    atv_add_regression(md freeze-electrode-atoms)
    atv_add_regression(md piston-electrode-pz)
    atv_add_regression(md const-pot-charge-observer)
    atv_add_regression(md capacitance-charge-fluctuation)
    atv_add_regression(md potential-profile-observer)
    atv_add_regression(md external-efield-md)
    atv_add_regression(md constant-force-pull-ion)
    atv_add_regression(md oscillating-efield-md)
    atv_add_regression(md four-site-tip4p-foothold)
    atv_add_regression(md virtual-site-multi-topology)
    atv_add_regression(md finite-field-single-electrode)
    atv_add_regression(md pim-induced-dipole-bulk)
    atv_add_regression(md thomas-fermi-electrode)
    atv_add_regression(md dipole-electrode-joint)
    # Ghost-charge transport after the const-pot Q solve under MPI DD:
    # self-relaunching binary (mpiexec -n 2), oracle=literature invariant.
    if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
        atv_add_regression(md const-pot-charge-halo-mpi)
        target_compile_definitions(reg_const-pot-charge-halo-mpi PRIVATE
            ATV_MPIEXEC="${MPIEXEC_EXECUTABLE}")
    endif()
endif()
if(ATV_BUILD_MC)
    atv_add_regression(mc gcmc_lj)
    atv_add_regression(mc npt_lj)
    atv_add_regression(mc widom_hs)
endif()
if(ATV_BUILD_PW)
    atv_add_regression(pw h_force_fd)
    atv_add_regression(pw h_force_analytic)
    atv_add_regression(pw h_stress_cubic)
    atv_add_regression(pw pw-lda-si)
    atv_add_regression(pw pw-analytic-stress-si)
    # In-tree ONCV Si UPF for NC-PW regressions (compile-time path pin).
    set(_atv_si_upf "${CMAKE_SOURCE_DIR}/tests/data/pw/Si_lda_nc.upf")
    if(TARGET reg_pw-lda-si)
        target_compile_definitions(reg_pw-lda-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    if(TARGET reg_pw-analytic-stress-si)
        target_compile_definitions(reg_pw-analytic-stress-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    if(ATV_BUILD_MINIMIZER)
        atv_add_regression(pw pw-variable-cell-si)
        if(TARGET reg_pw-variable-cell-si)
            target_compile_definitions(reg_pw-variable-cell-si PRIVATE
                ATV_PW_UPF_PATH="${_atv_si_upf}")
        endif()
    endif()
    # DFPT phonons at Γ (response): PW LSDA Si + Sternheimer response.
    atv_add_regression(response dfpt-phonon-si)
    if(TARGET reg_dfpt-phonon-si)
        target_compile_definitions(reg_dfpt-phonon-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Finite-q DFPT Sternheimer at X (response): same Si twin, q = (1/2,0,0).
    atv_add_regression(response dfpt-finite-q-si)
    if(TARGET reg_dfpt-finite-q-si)
        target_compile_definitions(reg_dfpt-finite-q-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Incommensurate-q DFPT: direct Sternheimer on primitive cell, q=(0.3,0,0).
    atv_add_regression(response dfpt-incommensurate-si)
    if(TARGET reg_dfpt-incommensurate-si)
        target_compile_definitions(reg_dfpt-incommensurate-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Γ DFPT with pure GGA (PBE) XC kernel (response).
    atv_add_regression(response dfpt-gga-si)
    if(TARGET reg_dfpt-gga-si)
        target_compile_definitions(reg_dfpt-gga-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Γ DFPT with collinear spin (LSDA, multiplicity=1 Si twin).
    atv_add_regression(response dfpt-spin-si)
    if(TARGET reg_dfpt-spin-si)
        target_compile_definitions(reg_dfpt-spin-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Γ DFPT with metallic occupation smearing (FD, restricted Si twin).
    atv_add_regression(response dfpt-metallic-si)
    if(TARGET reg_dfpt-metallic-si)
        target_compile_definitions(reg_dfpt-metallic-si PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Finite-q DFPT at X with metallic occupation smearing (FD).
    atv_add_regression(response dfpt-finite-q-smearing)
    if(TARGET reg_dfpt-finite-q-smearing)
        target_compile_definitions(reg_dfpt-finite-q-smearing PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Γ DFPT collinear spin + metallic smearing (FD, Si twin).
    atv_add_regression(response dfpt-spin-smearing)
    if(TARGET reg_dfpt-spin-smearing)
        target_compile_definitions(reg_dfpt-spin-smearing PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Finite-q DFPT collinear spin + metallic smear at X (FD).
    atv_add_regression(response dfpt-finite-q-spin-smearing)
    if(TARGET reg_dfpt-finite-q-spin-smearing)
        target_compile_definitions(reg_dfpt-finite-q-spin-smearing PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    # Incommensurate-q DFPT with collinear spin (direct Sternheimer).
    atv_add_regression(response dfpt-spin-incommensurate)
    if(TARGET reg_dfpt-spin-incommensurate)
        target_compile_definitions(reg_dfpt-spin-incommensurate PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}")
    endif()
    atv_add_regression(pw pw-pbe-gga-si)
    atv_add_regression(pw pw-spin-polarized-o2)
    atv_add_regression(pw pw-spin-gga)
    atv_add_regression(pw bands-dos-si)
    atv_add_regression(pw bands-spin-si)
    # Metallic occupation DOS (FD smear on diamond Si; integral to E_F).
    atv_add_regression(pw dos-metallic-si)
    # Metallic occupation bands (FD smear + high-symmetry path).
    atv_add_regression(pw bands-metallic-si)
    # Collinear spin + metallic smear bands (FD + nspin=2 diamond Si).
    atv_add_regression(pw bands-spin-metallic-si)
    # Collinear-spin DOS (fixed occ M=0 diamond Si; integral to joint VBM).
    atv_add_regression(pw dos-spin-si)
    # Collinear spin + metallic FD smear DOS (diamond Si; integral to E_F).
    atv_add_regression(pw dos-spin-metallic-si)
    atv_add_regression(pw berry-polarization)
    target_compile_definitions(reg_berry-polarization PRIVATE
        ATV_PW_UPF_PATH="${_atv_si_upf}")
    atv_add_regression(pw berry-spin-si)
    target_compile_definitions(reg_berry-spin-si PRIVATE
        ATV_PW_UPF_PATH="${_atv_si_upf}")
    atv_add_regression(pw periodic-wannier-si)
    target_compile_definitions(reg_periodic-wannier-si PRIVATE
        ATV_PW_UPF_PATH="${_atv_si_upf}")
    atv_add_regression(pw periodic-wannier-uks)
    target_compile_definitions(reg_periodic-wannier-uks PRIVATE
        ATV_PW_UPF_PATH="${_atv_si_upf}")
    atv_add_regression(pw pw-uks-kmesh-fe)
    atv_add_regression(pw pw-uks-fe-smearing)
    if(TARGET reg_pw-uks-fe-smearing)
        target_compile_definitions(reg_pw-uks-fe-smearing PRIVATE
            ATV_SOURCE_DIR="${CMAKE_SOURCE_DIR}")
    endif()
    # MPI k-point round-robin scaling gate: the binary re-launches itself
    # through the configured mpiexec at np=1 and np=4 (parity + wall ratio).
    if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
        atv_add_regression(pw mpi-pw-scaling-4rank)
        target_compile_definitions(reg_mpi-pw-scaling-4rank PRIVATE
            ATV_PW_UPF_PATH="${_atv_si_upf}"
            ATV_MPIEXEC="${MPIEXEC_EXECUTABLE}")
    endif()
    # pw-smearing-al merged when branch lands
    if(ATV_BUILD_MD)
        atv_add_regression(pw pw-bomd-h)
    endif()
    if(ATV_BUILD_MINIMIZER)
        atv_add_regression(pw h_cell_isotropic)
    atv_add_regression(pw pw-smearing-al)
    endif()
endif()
# FD supercell phonons (response): SW Si twin needs MD potentials; API is vib::Phonon.
if(ATV_BUILD_VIB AND ATV_BUILD_MD)
    atv_add_regression(pw fd-phonon-supercell)
    # Fourier ω(q) path from supercell IFCs (phonopy-style); same SW twin.
    atv_add_regression(response phonon-dispersion-si)
endif()
if(ATV_BUILD_GPW)
    atv_add_regression(gpw gpw-pbe-d3bj-ice)
    atv_add_regression(gpw gpw-d3-atm)
    atv_add_regression(gpw gpw-d4-water)
    atv_add_regression(gpw gpw-d4-atm)
    atv_add_regression(gpw gpw-d4-stress)
    atv_add_regression(gpw gpw-d4-forces)
    atv_add_regression(gpw gpw-ts-ice)
    atv_add_regression(gpw gpw-ts-forces)
    atv_add_regression(gpw gpw-xdm-ice)
    atv_add_regression(gpw gpw-xdm-forces)
    atv_add_regression(gpw gpw-mbd-ice)
    atv_add_regression(gpw gpw-mbd-periodic)
    atv_add_regression(gpw gpw-mbd-forces)
    atv_add_regression(gpw gpw-mbd-kpoint-forces)
    atv_add_regression(gpw gpw-vdw-df-ice)
    atv_add_regression(gpw gpw-vdw-df2-ice)
    atv_add_regression(gpw gpw-rvv10-ice)
    atv_add_regression(gpw gpw-vdw-scf-potential)
    atv_add_regression(gpw gpw-optb88-vdw)
    atv_add_regression(gpw gpw-mbd-hirshfeld)
    atv_add_regression(gpw gpw-hirshfeld-kpoints)
    atv_add_regression(gpw gpw-hirshfeld-uks)
    atv_add_regression(gpw ot-minimizer-water)
    atv_add_regression(gpw ot-uks-o2)
    atv_add_regression(gpw ot-kpoints)
    atv_add_regression(gpw ot-smearing)
    atv_add_regression(gpw ot-uks-smearing)
    atv_add_regression(gpw ot-uks-kpoints)
    atv_add_regression(gpw ot-dftu)
    atv_add_regression(gpw ot-dftu-kpoints)
    atv_add_regression(gpw ot-uks-dftu)
    atv_add_regression(gpw ot-dftu-smearing)
    atv_add_regression(gpw ot-uks-dftu-smearing)
    atv_add_regression(gpw gpw-water-pbe-mt)
    atv_add_regression(gpw gpw-blyp-water)
    atv_add_regression(gpw gpw-rpbe-water)
    atv_add_regression(gpw gpw-forces-water)
    atv_add_regression(gpw gpw-stress-water)
    atv_add_regression(gpw gpw-analytic-stress)
    atv_add_regression(gpw gpw-uks-o2-gamma)
    atv_add_regression(gpw gpw-uks-pbe)
    atv_add_regression(gpw gpw-uks-blyp)
    atv_add_regression(gpw gpw-uks-revpbe)
    atv_add_regression(gpw gpw-uks-rpbe)
    atv_add_regression(gpw gpw-uks-forces)
    atv_add_regression(gpw gpw-uks-stress)
    atv_add_regression(gpw gpw-uks-multigrid)
    if(ATV_BUILD_MD)
        atv_add_regression(gpw gpw-bomd-water)
        # Extended-Lagrangian BOMD (dynamics domain): XL auxiliary-density
        # propagation on the same H2O NVE twin.
        atv_add_regression(dynamics xl-bomd-water)
        # XL-BOMD on Γ-UKS (dual α/β auxiliary densities) — O2 triplet NVE.
        atv_add_regression(dynamics xl-bomd-uks)
        # Multi-k RKS XL-BOMD (complex P_k aux densities) — H2 NVE.
        atv_add_regression(dynamics xl-bomd-kpoints)
        atv_add_regression(dynamics cpmd-kpoints)
        # Multi-k UKS XL-BOMD (dual-spin P_k^sigma) — high-spin H2 NVE.
        atv_add_regression(dynamics xl-bomd-uks-kpoints)
        # Zero-SCF linearized-functional XL-BOMD forces (shadow BO potential).
        atv_add_regression(dynamics xl-bomd-linearized-forces)
        # Car–Parrinello MD (dynamics domain): fictitious electronic DOF on
        # the same H2O NVE twin family — not BOMD, not XL-BOMD.
        atv_add_regression(dynamics cpmd-water)
        # CPMD electronic RATTLE (constraint-preserving orbitals) — upgrade
        # of the Löwdin foothold on the same H2O NVE twin family.
        atv_add_regression(dynamics cpmd-rattle-water)
        # CPMD + Nosé–Hoover electronic thermostat on fictitious Ċ (H2O).
        atv_add_regression(dynamics cpmd-electronic-thermostat)
        # CPMD on Γ-UKS (dual α/β fictitious electronic DOF) — O2 triplet NVE.
        atv_add_regression(dynamics cpmd-uks)
        # Second-generation CP (Kühne ASPC + partial SCF corrector).
        atv_add_regression(dynamics second-gen-cpmd-water)
        # Second-generation CP on Γ-UKS (dual ASPC α/β + partial SCF).
        atv_add_regression(dynamics second-gen-cpmd-uks)
        # Second-generation CP with OT partial-SCF corrector (RKS Γ).
        atv_add_regression(dynamics second-gen-cpmd-ot)
        # Second-generation CP with OT partial-SCF corrector on Γ-UKS.
        atv_add_regression(dynamics second-gen-cpmd-uks-ot)
        # Multi-k second-generation CP (ASPC on complex P_k + partial SCF).
        atv_add_regression(dynamics second-gen-cpmd-kpoints)
    endif()
    atv_add_regression(gpw gpw-pbesol-si)
    atv_add_regression(gpw gpw-revpbe-si)
    atv_add_regression(gpw gpw-dftu-nio)
    atv_add_regression(gpw gpw-dftu-stress)
    atv_add_regression(gpw gpw-dftu-kpoints)
    atv_add_regression(gpw gpw-dftu-kpoint-stress)
    atv_add_regression(gpw gpw-kpoints-si)
    atv_add_regression(gpw periodic-wannier-gpw)
    # Dispersion ionic forces vs offline CP2K Quickstep D3 twins.
    atv_add_regression(gpw gpw-d3-forces)
    atv_add_regression(gpw gpw-d3bj-forces)
    # Gamma analytic cell stress folding periodic MBD@rsSCS.
    atv_add_regression(gpw gpw-mbd-stress)
    atv_add_regression(gpw gpw-metallic-smearing)
    atv_add_regression(gpw gpw-uks-smearing)
    atv_add_regression(gpw gpw-kpoint-forces)
    atv_add_regression(gpw gpw-kpoint-uks-forces)
    atv_add_regression(gpw gpw-kpoint-stress)
    atv_add_regression(gpw gpw-kpoint-uks-stress)
    atv_add_regression(gpw gpw-kpoint-uks-stress-gga)
    atv_add_regression(gpw gpw-kpoints-uks)
    # FD vibrational spectrum from GPW analytic forces (Si2 UKS triplet,
    # isolated, centred MT). Dual-gate vs 1D energy curvature at the
    # ATV-convergent window (8 A / 150 Ry); offline CP2K 10 A/400 Ry pin
    # stays in REFERENCE for a future mixer-overlap twin.
    if(ATV_BUILD_VIB)
        atv_add_regression(gpw gpw-fd-phonon-si)
    endif()
    atv_add_regression(gpw gpw-variable-cell-ice)
    atv_add_regression(gpw gpw-variable-cell-molopt)
    atv_add_regression(gpw periodic-mp2-si)
    atv_add_regression(gpw periodic-ump2-si)
    atv_add_regression(gpw periodic-ri-mp2-si)
    atv_add_regression(gpw periodic-mp2-kpoints)
    atv_add_regression(gpw periodic-ump2-kpoints)
    atv_add_regression(gpw periodic-rpa-si)
    atv_add_regression(gpw periodic-rpa-kpoints)
    atv_add_regression(gpw periodic-rpa-finite-q)
    atv_add_regression(gpw periodic-urpa-si)
    atv_add_regression(gpw periodic-urpa-kpoints)
    atv_add_regression(gpw periodic-gw-si)
    atv_add_regression(gpw periodic-ugw-si)
    atv_add_regression(gpw periodic-ugw-kpoints)
    # UGW@UKS RPBE spin XC residual (spectral G0W0; xc_grid_rpbe_spin).
    atv_add_regression(gpw ugw-uks-rpbe)
    atv_add_regression(gpw periodic-gw-kpoints)
    atv_add_regression(gpw periodic-gw-finite-q)
    # MPI k-point RR industrial gate (self-launches mpiexec when ATV_USE_MPI).
    if(ATV_USE_MPI AND MPIEXEC_EXECUTABLE)
        atv_add_regression(gpw mpi-gpw-scaling-4rank)
        target_compile_definitions(reg_mpi-gpw-scaling-4rank PRIVATE
            ATV_MPIEXEC="${MPIEXEC_EXECUTABLE}")
        # Parent self-launches mpiexec; link MPI for the --mpi-worker path.
        target_link_libraries(reg_mpi-gpw-scaling-4rank PRIVATE MPI::MPI_CXX)
    endif()
    # remaining green GPW e2e land with their branches
endif()
if(ATV_BUILD_HF)
    # Layout: SCF API hygiene — MolecularGrid ctor links + rys_roots throws
    # catchably (scf-review-05-api-hygiene). Needs the SCF headers, hence HF.
    atv_add_regression(layout scf_api_hygiene)
    atv_add_regression(scf h2o_rhf_sto3g)
    atv_add_regression(scf h2o_rhf_def2svp_sph)
    atv_add_regression(scf h2o_rhf_def2svp_grad)
    atv_add_regression(scf h2o_rhf_dipole)
    atv_add_regression(scf wannier-centers)
    atv_add_regression(scf static-polariz-h2o)
    atv_add_regression(scf h2o_embedding_field)
    atv_add_regression(scf embed_h2o_grad)
    atv_add_regression(scf field_h2o_grad)
    atv_add_regression(scf h2o_rhf_631gs_cpcm)
    atv_add_regression(scf ief-pcm-h2o)
    atv_add_regression(scf pcm-grad-h2o)
    atv_add_regression(scf h2o_charges_pyscf)
    atv_add_regression(scf n2_ionization)
    atv_add_regression(scf diffuse_anion)
    atv_add_regression(scf ecp_heavy_atom)
    atv_add_regression(scf ecp_i2_grad)
    atv_add_regression(scf h2o_scf_convergence)
    atv_add_regression(scf h2_rhf_stability)
    atv_add_regression(scf soscf_hard_case)
    atv_add_regression(scf o2_rohf_sto3g)
    atv_add_regression(scf o2_rohf_grad)
    atv_add_regression(scf roks-o2-sto3g)
    atv_add_regression(scf roks-grad-o2-sto3g)
    atv_add_regression(scf h2o_ri_jk)
    atv_add_regression(scf c10h22_cfmm_j)
    atv_add_regression(scf h3_ghf_sto3g)
    atv_add_regression(scf ghf-grad-h3-sto3g)
    atv_add_regression(scf h2o_rhf_freq)
    atv_add_regression(scf uhf-hessian-o2)
    atv_add_regression(scf rohf-hessian-o2)
    atv_add_regression(scf rks-hessian-h2o-sto3g)
    atv_add_regression(scf uks-hessian-o2-sto3g)
    atv_add_regression(scf h2o_rhf_ir)
    atv_add_regression(scf ch2o_tdhf)
    atv_add_regression(scf cisd-h2co-sto3g)
    atv_add_regression(scf adc2-h2co-sto3g)
    atv_add_regression(scf adc2x-h2co-sto3g)
    atv_add_regression(scf g0w0-h2-sto3g)
    atv_add_regression(scf direct-rpa-h2-sto3g)
    atv_add_regression(scf scf-review-01-linalg)
    atv_add_regression(scf gf2-h2-sto3g)
    atv_add_regression(scf bse-h2-sto3g)
    atv_add_regression(scf h2co_rks_tda)
    atv_add_regression(scf h2_rks_alda_tda)
    atv_add_regression(scf h2_full_grid_alda)
    atv_add_regression(scf h2-b3lyp-hybrid-tda)
    atv_add_regression(scf h2-b3lyp-rpa-tddft)
    # scf-review-06-scf-perf: restricted B3LYP/def2-SVP H2O single point through
    # the restricted-hybrid Fock build (the fused one-sweep jk_rhf path); pinned
    # PySCF 2.14.0 total energy, no live oracle at ctest time.
    atv_add_regression(scf b3lyp-fused-jk-h2o)
    atv_add_regression(scf uks-tda-o2-lsda)
    atv_add_regression(scf uks-full-lr-o2-lsda)
    atv_add_regression(scf spin-flip-tda-o2)
    atv_add_regression(scf sf-tddft-o2)
    atv_add_regression(scf delta-scf-n2)
    atv_add_regression(scf full-lr-tddft-h2-pbe)
    atv_add_regression(scf blyp-water-sto3g)
    atv_add_regression(scf b2plyp-h2o)
    atv_add_regression(scf d3bj-molecular-pbe)
    atv_add_regression(scf dft-d4-h2o-pbe)
    atv_add_regression(scf h2_analytic_tda_force)
    atv_add_regression(scf h2_tda_force)
    atv_add_regression(scf h2_tda_nac)
    atv_add_regression(scf h2_nac_vectors)
    atv_add_regression(scf h2_alda_vs_cis)
    atv_add_regression(scf h2o_mp2_sto3g)
    # scf-review-07-posthf-perf: H2O/STO-3G plain RMP2 end-to-end through the
    # block-assembled dense ERI tensor (IntegralPack::eri_dense); pinned PySCF
    # 2.14.0 e_hf/e_corr, no live oracle at ctest time.
    atv_add_regression(scf posthf-eri-block-mp2)
    atv_add_regression(scf ccsd-h2-sto3g)
    atv_add_regression(scf ccsd-t-h2o-sto3g)
    atv_add_regression(scf mp4-h2-sto3g)
    atv_add_regression(scf cc2-h2-sto3g)
    atv_add_regression(scf casci-h2-sto3g)
    atv_add_regression(scf ci-engine-equiv-h2-sto3g)
    # scf-review-08-ccsd-spatial-diis: closed-shell CCSD spatial-vs-spin-orbital
    # E_corr parity + DIIS-vs-Jacobi fixed point on H2O/STO-3G, vs PySCF golden.
    atv_add_regression(scf ccsd-spatial-so-parity)
    atv_add_regression(scf casscf-h2-sto3g)
    atv_add_regression(scf nevpt2-h2-sto3g)
    atv_add_regression(scf eom-ee-ccsd)
    atv_add_regression(scf ea-adc2-h2-sto3g)
    atv_add_regression(scf ip-adc2-h2o-sto3g)
    atv_add_regression(scf eom-cc2-h2)
    atv_add_regression(scf eom-ea-ccsd)
    atv_add_regression(scf eom-ip-ccsd)
    atv_add_regression(scf eom-sf-ccsd-o2)
    atv_add_regression(scf h2o_ri_mp2)
    atv_add_regression(scf ump2-o2-sto3g)
    atv_add_regression(scf romp2-o2-sto3g)
    atv_add_regression(scf frozen-core-mp2-h2o)
    atv_add_regression(scf scs-mp2-h2o-sto3g)
    atv_add_regression(scf sos-mp2-h2o-sto3g)
    atv_add_regression(scf mp2-grad-h2o-sto3g)
    atv_add_regression(scf ccsd-grad-h2o-sto3g)
    atv_add_regression(scf ump2-grad-o2-sto3g)
    atv_add_regression(scf mp3-h2-sto3g)
    atv_add_regression(scf h2o_dimer_cp_hf631gs)
    atv_add_regression(scf x2c_argon)
    # AO-basis Car–Parrinello MD (dynamics domain, esmm-cpmd-04): fictitious
    # electronic DOF on the Gaussian SCF methods. Needs the MD driver as well
    # as HF — the electrons are propagated on the integrator's clock, so these
    # only build when both product modules are on.
    if(ATV_BUILD_MD)
        # RHF/STO-3G H2 — seed / zero-iteration production / kConserved wiring.
        atv_add_regression(dynamics cpmd-rhf-h2)
        # RKS(LSDA)/STO-3G H2O — fixed-density KS branch + XC gradient.
        atv_add_regression(dynamics cpmd-rks-water)
        # UHF triplet O2/STO-3G — dual-channel (alpha + beta) fictitious DOF.
        atv_add_regression(dynamics cpmd-uhf-o2)
        # Extended-Lagrangian BOMD on the same seam (esmm-cpmd-05): the
        # auxiliary density is the SCF guess instead of the propagated
        # coefficients. RHF/STO-3G H2O — warm-step iteration reduction vs a
        # cold-start twin and a plain-BOMD (last-density reuse) twin.
        atv_add_regression(dynamics xl-bomd-rhf-water)
    endif()
endif()
if(ATV_BUILD_PATH)
    atv_add_regression(path neb_muller_brown)
    atv_add_regression(path remd_harmonic)
    atv_add_regression(path pimd_free_particle)
    atv_add_regression(path dimer_saddle)
    atv_add_regression(path string_muller_brown)
    atv_add_regression(path irc_muller_brown)
    atv_add_regression(path pimd_free_particle_estimators)
endif()
atv_add_regression(architecture harmonic_driver_matrix)
atv_add_regression(architecture module_compile_defs)
if(ATV_BUILD_MINIMIZER)
    atv_add_regression(minimizer lj_dimer)
    atv_add_regression(minimizer fire_lj_dimer)
    atv_add_regression(minimizer cell_lj_pressure)
    atv_add_regression(minimizer joint_lj_zero_p)
    atv_add_regression(minimizer aniso_lj_box)
    atv_add_regression(minimizer triclinic_lj_shear)
endif()
if(ATV_BUILD_POTENTIALS)
    atv_add_regression(forcefield tip3p_two_water)
    atv_add_regression(forcefield gaff_peo_fragment)
    atv_add_regression(forcefield lj_soft_dimer)
    atv_add_regression(forcefield coul_soft_dimer)
    atv_add_regression(forcefield yukawa_dimer)
    atv_add_regression(forcefield buckingham-pair-lammps)
    atv_add_regression(forcefield born-repulsion-foothold)
    atv_add_regression(forcefield fumi-tosi-nacl)
    atv_add_regression(forcefield steele-wall-argon)
endif()
if(ATV_BUILD_PATH)
    atv_add_regression(sampling umbrella_harmonic)
    atv_add_regression(sampling umbrella_angle)
    atv_add_regression(sampling umbrella_dihedral)
    atv_add_regression(sampling fep_soft_dimer)
    atv_add_regression(sampling metad_double_well)
    atv_add_regression(sampling metad_2d)
    atv_add_regression(sampling wham_double_well)
    atv_add_regression(path hremd_soft_lambda)
endif()
if(ATV_BUILD_MD)
    atv_add_regression(md brownian_einstein)
    atv_add_regression(md settle_tip3p)
    atv_add_regression(md fene_wca_chain)
    atv_add_regression(md dpd_groot_warren)
    atv_add_regression(md dpd_ab_mix)
endif()
if(ATV_BUILD_SEMIEMPIRICAL)
    atv_add_regression(semiempirical am1_water_force)
    atv_add_regression(semiempirical am1_water_minimize)
    atv_add_regression(semiempirical am1_water_bomd)
    atv_add_regression(semiempirical rm1_water_energy)
    atv_add_regression(semiempirical h2_dftb2)
    atv_add_regression(semiempirical h2o_dftb2_mio)
    atv_add_regression(semiempirical dftb2-3ob-organic)
    atv_add_regression(semiempirical dftb2-analytic-forces)
    atv_add_regression(semiempirical h2o_gfn0)
    atv_add_regression(semiempirical spin-dftb2-o2)
    atv_add_regression(multiscale qmmm_h2_ar)
    atv_add_regression(multiscale qmmm_es_h2)
    atv_add_regression(multiscale qmmm_link_ch3oh)
    atv_add_regression(multiscale qmmm_pol_h2o)
endif()
if(ATV_BUILD_KSPACE AND ATV_BUILD_MD)
    atv_add_regression(md npt_pme_nacl)
    atv_add_regression(md slab-2d-ewald-foothold)
    atv_add_regression(md const-pot-md-nvt-short)
    atv_add_regression(md global-electroneutrality-bo)
    atv_add_regression(md electrode-atom-forces)
    atv_add_regression(md maze-charge-dynamics)
    atv_add_regression(md gaussian-coul-neighborlist)
    atv_add_regression(md born-buck-neighborlist-bulk)
    atv_add_regression(md displacement-dfield-md)
    atv_add_regression(md const-pot-tip4p-edl)
    atv_add_regression(md charge-charge-autocorr)
    atv_add_regression(md elastance-condition-observer)
    atv_add_regression(md constant-velocity-smd)
endif()
# FMM operator math is always-on (facade longrange/); no product-module gate.
atv_add_regression(longrange fmm_point_charges)

# ── pure-Python import (python-on-demand-install-02) ──────────────────────
# No C++ link: isolated pure sources only. Needs a Python interpreter.
find_package(Python3 COMPONENTS Interpreter QUIET)
if(Python3_Interpreter_FOUND)
    add_custom_target(run_reg_python_pure_import
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/python/pure_import/reg_pure_import.sh
                ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/pure_import
        COMMENT "python pure_import: soft import without C++ extension"
        VERBATIM)
    add_dependencies(regressions run_reg_python_pure_import)

    add_custom_target(run_reg_python_install_cli
        COMMAND bash ${CMAKE_CURRENT_SOURCE_DIR}/python/install_cli/reg_install_cli.sh
                ${Python3_EXECUTABLE}
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/install_cli
        COMMENT "python install_cli: help / mutual exclusion / hf / version lock"
        VERBATIM)
    add_dependencies(regressions run_reg_python_install_cli)
endif()

# ── launcher parity (mpi-launcher-embedded-python-04) ─────────────────────
# Not an atv_add_regression() case: it needs mpirun + a Python environment
# with atomiverse/molrs importable (the ~/work/.x86_64 venv per the spec).
# The launcher binary and interpreter come from the configured build.
if(ATV_BUILD_LAUNCHER)
    # regressions/ may be processed before app/atomiverse/ — resolve the
    # interpreter here too (cached, so this is a no-op when already found).
    find_package(Python3 COMPONENTS Interpreter QUIET)
    add_custom_target(run_reg_launcher_lj_parity
        COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/launcher/lj_parity/reg_lj_parity.sh
                ${Python3_EXECUTABLE}
                $<TARGET_FILE:atomiverse_launcher>
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/launcher/lj_parity
        COMMENT "launcher lj_parity: serial-python vs 2-rank launcher (LJ + AM1)"
        VERBATIM)
    add_dependencies(run_reg_launcher_lj_parity atomiverse_launcher)
    add_dependencies(regressions run_reg_launcher_lj_parity)
endif()
