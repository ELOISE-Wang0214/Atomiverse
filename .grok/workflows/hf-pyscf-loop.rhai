// hf-pyscf-loop — Gaussian HF debt-first, PySCF-sourced, CPU + CUDA-only.
// Not an instance of oracle-capability-loop-template: that template keeps
// CUDA out of scope and a different phase order. This script is the order
// the operator set: census → implement from PySCF → test → classic systems
// → simplify → manual, debt before new capability.
//
// One slug per wave. HF sources collide, so there is no parallel staff.

let meta = #{
    name: "hf-pyscf-loop",
    description: "HF gaps vs PySCF: debt first, then CPU + CUDA-only, test, classic systems, simplify, manual",
    when_to_use: "Close Gaussian HF debt and missing PySCF capabilities on CPU and CUDA together",
    phases: [
        #{ title: "Census", detail: "debt before missing capability" },
        #{ title: "Implement", detail: "read PySCF; CPU + CUDA-only" },
        #{ title: "Test", detail: "scoped unit tests" },
        #{ title: "Validate", detail: "classic molecules" },
        #{ title: "Simplify", detail: "diff-only cut" },
        #{ title: "Manual", detail: "roadmap + tutorial" },
        #{ title: "Land", detail: "merge green slug, wipe trees" },
    ],
};

let mode = "burst";
let max_waves = 1;
let target_branch = "dev";
let slug_arg = "";
if args != () {
    if args.mode != () { mode = args.mode; }
    if args.max_waves != () { max_waves = args.max_waves; }
    if args.target_branch != () { target_branch = args.target_branch; }
    if args.slug != () { slug_arg = args.slug; }
}
if mode != "exhaust" && mode != "burst" { mode = "burst"; }
if mode == "exhaust" && (args == () || args.max_waves == ()) { max_waves = 8; }
if max_waves < 1 { max_waves = 1; }
if max_waves > 20 { max_waves = 20; }

let census_schema = #{
    "type": "object",
    "required": ["slug", "kind", "summary"],
    "properties": #{
        "slug": #{ "type": "string" },
        "kind": #{ "type": "string" },
        "pyscf_files": #{ "type": "string" },
        "classic_systems": #{ "type": "string" },
        "summary": #{ "type": "string" },
        "gaps_remaining": #{ "type": "integer" },
    },
};
let stage_schema = #{
    "type": "object",
    "required": ["status", "summary"],
    "properties": #{
        "status": #{ "type": "string" },
        "summary": #{ "type": "string" },
        "worktree_path": #{ "type": "string" },
        "branch": #{ "type": "string" },
        "head_sha": #{ "type": "string" },
    },
};

let law = "";
law += "REPO: Atomiverse. Active branch is dev. Never commit to master. ";
law += "Read CLAUDE.md and .claude/notes/es-architecture-debt.md before editing. ";
law += "GOLDEN SOURCE: implement by reading the named PySCF files on disk. ";
law += "Do not re-derive the algorithm from memory. Pin numbers come from tools/pyscf_oracle (offline), never from ATV's own measurement. ";
law += "DEBT FIRST: a wave closes one architecture-debt slug before any new capability. ";
law += "Do not add a method, basis, or CUDA kernel on top of an open red row in es-architecture-debt.md that this slug touches. ";
law += "SPLIT: one type, one job. MeanFieldMethod keeps SCF run/forces only. Post-HF, solvent, and response stay out of new methods on that class. ";
law += "CLOSED, do not re-open: es-cuda-rhf-driver (cuda::scf::RHF is cuda::Driver<scf::Phase>); the mean-field split (ADC, DirectRPA, BSE, GF2, G0W0, CISD, TDA/TDHF, Stability, CCSD EOM-EE/IP/EA, CC2 EOM, CCSD EOM-SF, SfTDA/SfTDHF, Hessian, Polarizability, WithSolvent, MP2 are separate types). ghost_atoms stays a method member. ";
law += "OUT OF SCOPE even if the ledger names HF: es-kind-model, es-kind-name-hf, es-bind-dialects, es-se-driver, es-cuda-pw-driver, es-cuda-gpw-frame-io, es-ghost-atoms. Closing es-kind-model means editing Semiempirical, PlaneWave, and GPW kinds. Do not pick those rows. ";
law += "CUDA-ONLY: after setup/build, the production step (Fock, density, DIIS, diagonalize, forces) runs on device. ";
law += "Host code may pack geometry, basis, and sizes once at setup. It must not compute ERIs, the Fock matrix, the density, or the nuclear forces for a path the device claims. ";
law += "No host-fill-then-upload fallback for a covered angular momentum. If the kernel does not exist, throw and name the missing kernel. Do not call cpu:: native integrals from the step loop. ";
law += "CUDA compiles only on the GH200 via sbatch tools/gpu_job.sh. Never #SBATCH --exclusive. Never request a full node. ";
law += "Login nodes have no CUDA toolkit: do not pretend a CPU run validates the device path. ";
law += "TESTS: tests/ are single-function unit tests. End-to-end classic molecules go to regressions/scf/<slug>/ as an out-of-library consumer. ";
law += "Build only atomiverse plus the test binaries and reg_<slug> this slug touches. Never bare cmake --build build. Never unfiltered ctest. ";
law += "CPU configure: bash -lc with module load buildenv-gcc/2026.03-mpich, OPENBLAS_ROOT, LD_LIBRARY_PATH mirror, ATV_BUILD_HF=ON, ATV_BUILD_CUDA=OFF, ATV_BUILD_TESTS=ON, ATV_BUILD_REGRESSIONS=ON. ";
law += "MANUAL: update the matching row in docs/roadmap/hartree-fock.md. Add or adjust one tutorial page under docs/tutorial/ only when the public call shape changed. ";
law += "Do not write wave ids, slugs, or job ids into docs/. ";

fn stage_ok(r) {
    if r == () || r.success != true || r.output == () || r.output.status == () {
        return false;
    }
    let s = r.output.status;
    s == "green" || s == "pass" || s == "ok" || s == "merged" || s == "success"
}

let pyscf_root = "/nobackup/proj/disk/teoroo/personal/jicli594/work/molcrafts";

let wave = 0;
let landed = [];
let last_slug = "";

while wave < max_waves {
    wave += 1;
    phase("Census");
    let census_prompt = law;
    census_prompt += "You are GapCensus for Gaussian HF, wave " + wave.to_string() + ". ";
    census_prompt += "Use grep and read_file. Read docs/roadmap/hartree-fock.md, .claude/notes/es-architecture-debt.md, include/atomiverse/cpu/scf/rhf.h, include/atomiverse/cuda/scf/rhf.cuh. ";
    census_prompt += "Find PySCF sources with a filesystem search under " + pyscf_root + " (pyscf/scf, pyscf/df, pyscf/hessian). Name the files you actually opened in pyscf_files. ";
    census_prompt += "Pick ONE open Gaussian-HF capability from docs/roadmap/hartree-fock.md (a ○ cell: k-point HF, CUDA l>1, UHF on device, incremental Fock, x2c forces, ROHF CPHF, complex/spinor HF). ";
census_prompt += "Do not pick es-kind-model, es-kind-name-hf, es-ghost-atoms, es-bind-dialects, es-se-driver, es-cuda-pw-driver, es-cuda-gpw-frame-io, es-cuda-rhf-driver, es-mean-field-split, or mean-field-method-split. ";
census_prompt += "Do not rename a closed slug. kind=capability unless you find a red debt row whose files are only under cpu/scf, cuda/scf, or docs/roadmap/hartree-fock.md. ";
    census_prompt += "A foothold kernel is not a closed capability. gaps_remaining is the count of open HF debt rows plus roadmap ○ cells you saw. ";
    if slug_arg != "" {
        census_prompt += "The operator forced slug " + slug_arg + ". Use that slug. ";
    }
    census_prompt += "Return slug, kind (debt or capability), pyscf_files, classic_systems (molecules that exercise this slug, e.g. H2 STO-3G and H2O 6-31G), summary, gaps_remaining.";

    let census = agent(census_prompt, #{
        label: "census:" + wave.to_string(),
        capability_mode: "read-only",
        output_schema: census_schema,
    });
    if census == () || census.success != true || census.output == () || census.output.slug == () || census.output.slug == "" {
        log("census failed on wave " + wave.to_string());
        complete(#{ status: "blocked", wave: wave, landed: landed, summary: "census failed" });
    }
    let slug = census.output.slug;
    if slug == last_slug && slug_arg == "" {
        log("census repeated " + slug + "; stop");
        complete(#{ status: "stall", wave: wave, landed: landed, summary: "same slug twice" });
    }
    last_slug = slug;
    let kind = "capability";
    if census.output.kind != () { kind = census.output.kind; }
    let pyscf_files = "";
    if census.output.pyscf_files != () { pyscf_files = census.output.pyscf_files; }
    let classic = "H2 STO-3G, H2O";
    if census.output.classic_systems != () { classic = census.output.classic_systems; }
    log("wave " + wave.to_string() + " " + kind + " " + slug);

    phase("Implement");
    let impl_prompt = law;
    impl_prompt += "Implement slug " + slug + " (kind " + kind + ") on a new branch wf/hf-" + slug + " from " + target_branch + ". ";
    impl_prompt += "Work in the isolation worktree. Read these PySCF files first and cite the functions you transplanted: " + pyscf_files + ". ";
    impl_prompt += "CPU and CUDA both change when the slug is a product path. CUDA step loop is device-only; host work stays in setup. ";
    impl_prompt += "If kind=debt, refactor boundaries and delete the dialect. Do not add a new method in the same commit. ";
    impl_prompt += "Spec belongs in .claude/specs/ only if you need acceptance criteria; do not expand scope. ";
    impl_prompt += "Commit on the worktree branch. status=green only if the code compiles for the CPU configure. status=blocked if PySCF has no corresponding source or a named debt row forbids the change. ";
    impl_prompt += "Return status, summary, worktree_path, branch, head_sha.";
    let impl = agent(impl_prompt, #{
        label: "impl:" + slug,
        capability_mode: "execute",
        isolation_worktree: true,
        output_schema: stage_schema,
    });
    if !stage_ok(impl) {
        let why = "implement failed";
        if impl != () && impl.output != () && impl.output.summary != () { why = impl.output.summary; }
        log(why);
        phase("Land");
        let drop_prompt = "Remove the isolation worktree and branch wf/hf-" + slug + " if they exist. Do not reset " + target_branch + ". Do not delete unrelated worktrees. Wipe that worktree's build/ only. Return status and summary.";
        agent(drop_prompt, #{ label: "drop:" + slug, capability_mode: "execute", output_schema: stage_schema });
        if mode == "burst" {
            complete(#{ status: "blocked", slug: slug, summary: why, landed: landed });
        }
    } else {
        let tree = impl.output.worktree_path;
        if tree == () { tree = ""; }
        let branch = impl.output.branch;
        if branch == () { branch = "wf/hf-" + slug; }

        phase("Test");
        let test_prompt = law;
        test_prompt += "In worktree " + tree + " branch " + branch + ", run scoped CPU unit tests for slug " + slug + ". ";
        test_prompt += "cd to that worktree for every command. Add a failing unit test first if the slug is a bugfix; do not put an SCF solve in tests/. ";
        test_prompt += "status=green only if those ctests passed. Return status, summary, worktree_path, branch, head_sha.";
        let tested = agent(test_prompt, #{
            label: "test:" + slug,
            capability_mode: "execute",
            output_schema: stage_schema,
        });

        let ok = stage_ok(tested);
        if ok {
            phase("Validate");
            let val_prompt = law;
            val_prompt += "In worktree " + tree + ", validate slug " + slug + " on classic systems: " + classic + ". ";
            val_prompt += "Add or run regressions/scf/" + slug + " against an offline PySCF pin (tools/pyscf_oracle). ";
            val_prompt += "CPU energy or forces must match the pin within the tolerance already used for that observable. ";
            val_prompt += "Device path: submit tools/gpu_job.sh without --exclusive only when the slug touches CUDA. If this machine has no GPU, status=blocked and say cuda_not_run. Do not mark green from the host ERI path. ";
            val_prompt += "Return status, summary, worktree_path, branch, head_sha.";
            let val = agent(val_prompt, #{
                label: "validate:" + slug,
                capability_mode: "execute",
                output_schema: stage_schema,
            });
            ok = stage_ok(val);
        }

        if ok {
            phase("Simplify");
            let sim_prompt = law;
            sim_prompt += "In worktree " + tree + ", simplify only the diff of branch " + branch + " against " + target_branch + ". ";
            sim_prompt += "Delete duplicated host/device loops and pass-through wrappers. Keep public run/forces. Re-run the same scoped tests. ";
            sim_prompt += "status=green only if tests still pass. Return status, summary, worktree_path, branch, head_sha.";
            let sim = agent(sim_prompt, #{
                label: "simplify:" + slug,
                capability_mode: "execute",
                output_schema: stage_schema,
            });
            ok = stage_ok(sim);
        }

        if ok {
            phase("Manual");
            let man_prompt = law;
            man_prompt += "In worktree " + tree + ", document slug " + slug + ". ";
            man_prompt += "Update docs/roadmap/hartree-fock.md so the row matches the code you just read. ";
            man_prompt += "If the user-facing call changed, add one short tutorial section with a Frame in / Frame out example. ";
            man_prompt += "Commit on " + branch + ". status=green when the docs match the symbols in the headers. Return status, summary, branch, head_sha.";
            let man = agent(man_prompt, #{
                label: "manual:" + slug,
                capability_mode: "execute",
                output_schema: stage_schema,
            });
            ok = stage_ok(man);
        }

        phase("Land");
        if ok {
            let land_prompt = law;
            land_prompt += "Merge branch " + branch + " from worktree " + tree + " into " + target_branch + " with a merge commit. ";
            land_prompt += "No force-push. If the merge conflicts, status=blocked and leave the branch. ";
            land_prompt += "On success, delete that worktree and the merged local branch, and remove only that worktree's build directory. ";
            land_prompt += "Return status and summary.";
            let land = agent(land_prompt, #{
                label: "land:" + slug,
                capability_mode: "execute",
                output_schema: stage_schema,
            });
            if stage_ok(land) {
                landed.push(slug);
                log("landed " + slug);
            } else {
                log("merge failed for " + slug);
                if mode == "burst" {
                    complete(#{ status: "blocked", slug: slug, landed: landed, summary: "merge failed" });
                }
            }
        } else {
            let drop_prompt = "Tests or validation failed for " + slug + ". Remove worktree " + tree + " and branch " + branch + " only if the branch has no commit you were told to keep. Prefer leaving the branch and deleting the worktree's build/. Do not reset " + target_branch + ". Return status and summary.";
            agent(drop_prompt, #{ label: "drop:" + slug, capability_mode: "execute", output_schema: stage_schema });
            if mode == "burst" {
                complete(#{ status: "blocked", slug: slug, landed: landed, summary: "stage failed" });
            }
        }
    }
}

complete(#{
    status: "ok",
    waves: wave,
    landed: landed,
    summary: landed.len().to_string() + " slugs landed",
});
