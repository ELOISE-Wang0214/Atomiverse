// cuda-kokkos-loop — classical CUDA and the Kokkos backend.
// Runs beside hf-pyscf-loop. That loop owns Gaussian HF, including cuda/scf.
// This loop never edits those files. One slug per wave; device sources collide.

let meta = #{
    name: "cuda-kokkos-loop",
    description: "Advance classical CUDA and the Kokkos backend, one slug per wave, beside the HF loop",
    when_to_use: "Parallel with hf-pyscf-loop: CUDA MD/neighbor/longrange gaps and a real Kokkos backend",
    phases: [
        #{ title: "Census", detail: "one CUDA or Kokkos gap" },
        #{ title: "Implement", detail: "worktree wf/ck-*" },
        #{ title: "Test", detail: "scoped unit tests" },
        #{ title: "Land", detail: "merge into dev after tests pass" },
    ],
};

let mode = "burst";
let max_waves = 1;
let target_branch = "dev";
let slug_arg = "";
if args != () {
    if args.mode != () { mode = args.mode; }
    if args.target_branch != () { target_branch = args.target_branch; }
    if args.slug != () { slug_arg = args.slug; }
    if args.max_waves != () { max_waves = args.max_waves; }
}
if mode != "exhaust" && mode != "burst" { mode = "burst"; }
if mode == "exhaust" && (args == () || args.max_waves == ()) { max_waves = 4; }
if max_waves < 1 { max_waves = 1; }
if max_waves > 8 { max_waves = 8; }

let census_schema = #{
    "type": "object",
    "required": ["slug", "kind", "summary"],
    "properties": #{
        "slug": #{ "type": "string" },
        "kind": #{ "type": "string" },
        "summary": #{ "type": "string" },
        "classic_systems": #{ "type": "string" },
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
law += "Read CLAUDE.md before editing. CUDA classical MD is the OpenMM cluster/tile engine. Do not add a per-atom Verlet production force path. ";
law += "FILE FENCE: do not create or edit files under cpu/scf, include/atomiverse/cpu/scf, cuda/scf, include/atomiverse/cuda/scf, cpu/gpw, cuda/gpw, cpu/planewave, or cuda/planewave. ";
law += "Those trees belong to hf-pyscf-loop and the CP2K loops. Before the first edit, run git worktree list. If any wf/hf-* worktree diffs a file you need, status=blocked and name the file. Do not edit it. ";
law += "Do not edit the live x2c branch wf/hf-hf-x2c-forces. ";
law += "KOKKOS: there is no kokkos backend yet. docs/roadmap/index.md says Kokkos is planned and is not a table column. A Kokkos slug must add a compiling target, not a comment. ";
law += "Reference for parallel structure is the LAMMPS checkout at ../../lammps-src. Read the Kokkos pair or comm file you port and stamp that path in the provenance header. Do not invent a second communication pattern. ";
law += "CUDA compiles only on the GH200 via sbatch tools/gpu_job.sh. Never #SBATCH --exclusive. Never request a full node. Login nodes have no CUDA toolkit. ";
law += "A CPU fallback does not close a CUDA cell. If the kernel does not exist, throw and name it. ";
law += "TESTS: tests/ are single-function unit tests. End-to-end goes to regressions/<domain>/<slug>/ as an out-of-library consumer. ";
law += "Build only the targets this slug touches. Never bare cmake --build of the whole tree. Never unfiltered ctest. ";
law += "CPU configure: bash -lc with module load buildenv-gcc/2026.03-mpich, OPENBLAS_ROOT, LD_LIBRARY_PATH mirror, ATV_BUILD_CUDA=OFF, ATV_BUILD_HF=OFF, ATV_BUILD_TESTS=ON. ";
law += "LAND GATE: merge into dev as soon as the scoped unit tests pass. Do not wait for a GPU job, a simplify pass, or a tutorial. ";
law += "Update the matching roadmap row in the same commit as the code when the public symbol changed. Do not write wave ids or job ids into docs/. ";

fn stage_ok(r) {
    if r == () || r.success != true || r.output == () || r.output.status == () {
        return false;
    }
    let s = r.output.status;
    s == "green" || s == "pass" || s == "ok" || s == "merged" || s == "success" || s == "done"
}

let wave = 0;
let landed = [];
let last_slug = "";

while wave < max_waves {
    wave += 1;
    phase("Census");
    let census_prompt = law;
    census_prompt += "You are GapCensus for classical CUDA and Kokkos, wave " + wave.to_string() + ". ";
    census_prompt += "Use grep and read_file. Read docs/roadmap/index.md, docs/roadmap/molecular-dynamics.md, docs/roadmap/force-fields.md, and list src/cuda and include/atomiverse. ";
    census_prompt += "If no directory kokkos/ exists under src/ or include/atomiverse/, the slug must be kokkos-foothold and kind=kokkos: one ATV_BUILD_KOKKOS option, one CPU-execution-space kernel whose result matches an existing CPU LJ or axpy unit test, and a header that names the LAMMPS Kokkos file you read. Do not vendor all of Kokkos if find_package fails; status stays for the implementer. ";
    census_prompt += "Otherwise pick ONE open CUDA cell in molecular-dynamics.md or force-fields.md whose sources are under src/cuda/dynamics, neighbor, potentials, longrange, communication, device, frame, or minimize. kind=cuda. ";
    census_prompt += "Do not pick a Gaussian HF, GPW, or plane-wave slug. Do not rename a closed slug. ";
    if slug_arg != "" {
        census_prompt += "The operator forced slug " + slug_arg + ". Use that slug if it is inside the file fence. ";
    }
    census_prompt += "Return slug, kind (cuda or kokkos), classic_systems, summary, gaps_remaining.";

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
    let kind = "cuda";
    if census.output.kind != () { kind = census.output.kind; }
    let classic = "a small LJ box or one pair";
    if census.output.classic_systems != () { classic = census.output.classic_systems; }
    log("wave " + wave.to_string() + " " + kind + " " + slug);

    phase("Implement");
    let impl_prompt = law;
    impl_prompt += "Implement slug " + slug + " (kind " + kind + ") on a new branch wf/ck-" + slug + " from " + target_branch + ". ";
    impl_prompt += "Work in the isolation worktree. Stay inside the file fence. If the change needs a forbidden path, status=blocked and do not commit that edit. ";
    impl_prompt += "Commit on the worktree branch. status=green only if the CPU configure builds the targets you added. ";
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
        let drop_prompt = "Remove the isolation worktree and branch wf/ck-" + slug + " if they exist. Do not reset " + target_branch + ". Do not delete wf/hf-* worktrees. Wipe only that worktree's build directory. Return status and summary.";
        agent(drop_prompt, #{ label: "drop:" + slug, capability_mode: "execute", output_schema: stage_schema });
        if mode == "burst" {
            complete(#{ status: "blocked", slug: slug, summary: why, landed: landed });
        }
    } else {
        let tree = impl.output.worktree_path;
        if tree == () { tree = ""; }
        let branch = impl.output.branch;
        if branch == () { branch = "wf/ck-" + slug; }

        phase("Test");
        let test_prompt = law;
        test_prompt += "In worktree " + tree + " branch " + branch + ", run scoped CPU unit tests for slug " + slug + ". ";
        test_prompt += "cd to that worktree for every command. Do not register an end-to-end MD run under tests/. ";
        test_prompt += "status=green only if those ctests passed. Return status, summary, worktree_path, branch, head_sha.";
        let tested = agent(test_prompt, #{
            label: "test:" + slug,
            capability_mode: "execute",
            output_schema: stage_schema,
        });

        let ok = stage_ok(tested);

        phase("Land");
        if ok {
            let land_prompt = law;
            land_prompt += "Tests passed for slug " + slug + " on branch " + branch + " in worktree " + tree + ". Merge that branch into " + target_branch + " now. ";
            land_prompt += "The parent checkout of dev is dirty with unrelated molrs and I/O edits. Do not stash, reset, or overwrite those files. ";
            land_prompt += "Make the merge commit in a detached worktree: git worktree add --detach /tmp/ck-land-" + slug + " " + target_branch + ", then git merge --no-ff " + branch + " there. ";
            land_prompt += "Then, in the parent repo, git merge --ff-only that merge commit. If ff-only would touch a dirty file, status=blocked and leave both branches. ";
            land_prompt += "No force-push. If the merge conflicts, status=blocked and leave " + branch + ". ";
            land_prompt += "On success, delete worktree " + tree + ", /tmp/ck-land-" + slug + ", and the merged local branch. Remove only those worktrees' build directories. ";
            land_prompt += "Do not delete wf/hf-* worktrees or wf/ck-kokkos-foothold unless that is the branch you just merged. Do not wipe the parent build, build-py, or .cache. Return status and summary.";
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
            let drop_prompt = "Unit tests failed for " + slug + ". Remove worktree " + tree + " only after deleting its build/. Leave branch " + branch + " if it has commits. Do not reset " + target_branch + ". Do not delete wf/hf-* worktrees. Do not wipe the parent build. Return status and summary.";
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
    summary: "cuda-kokkos-loop finished",
});
