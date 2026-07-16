---
slug: frame-block-04-delete-compatibility
criteria:
  - id: ac-001
    summary: CPU Frame physically owns only blocks, typed meta, and simbox
    type: code
    pass_when: |
      cpu::Frame's storage is exactly backend-local BlockMap + typed MetaMap +
      simbox. Frame has no atoms, scalars, or box member; meta entries do not
      proxy another aggregate; Columns own the xtensor buffers. cpu::Atom is no
      longer a Frame storage container and there is no dual ownership. Mutable
      typed meta references cover counts, all energies, and both stress arrays.
    status: verified
    last_checked: 2026-07-16
    notes: |
      cpu::Frame has exactly BlockMap blocks, MetaMap meta, and Box simbox as
      physical state. cpu::Atom and SystemScalars storage are deleted; typed
      counts, energies, stress arrays, and the backend-local tag index live in
      MetaMap, while Column buffers own per-row storage.

  - id: ac-002
    summary: old Frame and Box source surfaces are deleted
    type: code
    pass_when: |
      A source gate and negative compile fixtures prove that frame.atoms,
      frame.scalars, frame.box, ATV_FRAME_META_*_VIEWS, ATV_FIELD_BIND,
      ExtraColumn, the Frame-to-Atom insert forwarder, Box::operator()(i,j),
      Box::operator=(Mat3), and every inventoried migration accessor/alias are
      absent. No deprecated declaration, forwarding replacement, or const_cast
      resolver exists.
    status: pending

  - id: ac-003
    summary: every in-tree CPU consumer uses typed mapping handles
    type: code
    pass_when: |
      All real cpu::Frame consumers in CPU production headers/sources, app,
      Python bindings, tests, regressions, examples, and docs use typed
      span/meta handles and simbox. Build may cache logical slots but not raw
      handles across structural mutation; raw handles resolve/refresh at the
      safe launch boundary. No particle, pair, grid, or topology inner loop
      resolves a key.
    status: verified
    last_checked: 2026-07-16
    notes: |
      The CUDA-off CPU tree was migrated to typed span/meta handles and
      simbox. compatibility_surface_scrub.sh, direct_access_scrub.sh, and
      field_list_scrub.sh pass; the resolution-window unit deliberately puts a
      lookup in an inner loop and proves the Driver rejects it.

  - id: ac-004
    summary: the registry is the sole storage and movement schema
    type: runtime
    pass_when: |
      Block construction, capacity change, copy_row, memory accounting,
      transaction rollback, exchange/halo/migration packing, materialize, and
      write_back all derive from core/frame_schema.h semantic tags. Tests add
      a runtime column and prove it follows those operations without editing a
      second field list or entering a separate ExtraColumn path. Resolved
      column/pack plans contain no string/hash/variant dispatch in per-row
      loops and are invalidated/rebuilt after capacity or topology re-stride.
    status: pending

  - id: ac-005
    summary: topology and Frame round-trip preserve the current contract
    type: runtime
    pass_when: |
      Atomic and bonded Frames round-trip x/v/f/charge/image/tag and the current
      write-only scalars/stress contract. CPU per-atom bond/angle/dihedral
      columns retain runtime strides; the molrs flat conversion occurs only at
      materialize/write_back. Tag narrowing still throws on overflow.
    status: verified
    last_checked: 2026-07-16
    notes: |
      CPU materialize/write_back and MPI Frame tests cover atomic and bonded
      topology, x/v/f/charge/image/tag, write-only scalars/stress, flat/per-atom
      topology conversion, and checked i32 tag narrowing.

  - id: ac-006
    summary: typed meta is one exact cross-repository contract
    type: runtime
    pass_when: |
      molrs core, views, serializers, C/Python/CXX APIs, Atomiverse, and molpy
      consume one exact typed-meta semantic schema and API/format version.
      Exchangeable kind/wire IDs match; the Atomiverse tag index is a typed
      local entry that is provably never serialized.
      HashMap<String,String>, typed-to-string downgrade, old-format decoder,
      permissive version warning, and alternate metadata payloads are absent;
      every non-current version fails before data access.
    status: pending

  - id: ac-007
    summary: compatibility-only CPU and Python APIs are gone
    type: code
    pass_when: |
      The checked inventory has zero remaining CPU declarations whose purpose
      is preserving an old spelling, namespace, overload, default behavior,
      import path, ABI, or payload. Known aliases and transient-solver shims are
      deleted; molpy.core.frame and molpy Frame/Block re-exports do not import.
      RHF/UHF/UKS adapter results, the implicit grid overload, UPF diagonal
      cache, Mat33, CPU flat core types/constants/log overloads, and TaskSpec
      variant-order IDs are absent. In-tree callers use only selected current
      APIs and explicit current ComponentId values.
    status: verified
    last_checked: 2026-07-16
    notes: |
      The checked Atomiverse compatibility inventory passes and pins deletion
      of Frame/Box aliases, solver shims, result adapters, flat numeric headers,
      old constants/overloads, and variant-order TaskSpec IDs. Molpy tests pin
      that molpy.core.frame and facade Frame/Block/Element imports are absent.

  - id: ac-008
    summary: old inputs fail rather than translate
    type: runtime
    pass_when: |
      Negative tests submit old Frame/meta keys and dtypes, old TaskSpec/Frame
      versions, old Python imports, and representative old C++ spellings. Each
      fails at compile/import/decode/validation time; instrumentation proves no
      fallback, converter, dual read, or dual write executes.
    status: pending

  - id: ac-009
    summary: native-dtype wire and canonical MPI schema are bit exact
    type: runtime
    pass_when: |
      Pack/broadcast/border/exchange/migration/gather preserve StorageReal,
      AccReal, int8/int32, and 32/64-bit TagInt without a double conversion
      funnel; 32-bit extrema and bigbig TagInt 2^53+1 round-trip bit-exactly.
      Ranks collectively compare a fingerprint of version plus sorted
      key/dtype/stride/flags before payload.
      Equivalent runtime columns inserted in different order produce identical
      wire; a true mismatch fails all ranks without deadlock. Registration
      rejects any column limit before pack; no payload is silently truncated.
    status: pending

  - id: ac-010
    summary: structural mutation cannot leave a stale raw handle
    type: runtime
    pass_when: |
      Build caches logical slots only. ASan tests trigger capacity growth,
      topology re-stride, migration, and runtime column/meta insertion between
      two launches. Storage generation invalidates every affected plan, and an
      old Span/MetaRef is refreshed or rejected before dereference; no UAF or
      stale row/stride result is possible.
    status: verified
    last_checked: 2026-07-16
    notes: |
      test_frame_registry: Span stales on grow/column-insert/topology restride;
      MetaRef stales on meta insert; operator[]/data() throw stale_handle_error.
      Full Linux ASan matrix (migration between launches) still desirable on CI
      but generation rejection is unit-gated and default-on.

  - id: ac-011
    summary: CPU precision ID MPI and hybrid-thread matrices are non-vacuous
    type: runtime
    pass_when: |
      CUDA-off frame/comm tests cover supported StorageReal float/double x
      WorkReal float/double x MPI off/on plus smallbig/bigbig ID modes; the
      double full-domain build enables RHF/PW/GPW/KSpace. OpenMP is REQUIRED,
      never a silent serial fallback. mpirun -n 2 with OMP_NUM_THREADS=1 and 4
      runs broadcast->partition->halo->migration->gather+LJ and matches
      tag-sorted state/energy/force; debug asserts MPI is called outside OpenMP
      parallel regions. Unit tests, examples, and public regressions pass.
    status: pending

  - id: ac-012
    summary: resolution is O(1) per operation and zero-cost in Release
    type: performance
    pass_when: |
      Every accepted hot component has an exact resolution manifest. Debug
      runs at N=1, N=64, and N>=100000, including nested replay/sub-driver,
      prove span/has/block/meta counts do not scale with atoms or pairs. A bad
      OpenMP component performs parallel lookups; the thread-safe counter is
      TSAN-clean and deterministically fails after leaving the parallel region.
      Release compile/object inspection proves begin/end/note add no hot-path
      work.
    status: pending

  - id: ac-013
    summary: affected classical steady paths allocate zero bytes
    type: performance
    pass_when: |
      A mandatory Linux acceptance runner intercepts thread-safe new/new[],
      malloc/calloc/realloc, aligned_alloc, and posix_memalign and reports count
      plus bytes. After warmup and within reserved capacity, one/many-thread
      LJ+neighbor+integrator reuse and forced-rebuild paths plus two-rank halo
      refresh/rebuild/migration allocate zero. First build/grow is a separate
      ledger. SCF/PW/GPW are not falsely declared allocation-free; the Frame
      cutover adds zero count/bytes versus their baselines.
    status: pending

  - id: ac-014
    summary: each reproducible LJ sentinel regresses by no more than 2 percent
    type: performance
    pass_when: |
      Pre-cutover baselines pin git SHA, host/CPU/ISA, toolchain/flags, N,
      xtensor/xsimd version and batch width, neighbor cadence, warmup, thread
      count/affinity, and BLAS threads=1. Noise is demonstrated below 1%; at
      least ten samples use the median. Neighbor-reuse and forced-reneighbor
      single-rank cases each fail independently above 2%; no geomean, 8%
      fallback, warning-only result, or missing-baseline skip can pass.
    status: pending

  - id: ac-015
    summary: Frame and communication lifecycle stages have independent gates
    type: performance
    pass_when: |
      Atomic/bonded materialize and write_back report ns/atom, GB/s,
      allocations, and peak bytes separately. Pack/unpack, broadcast, steady
      halo, halo rebuild, forced migration, and gather are separate 1/2-rank
      cases; plan rebuild is not averaged into steady work. Every case has a
      same-runner baseline and fails above 5%. The memory ledger separates
      blocks/meta/simbox, capacity/used, comm scratch, and source/destination
      and proves no second per-row mirror; RSS alone is insufficient.
    status: pending

  - id: ac-016
    summary: compiler SIMD and SCF sentinels retain their current performance
    type: performance
    pass_when: |
      compile_commands proves optimized CPU CXX TUs retain private,
      COMPILE_LANGUAGE:CXX -march=native and no target defines
      XTENSOR_USE_XSIMD. The exact validated xtensor/xsimd pair is pinned and
      reported; scalar/SIMD LJ agree, SIMD batches execute, and XC reshape_view
      passes. BM_GBuild_Block_Aspirin is <=950 ms and <=5% above baseline;
      BM_GBuild_Block_Aspirin_ccPVDZ is <=5% above baseline; the Rys d-shell
      cross-check stays within 1e-7.
    status: pending

  - id: ac-017
    summary: docs and releases state a hard break with no migration path
    type: docs
    pass_when: |
      Atomiverse, molrs, and molpy docs/changelogs show only the new API,
      identify removed symbols and the required exact versions, and explicitly
      state that old callers/data must be rewritten or regenerated. They do not
      publish a compatibility table, adapter, converter, or deprecation period.
    status: pending

  - id: ac-018
    summary: GPU and Torch physical cutover
    type: code
    pass_when: |
      DEFERRED. A future campaign must move CUDA and Torch to physical mapping
      ownership, migrate their consumers, delete their compatibility code, and
      establish cross-backend scientific/performance parity without reviving a
      CPU compatibility surface.
    status: deferred
    last_checked: 2026-07-15
    notes: |
      User-directed GPU deferral; excluded from this CPU link's exit gate.
---

# Acceptance — frame-block-04-delete-compatibility

`ac-001` through `ac-017` are one indivisible CPU hard-cutover gate. Partial
migration, compatibility scaffolding, and accepting old inputs cannot be
reported as progress toward closure. `ac-018` remains explicitly deferred.
