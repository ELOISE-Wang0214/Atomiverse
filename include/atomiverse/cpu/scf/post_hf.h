/// @file cpu/scf/post_hf.h
/// @brief Post-HF correlation types. They are not methods of `MeanFieldMethod`.
///
/// Boundary transplanted from PySCF (not re-derived):
/// - `pyscf.mp.mp2.MP2(mf)` / `kernel` — correlation on a mean-field reference
/// - `pyscf.cc.ccsd.CCSD` / `ccsd_t` — CCSD and the (T) correction live on the
///   CC object, not on the SCF class
/// - `pyscf.mcscf.CASCI` / `CASSCF(mf, ncas, nelecas)` and `mrpt.NEVPT(mc)`
///
/// Each type holds a reference to a concrete mean-field method (`RHF`, `UHF`,
/// `ROHF`, `RKS`, …) and drives `run` itself. `MeanFieldMethod` keeps SCF
/// `run` / `forces` only.
#pragma once

#include "atomiverse/cpu/frame/frame.h"
#include "atomiverse/cpu/scf/casci.h"
#include "atomiverse/cpu/scf/casscf.h"
#include "atomiverse/cpu/scf/cc2.h"
#include "atomiverse/cpu/scf/ccsd.h"
#include "atomiverse/cpu/scf/fci.h"
#include "atomiverse/cpu/scf/gradient_result.h"
#include "atomiverse/cpu/scf/mp2.h"
#include "atomiverse/cpu/scf/mp3.h"
#include "atomiverse/cpu/scf/mp4.h"
#include "atomiverse/cpu/scf/nevpt2.h"
#include "atomiverse/cpu/scf/scf_guess.h"

namespace ATV {
namespace cpu {
namespace scf {

/// @brief Canonical MP2 on a mean-field reference.
///
/// Mirrors `pyscf.mp.mp2.MP2(mf)`: `frozen` and the Grimme OS/SS scales are
/// parameters of this object (`kernel` → `run`). Closed-shell `RHF`/`RKS` →
/// RMP2 or RI-MP2; `UHF` → UMP2; `ROHF` → ROMP2. `forces` is the nuclear
/// gradient of E_HF + E_corr (PySCF `nuc_grad_method().kernel()`).
/// @tparam Method Concrete mean-field type (`RHF`, `UHF`, `ROHF`, `RKS`, …).
template <typename Method>
class MP2 {
  public:
    /// @brief Bind a mean-field reference. The reference must outlive `*this`.
    /// @param mf Mean-field method (not copied).
    explicit MP2(Method& mf);

    /// @brief Lowest occupied spatial MOs excluded from the correlation sum.
    ///        PySCF `MP2(mf, frozen=N)`. 0 = correlate every occupied MO.
    int frozen = 0;

    /// @brief Opposite-spin scale (dimensionless). 1 = plain MP2; SCS = 1.2;
    ///        SOS = 1.3; B2PLYP = `hybrid_pt2_fraction`.
    double os_scale = 1.0;

    /// @brief Same-spin scale (dimensionless). 1 = plain MP2; SCS = 1/3;
    ///        SOS = 0; B2PLYP = `hybrid_pt2_fraction`.
    double ss_scale = 1.0;

    /// @brief Fresh SCF on `mf`, then the MP2 correction.
    /// @param frame Geometry [Å].
    /// @param guess Optional density guess for the reference SCF.
    /// @return Correlation and total energy [Hartree].
    MP2Result run(Frame& frame, const SCFInitialGuess* guess = nullptr);

    /// @brief Analytic nuclear forces on E_HF + E_corr [eV/Å].
    /// @param frame Geometry; forces written to fx/fy/fz, energy to meta.
    /// @param guess Optional density guess for the reference SCF.
    /// @return Forces [eV/Å] and total MP2 energy [eV].
    GradientResult forces(Frame& frame, const SCFInitialGuess* guess = nullptr);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell MP3. Mirrors a post-HF object on an RHF reference
///        (PySCF 2.x has no `mp.MP3`; the energy kernel is unchanged).
template <typename Method>
class MP3 {
  public:
    explicit MP3(Method& mf);
    /// @brief Frozen occupied spatial MOs (0 = full). Same convention as `MP2`.
    int frozen = 0;
    /// @param frame Geometry [Å].
    /// @param guess Optional density guess for the reference SCF.
    /// @return MP3 correlation pieces [Hartree].
    MP3Result run(Frame& frame, const SCFInitialGuess* guess = nullptr);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell MP4(SDTQ) on an RHF reference.
template <typename Method>
class MP4 {
  public:
    explicit MP4(Method& mf);
    /// @brief Frozen occupied spatial MOs (0 = full).
    int frozen = 0;
    /// @param frame Geometry [Å].
    /// @param guess Optional density guess for the reference SCF.
    /// @return MP4 correlation pieces [Hartree].
    MP4Result run(Frame& frame, const SCFInitialGuess* guess = nullptr);

  private:
    Method* mf_ = nullptr;
};

/// @brief Tiny full CI on an RHF reference (`pyscf.fci.FCI(mf)` shape).
template <typename Method>
class FCI {
  public:
    explicit FCI(Method& mf);
    /// @param frame Geometry [Å]. Drives a fresh SCF, then the CI.
    /// @return Full-CI total and correlation energy [Hartree].
    FCIResult run(Frame& frame);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell CASCI. Mirrors `pyscf.mcscf.CASCI(mf, ncas, nelecas)`.
template <typename Method>
class CASCI {
  public:
    explicit CASCI(Method& mf);
    /// @param frame Geometry [Å].
    /// @param ncas Active spatial orbitals.
    /// @param nelecas Active electrons (even, closed-shell).
    /// @return CASCI total and correlation energy [Hartree].
    CASCIResult run(Frame& frame, int ncas, int nelecas);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell CASSCF. Mirrors `pyscf.mcscf.CASSCF(mf, ncas, nelecas)`.
template <typename Method>
class CASSCF {
  public:
    explicit CASSCF(Method& mf);
    /// @param frame Geometry [Å].
    /// @param ncas Active spatial orbitals.
    /// @param nelecas Active electrons (even, closed-shell).
    /// @return CASSCF total and correlation energy [Hartree].
    CASSCFResult run(Frame& frame, int ncas, int nelecas);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell SC-NEVPT2. Drives the same RHF→CASSCF→PT2 sequence as
///        the previous `MeanFieldMethod::nevpt2` entry (PySCF
///        `mrpt.NEVPT(mc).kernel()` on that CASSCF reference).
template <typename Method>
class NEVPT2 {
  public:
    explicit NEVPT2(Method& mf);
    /// @param frame Geometry [Å].
    /// @param ncas Active spatial orbitals.
    /// @param nelecas Active electrons (even, closed-shell).
    /// @return NEVPT2 total and PT2 correlation [Hartree].
    NEVPT2Result run(Frame& frame, int ncas, int nelecas);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell CCSD. Mirrors `pyscf.cc.ccsd.CCSD(mf)`: `run` is
///        `kernel`, `ccsd_t` is `CCSD.ccsd_t`, `forces` is the nuclear
///        gradient. `diis` selects amplitude DIIS (false = plain Jacobi).
template <typename Method>
class CCSD {
  public:
    explicit CCSD(Method& mf);
    /// @brief Amplitude DIIS. Default true. False = Jacobi, same fixed point.
    bool diis = true;
    /// @brief Frozen-core gate for `forces` (energy path does not freeze).
    ///        Nonzero `forces` throws — frozen-core CCSD gradients are not
    ///        implemented.
    int frozen = 0;

    /// @param frame Geometry [Å]. Drives a fresh SCF, then the amplitude solve.
    /// @return CCSD correlation and total energy [Hartree].
    CCSDResult run(Frame& frame);

    /// @brief (T) correction. Mirrors `pyscf.cc.ccsd.CCSD.ccsd_t`.
    /// @param frame Geometry [Å].
    /// @return CCSD(T) correlation pieces [Hartree].
    CCSDTResult ccsd_t(Frame& frame);

    /// @brief Nuclear forces on E_HF + E_CCSD [eV/Å].
    /// @param frame Geometry; forces written to fx/fy/fz, energy to meta.
    /// @param guess Optional density guess for the reference SCF.
    /// @return Forces [eV/Å] and total CCSD energy [eV].
    GradientResult forces(Frame& frame, const SCFInitialGuess* guess = nullptr);

  private:
    Method* mf_ = nullptr;
};

/// @brief Closed-shell CC2. Mirrors PySCF `cc.RCCSD` with `cc2=True`.
template <typename Method>
class CC2 {
  public:
    explicit CC2(Method& mf);
    /// @param frame Geometry [Å].
    /// @return CC2 correlation and total energy [Hartree].
    CC2Result run(Frame& frame);

  private:
    Method* mf_ = nullptr;
};

} // namespace scf
} // namespace cpu
} // namespace ATV
