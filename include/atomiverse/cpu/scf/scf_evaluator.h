/// @file cpu/scf/scf_evaluator.h
/// @brief Compatibility include — use concrete methods from methods.h.
///
/// Historical `cpu::scf::SCF` was a method-agnostic bag driven by `SCFConfig`.
/// The public API is now `RHF` / `UHF` / `RKS` / `UKS` (see methods.h).
/// `SCF` is a typedef of `RKS` (production default: closed-shell KS/PBE).
#pragma once

#include "atomiverse/cpu/scf/methods.h"
