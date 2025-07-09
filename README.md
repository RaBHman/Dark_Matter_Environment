# 🛰️ Post-Adiabatic Waveforms from Extreme Mass Ratio Inspirals in the Presence of Dark Matter

## 👥 Authors:
**Mostafizur Rahman**, **Takuya Takahashi**

This repository contains an accompanying notebook that supports the research presented in the above work. It includes symbolic and numerical computations for generating post-adiabatic waveforms from extreme mass ratio inspirals (EMRIs), accounting for the influence of surrounding dark matter.

The notebook implements the perturbative source terms, master variables, and coefficient structures required to evolve the system and extract physically meaningful waveforms.

---

## 📘 Source Term Coefficients and Master Variables

This code computes source terms and coefficients relevant to perturbative analyses in gravitational and fluid systems, particularly for axial and polar modes. The symbolic notation is simplified in the code for readability and consistency.

---

## 🔢 Notation Conventions

- **Coefficient Naming:**
  - `A0lm` → `Subsuperscript[A, 0, lm, (1,0)]`
  - `δA0lm` → `Subsuperscript[A, 0, lm, (1,1)]`

---

## 🧩 Master Variables

The following variables represent master functions used in the perturbative formalism:

| Variable     | Description                           | Symbolic Representation               |
|--------------|---------------------------------------|----------------------------------------|
| `R[r]`       | Axial master variable (order ε)       | `Subsuperscript[Ψ, R, (1,0)]`          |
| `δR[r]`      | Axial master variable (order εξ)      | `Subsuperscript[Ψ, R, (1,1)]`          |
| `Kh[r]`      | Polar master variable (order ε)       | `Subsuperscript[Ψ, Z, (1,0)]`          |
| `δKh[r]`     | Polar master variable (order εξ)      | `Subsuperscript[Ψ, Z, (1,1)]`          |
| `XF[r]`      | Polar fluid master variable (order εξ)| `Subsuperscript[Ψ, F, (1,1)]`          |

Additional parameters:

- `n1 = n`
- `r01v1 = Subscript[r, (0,1)]`

---

## 🧾 Source Term Matrix Entries

Below is the list of entries in the source term matrix. Each item corresponds to a specific perturbative source term or coefficient:

- **Source term at order ε (axial gravity):** `Subsuperscript[S, R, (1,0)]`
- **Source term at order εξ (axial gravity):** `Subsuperscript[S, R, (1,1)]`
- **Source term at order ε (polar gravity):** `Subsuperscript[S, Z, (1,0)]`
- **Source term at order εξ (polar fluid):** `Subsuperscript[S, F, (1,1)]`
- `Subsuperscript[C, 1, lm, (1,0)]`
- `Subsuperscript[C, 2, lm, (1,0)]`
- `Subsuperscript[C, 1, lm, (1,1)]`
- `Subsuperscript[C, 2, lm, (1,1)]`
- **Source term at order εξ (polar gravity):** `Subsuperscript[S, Z, (1,1)]`

