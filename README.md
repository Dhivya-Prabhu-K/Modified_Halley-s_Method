# Comparison Codes for MHM

This repository contains the MATLAB codes used for comparing the performance of the proposed MHM with other existing methods for computing the zeros of orthogonal polynomials.

## Folder Structure

### `Comparison_codes_MHM`

This folder contains two subfolders:

---

## 1. `Table comparison`

This folder contains the codes used to generate tables comparing three different methods, including MHM, based on:

- CPU run time
- Relative error
- Total number of iterations

The comparisons are carried out separately for:

- Legendre polynomials
- Hermite polynomials

Accordingly, this folder contains two separate subfolders corresponding to these polynomial families.

---

## 2. `cpu_comparison`

This folder contains the MATLAB codes used to compare the CPU run time of three different methods, including MHM, for computing the zeros of:

- Legendre polynomials
- Hermite polynomials

The comparison focuses on the computational efficiency of the methods for large-scale computations.

---

## Software and Hardware Details

All codes were executed using MATLAB R2024b (64-bit) on a 64-bit Windows PC equipped with a 12th-generation Intel Core i7 processor and 16 GB of RAM.

---

## Notes

- The scripts are organized according to the type of comparison and polynomial family.
- Each folder contains the required MATLAB files to reproduce the numerical experiments reported in the work.
