# DGmini, a minimal 1D discontinuous Galerkin solver

Light-weight 1D DG reference C++ implementation intended for experimentation and understanding of DG building blocks.

---

## Build & Run

Manually build and run with:

```bash
cmake -S . -B build
cmake --build build
./build/dgmini            # To run the solver
./build/dgmini_unittest   # To run unittests
```

**VS Code interface**: Following the definitions in `.vscode/tasks.json`, one can also:

1) Clone repo
2) Open `DGmini` folder in VS Code
3) Press F1, select `Tasks: Run Task`, select `CMake Configure`
4) Press F1, select `Tasks: Run Task`, select `Build dgmini`
5) Press F1, select `Tasks: Run Task`, select `Run dgmini`
6) Press F1, select `Tasks: Run Task`, select `Build dgmini_unittest`
7) Press F1, select `Tasks: Run Task`, select `Run dgmini_unittest`

---