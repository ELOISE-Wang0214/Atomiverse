# Offline pin — finite-field-single-electrode

No live MetalWalls run. Analytic pin documented in
`regressions/md/finite-field-single-electrode/REFERENCE.md`.

```python
import math
coulomb = 14.399645
eta, Ez, z, L = 2.0, 0.1, 10.0, 20.0
Aii = coulomb * eta * math.sqrt(2 / math.pi)
Q = (Ez * z) / Aii
Pz = Q * z / (L ** 3)
print("Aii", Aii)
print("Q", Q)
print("Pz", Pz)
```

Expected:

- `Aii = 22.978508853084353` V/e  
- `Q = 0.04351892485250505` e  
- `Pz = 5.439865606563132e-05` e/Å²  
