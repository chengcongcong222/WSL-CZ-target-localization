# MATURE_MODE_SOLVER_NOTE

UTC: 2026-09-23T07:32:46.685640+00:00

- Tool: Acoustics Toolbox KRAKEN (OALib atWin10_2018_7 prebuilt)
- Path: `C:\Users\ccc\XiaomiMiMoProjects\WSL系统汇聚区目标定位\tools\acoustics_toolbox\atWin10\at\bin\kraken.exe`
- Self-built FEM solver: **NOT used / NOT ADMISSIBLE for final claims**
- Scene: Munk-like SSP, H=5 km, zs in {180,200,220}, zr=200, f in {201,235,283,338} Hz
- BC: vacuum top (V), fluid/rigid bottom (documented in env)
- Mode selection: propagating modes from KRAKEN .mod (complex k, phi)
- Yang formulas: **YANG_FORMULA_RECOVERY_PARTIAL** (SA from mode sum, modal peak→depth via φ(z)φ(zr), δ=0)
- δ=0 paper-like ideal this round
