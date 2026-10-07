# Analysis

**DrawSpectra.h** — сохранение PDF со спектрами и кривыми BlastWave.

## Фит + PDF (всё в `run_fit.C`)

```bash
root -l 'macros/run_fit.C("config/phenix_uu.json", "handNoFit")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "handFit")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "global")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "full")'
```

PDF пишутся в `output/<system>/`.

## Только перерисовка из txt

```bash
root -l 'macros/run_fit.C("config/phenix_uu.json", "draw")'        # BWparams.txt
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawGlobal")'  # GlobalBWparams.txt
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams")'   # T(beta), T(Npart), beta(Npart)
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:individual:pos")'
```

## Графики параметров (`DrawParams.h`)

Режим **`drawParams`** читает `BWparams.txt` и/или `GlobalBWparams.txt` и строит три PDF:

- `BlastWave_T_beta.pdf` — T<sub>kin</sub>(β<sub>T</sub>)
- `BlastWave_T_Npart.pdf` — T<sub>kin</sub>(N<sub>part</sub>), log x
- `BlastWave_beta_Npart.pdf` — β<sub>T</sub>(N<sub>part</sub>), log x

Global fit — крупные красные кружки, подпись «Global fit». Легенда частиц — 2 столбца (π⁺|π⁻, K⁺|K⁻, p|p̄). Название системы и `#sqrt{s_{NN}}` — отдельными строками внизу слева.

Подрежимы (через двоеточие после `drawParams`):

| Спецификация | Смысл |
|---|---|
| *(пусто)* / `both` | individual + global, все 6 частиц |
| `individual` | только per-particle fit |
| `global` | только global fit |
| `:pos` | π⁺, K⁺, p |
| `:neg` | π⁻, K⁻, p̄ (+ global negative) |
| `:pi`, `:k`, `:p` | только выбранная пара частиц |
| `:nopredict` | без кривой T(β) penalty |

Примеры:

```bash
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:global")'
root -l 'macros/run_fit.C("config/phenix_uu.json", "drawParams:both:pos")'
```

Для T(N<sub>part</sub>) и β(N<sub>part</sub>) в JSON нужен массив **`npart`** той же длины, что и `centrality` (значения из PHENIX `def.h`).

## Важно

- `handNoFit` — кривые из **`hand_params` в JSON**, не из `BWparams.txt`
- `handFit` / `full` — кривые из **результатов фита** (`pipeline.IndividualResults()`)
- `draw` — перечитать `BWparams.txt` и нарисовать заново
- Лимиты оси Y на PDF — секция **`draw`** в config (`y_min`, `y_max`, 6 значений: π⁺, π⁻, K⁺, K⁻, p, p̄)
- Множитель по центральности при отрисовке: **×10^(N-1-c)** — 0–20% (c=0) сверху; фит и точки умножаются одинаково

Планируется: T(centrality), T(beta), контуры χ².
