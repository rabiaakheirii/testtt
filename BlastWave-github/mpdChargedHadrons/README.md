## Установить pyROOT и запустить postprocess

```
PYTHONPATH=/home/dasha/NICA/mpdGit/venv/bin/python
PYTHONPATH=$PYTHONPATH:$ROOTSYS/lib
source thisroot.sh (root6)
pip install rootpy
python postprocess.py -i 'input/mpdpid10.root' -d output -s input/settings.json
```

## Файлы

 __input/train_output/\*__ - выходные файлы паровоза Киреева, взаты из /junk/kireev/paper_mpdpid*.root

**input/postprocess_mpdpid10.root** - постпроцесс файлов паровоза Киреева из train_output

**input/nuclei_spectra** - то, что я сама запускала раньше на маленькой статистике (на 1 файле.)

**ChemicalPotential** - строит фазовую диаграмму и отношения

**BlastWaveGlobal** - строит BlastWave GlobalFit. 
output - BlastWaveGlobalFit.pdf, GlobalBWparams.txt

**BlastWaveGlobalAllParts** - строит BlastWave GlobalFit одновременно по всем частицам. Плохо работает

**BlastWave.C** - строит отдельные BlastWave фиты для каждой частицы, а в качетсве начальных параметров аппроксимации задаются параметры, полученные в ГлобалФите (output/txtParams/GlobalBWparams.txt)

## BlastWave  12.2024

**BlastWaveGlobal** - строит Глобальный (общий) фит по pi, K, p одновременно. 

* output/txtParams/BlastWaveGlobalFit.pdf - спектры с глобал фитом
* output/txtParams/GlobalBWparams.txt - параметаы глобал фита (charge, centr, T, constPi, constK, constP)

**BlastWave.C** - строит отдельные BlastWave фиты для каждой частицы, а в качетсве начальных параметров аппроксимации задаются параметры, полученные в ГлобалФите (output/txtParams/GlobalBWparams.txt)

* output/BlastWave.pdf - спектры с отдальными фитами
* output/txtParams/BWparams.txt - параметаы фитов (part, centr, const, T, Terr, ut, ut_err)
* output/txtParams/BlastWave_contour.pdf - контурные графики (не очень работают)

**BlastWaveFit.h** - фитирует БластВейвом различными способами (разные варианты задания начальных параметров). По умолчанию - 0.

**BlastWaveSyst.C** - строит систематические погрешности для фитирования.
* output/txtParams/BWparamsSyst.txt 

**BWDrawParams** - строит T и ut как функции от центральности по результатам output/txtParams/BWparams.txt

* output/BWparam_T.pdf 
* output/BWparam_ut.pdf

**WriteReadFiles.h** - запись считывание из txt файлов 

## Cumulative 12.12.2024

**Cumulative.h** - void DrawCumulativeBorder( int part, double pad_min, double pad_max ) - строит вертикальную линию по кумулятивной границе. Эта граница добавлена в **spectra.C** и на spectra_postprocess_mpdpid10.pdf