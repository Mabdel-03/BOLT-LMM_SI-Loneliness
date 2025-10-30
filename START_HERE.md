# 🚀 START HERE: BOLT-LMM Social Isolation/Loneliness Analysis

**Quick start guide for running the complete BOLT-LMM analysis pipeline**

---

## What This Analysis Does

- **3 social isolation phenotypes**: Loneliness, FreqSoc (frequency of social contact), AbilityToConfide
- **3 population stratifications**: EUR_MM, EUR_Male, EUR_Female
- **1 covariate model**: Day_NoPCs (age, sex, array)
- **Total**: 9 GWAS analyses

---

## Three-Step Workflow

### Step 1: Filter to Populations (~10 min)

```bash
cd /home/mabdel03/data/files/Isolation_Genetics/GWAS/Scripts/ukb21942/BOLT-LMM_SI-Loneliness
sbatch 0a_filter_populations.sbatch.sh
tail -f 0a_filter.out  # Monitor
```

**Creates**: 6 filtered files (3 phenotype + 3 covariate)

### Step 2: Test Run (~1-2 hours)

```bash
sbatch 0b_test_run.sbatch.sh
tail -f 0b_test.out  # Monitor
```

**Tests**: Loneliness phenotype with EUR_MM population  
**⚠️ MUST PASS before Step 3**

### Step 3: Full Analysis (~2 hours, 9 jobs)

```bash
sbatch 1_run_bolt_lmm.sbatch.sh
squeue -u mabdel03  # Monitor all jobs
```

**Runs**: All 9 phenotype-population combinations

---

## Job Mapping

```
Task 1-3: Loneliness        × [EUR_MM, EUR_Male, EUR_Female]
Task 4-6: FreqSoc           × [EUR_MM, EUR_Male, EUR_Female]
Task 7-9: AbilityToConfide  × [EUR_MM, EUR_Male, EUR_Female]
```

---

## Expected Output

```
results/Day_NoPCs/
├── EUR_MM/
│   ├── bolt_Loneliness.Day_NoPCs.stats.gz
│   ├── bolt_Loneliness.Day_NoPCs.log.gz
│   ├── bolt_FreqSoc.Day_NoPCs.stats.gz
│   ├── bolt_FreqSoc.Day_NoPCs.log.gz
│   ├── bolt_AbilityToConfide.Day_NoPCs.stats.gz
│   └── bolt_AbilityToConfide.Day_NoPCs.log.gz
    ├── EUR_Male/
    │   └── [same 6 files]
    └── EUR_Female/
        └── [same 6 files]
```

**Total**: 18 files (9 stats.gz + 9 log.gz)

---

## Quick Verification

```bash
# Check all jobs completed
ls -1 results/Day_NoPCs/*/bolt_*.stats.gz | wc -l
# Should be: 9

# Quick QC for one phenotype
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep "Analyzing"
zcat results/Day_NoPCs/EUR_MM/bolt_Loneliness.Day_NoPCs.log.gz | grep "h2:"
```

---

## Expected Results (Day et al. 2018)

- **Loneliness**: h² ≈ 3-5%, multiple genome-wide significant loci
- **FreqSoc**: h² ≈ 4-6%, associations with personality genes
- **AbilityToConfide**: h² ≈ 3-5%, overlap with mental health loci

---

## Key Difference from Original Study

**Sex Stratification**: This analysis includes separate male and female analyses to detect sex-specific genetic effects on social isolation traits.

---

## Need Help?

- **Full documentation**: See `README.md`
- **Troubleshooting**: Check SLURM logs (`*.out` and `*.err` files)
- **QC guidelines**: See README.md "Quality Control" section

---

**Ready?** Run Step 1: `sbatch 0a_filter_populations.sbatch.sh`

*Last Updated: October 30, 2025*

