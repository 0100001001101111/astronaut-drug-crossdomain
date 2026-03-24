# Astronaut Drug Cross-Domain Enrichment Analysis

ISS medical kit drugs produce side effects that overlap with the physiological systems already stressed by spaceflight. This analysis quantifies that overlap.

## What this found

We scored 10 drugs from the ISS medical kit against 1.2 million FDA adverse event reports (FAERS 2023Q4) across seven space-medicine domains: vestibular, neuro-ocular, cardiovascular, immune, sleep, cognitive, and bone/muscle.

**Astronaut drugs are enriched in 6.3 of 7 space-relevant domains on average.** Random drug sets hit 4.2 domains (z = 2.73, p = 0.002). The result survives case-count matching (Wilcoxon p = 0.002).

Key findings:
- **Pseudoephedrine** (decongestant): 16.8x cognitive enrichment, 13.5x bone/muscle — worst offender
- **Promethazine** (anti-nausea, most commonly used space drug): enriched in all 7 domains at 2–3x each
- **4 of 10 drugs** hit the cognitive domain hardest; cognitive is the highest-enriched domain overall (mean 4.3x)
- **Modafinil** (wakefulness): 4.8x cognitive, 3.5x sleep enrichment

## Why this matters

Current ISS medical kit selection focuses on primary indication efficacy. No systematic screen exists for whether a drug's side-effect profile overlaps with spaceflight-stressed systems. As missions extend to Mars (2.5 years), cumulative multi-system drug effects become a crew safety question.

## Reproduce every number

```bash
pip install -r requirements.txt
python verify_and_analyze.py
```

This single script:
- Defines all 10 astronaut drugs with FAERS enrichment profiles (inline, no external files)
- Defines MedDRA preferred terms for each space-medicine domain
- Runs permutation test (100,000 iterations, seed 42) for the z-score
- Runs case-count-matched Wilcoxon signed-rank test
- Generates Figure 1 (heatmap + null distribution) and enrichment table
- Prints every claim-relevant number to stdout

## Methods

### Data source
FDA Adverse Event Reporting System (FAERS) 2023Q4 quarterly extract. 1,213,478 reports with adverse event and drug information.

### Domain mapping
MedDRA preferred terms were mapped to seven NASA Human Research Program domains. For example, "vestibular" includes VERTIGO, DIZZINESS, SYNCOPE, BALANCE DISORDER, etc. Full term lists are defined in the script.

### Enrichment calculation
For each drug-domain pair: (rate of domain AEs among drug reports) / (rate of domain AEs across all FAERS reports). Values >1.5x are considered enriched (standard pharmacovigilance threshold).

### Statistical tests
1. **Permutation test**: 10,000 random drug sets of size 10 drawn from full FAERS drug list. Mean enriched domains per set forms null distribution. Astronaut drugs z = 2.73, p = 0.002.
2. **Case-count matching**: Each astronaut drug matched to random drug with similar FAERS report count. Controls for drugs with more reports having more detectable enrichment. Wilcoxon p = 0.002.

### What this analysis does NOT do
- This is not a clinical study of astronauts
- FAERS enrichment ≠ causation in microgravity
- Individual crew data (LSAH) is not publicly available
- The analysis identifies systematic overlap; validation requires crew symptom data

## Output

- `output/figure1.png` — Enrichment heatmap + permutation null distribution (300 DPI)
- `output/figure1.pdf` — Vector format
- `output/enrichment_table.md` — Complete drug × domain enrichment table

## Software

Tested with Python 3.14.3, scipy 1.17.1, numpy 2.4.3. Compatible with Python 3.10+.

## License

CC-BY 4.0

## Contact

Bo Brothers
Project Aletheia
bo@projectaletheia.org
https://project-aletheia.vercel.app/
