# Your Astronauts' Medications Are Hitting Every System You Study

**Bo Brothers** | Project Aletheia | bo@projectaletheia.org
ORCID: 0009-0002-7449-5459

---

## The finding

We scored 10 drugs from the ISS medical kit against 1.2 million FDA adverse event reports (FAERS) across seven space-medicine domains: vestibular, neuro-ocular, cardiovascular, immune, sleep, cognitive, and bone/muscle.

**Astronaut drugs are enriched in 6.3 of 7 space-relevant domains on average.** Random drug sets of the same size hit 4.2 domains (z = 2.73, p = 0.002). The result survives case-count matching (Wilcoxon p = 0.002).

This means the side-effect profiles of ISS medications overlap substantially with the physiological systems already stressed by spaceflight.

## Worst offenders

| Drug | Domains enriched | Highest enrichment |
|:---|:---:|:---|
| Pseudoephedrine | 7/7 | Cognitive 16.8x, bone/muscle 13.5x |
| Promethazine | 7/7 | Vestibular 2.8x, cardiovascular 2.7x |
| Modafinil | 7/7 | Cognitive 4.8x, sleep 3.5x |
| Ibuprofen | 7/7 | Cognitive 4.5x, bone/muscle 3.2x |
| Acetaminophen | 7/7 | Cognitive 3.2x, bone/muscle 3.0x |
| Zolpidem | 6/7 | Bone/muscle 3.9x, cognitive 3.6x |

Four of ten drugs hit the cognitive domain hardest (pseudoephedrine, ibuprofen, acetaminophen, modafinil). Cognitive is the highest-enriched domain overall at mean 4.3x. These are medications used routinely on ISS for congestion, pain, and fatigue.

## What this is not

This is not a clinical study. FAERS enrichment does not prove causation in microgravity. It identifies a systematic overlap between known terrestrial side-effect profiles and space-medicine domains that warrants controlled investigation.

## What we are proposing

A collaboration to test whether FAERS-predicted multi-domain drug effects correlate with actual crew symptom reports. The analysis requires:

1. **De-identified crew medication logs and symptom data** — to test whether drugs with higher FAERS cross-domain enrichment associate with more multi-system complaints in flight
2. **Pharmacokinetic modeling integration** — microgravity alters drug absorption, distribution, and metabolism; our terrestrial enrichment scores could serve as priors for Bayesian PK/PD models

The computational infrastructure (1.2M FAERS cases indexed by drug, domain, and co-occurrence) is built and available for collaboration.

## Why this matters

Current ISS medical kit selection focuses on primary indication efficacy. No systematic screen exists for whether a drug's side-effect profile overlaps with spaceflight-stressed physiological systems. As missions extend to Mars (2.5 years), cumulative multi-system drug effects become a crew safety question.

A drug that is 2-3x enriched in vestibular side effects on Earth could be significantly worse when the vestibular system is already destabilized by microgravity. The FAERS data lets us identify these overlaps before they become in-flight incidents.

## Contact

Bo Brothers
Project Aletheia (https://project-aletheia.vercel.app/)
bo@projectaletheia.org

All analysis code and data available on request. Methodology paper in preparation.
