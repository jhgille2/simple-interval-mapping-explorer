# Simple Interval Mapping Explorer

An interactive, dependency-free web application that visualizes the statistical
methodology of **simple interval mapping** (Haley–Knott regression) for QTL
detection in a simulated recombinant inbred line (RIL) population on a single
chromosome.

**Live app:** https://jhgille2.github.io/simple-interval-mapping-explorer/

## What it shows

1. **RIL marker genotypes** — one row per line, one column per marker, in
   centiMorgan position. Click a row to select that RIL.
2. **QTL genotype weights from flanking markers** — the core of interval
   mapping: *P*(QTL = Parent A allele | observed marker data) scanned along the
   chromosome. At an observed marker the probability is 0 or 1; between markers
   it is a weighted interpolation set by the flanking genotypes and their
   recombination distances. The selected RIL is highlighted, the population mean
   |2*P*−1| ("information content") is overlaid, and the RIL's true (unobserved)
   QTL genotype is marked.
3. **Haley–Knott LOD scan** — at each 1 cM grid position the phenotype is
   regressed on the expected QTL genotype *x̂* = 2*P*−1; LOD =
   (*n*/2)·log₁₀(RSS₀/RSS₁). An optional Churchill–Doerge permutation
   threshold gives the genome-wide significance line.
4. **Scan summary** — peak position/LOD, effect estimate ± SE vs. the true
   effect, R², and information content at the peak.

Controls cover population size (50–1000 RILs), chromosome length, marker
spacing (density), marker missing rate, QTL position, additive effect,
residual variance, and the random seed (every simulation is reproducible).

## Statistical model

- RILs by repeated selfing; genotype along the chromosome simulated as a
  Markov chain with switch probability *R* = 2*r*/(1+2*r*) between loci at
  per-meiosis recombination fraction *r* (Haldane–Waddington map expansion;
  Haldane map, no interference).
- Phenotype: *yᵢ* = μ + *a·xᵢ* + εᵢ, *xᵢ* ∈ {−1, +1}, εᵢ ~ N(0, σ²).
- Interval weights from the Lander–Botstein conditional probabilities; e.g.
  *P*(QTL=A | A,A) = (1−*R*₁)(1−*R*₂)/(1−*R*₁₂).
- Haley–Knott regression scan; permutation threshold per Churchill & Doerge.

See the "Statistical methodology" section in the app for the full equations.

## Run locally

Any static file server works, e.g.:

```sh
python3 -m http.server 8000
# open http://localhost:8000/
```

No build step, no dependencies.

## Files

- `index.html` — page structure and controls
- `css/styles.css` — styling
- `js/sim-core.js` — pure statistical core (simulation, weights, scan,
  permutations); also loadable in Node for testing
- `js/app.js` — UI wiring and canvas rendering

## References

- Lander ES, Botstein D (1989) Mapping Mendelian factors underlying
  quantitative traits using RFLP linkage maps. *Genetics* 121:185–199.
- Haley CS, Knott SA (1992) A simple regression method for mapping
  quantitative trait loci in line crosses using flanking markers.
  *Heredity* 69:315–324.
- Churchill GA, Doerge RW (1994) Empirical threshold values for quantitative
  trait mapping. *Genetics* 138:963–971.
- Broman KW, Sen Ś (2009) *A Guide to QTL Mapping in Experimental Crosses*.
  Springer.

## License

MIT.
