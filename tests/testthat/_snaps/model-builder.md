# model printout snapshots

    Code
      print(buildCacoaModel(meta, ~ group + batch, test = "group"))
    Output
      Model: ~group + batch      dispersion: ~group
      Test: group: B vs A  (reference 'A': most frequent level; to change: test = "group: A vs B")
             adjusted for batch; permutations: block within 2 strata (400 distinct)
             shift > 0: B samples differ from A samples in a common direction, beyond within-group variability
      Samples: 12 used.   Issues: 1 warning (see $issues)

---

    Code
      print(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: G2 vs G1"))
    Output
      Model: ~Group + Batch      dispersion: ~Group
      Test: Group: G2 vs G1  (reference 'G1': as requested)
             adjusted for Batch; permutations: block within 2 strata (400 distinct)
             shift > 0: G2 samples differ from G1 samples in a common direction, beyond within-group variability
      Samples: 18 used.   Issues: none

---

    Code
      print(buildCacoaModel(meta3, ~ Group + Batch, test = "Group: all"))
    Output
      Model: ~Group + Batch      dispersion: ~Group
      Test: Group (3 levels)
             adjusted for Batch; permutations: block within 2 strata (2,822,400 distinct)
             location: the levels of Group differ in where their samples sit (a common direction per level), beyond within-level variability
      Samples: 18 used.   Issues: none

---

    Code
      print(buildCacoaModel(meta, ~ group * batch, test = "group"))
    Output
      Model: ~group * batch      dispersion: ~group
      Test: group: B vs A  (reference 'A': most frequent level; to change: test = "group: A vs B")
             adjusted for batch; permutations: block within 2 strata (400 distinct)
             shift > 0: B samples differ from A samples in a common direction, beyond within-group variability
             note: group interacts with batch: compared marginally (equal weights over batch levels); use a structured test with at = or over = to change
      Samples: 12 used.   Issues: 1 warning, 1 note (see $issues)

---

    Code
      print(buildCacoaModel(meta, ~ age + batch, test = "age"))
    Output
      Model: ~age + batch      dispersion: ~age
      Test: age: per 1 unit
             adjusted for batch; permutations: block within 2 strata (518,400 distinct)
             shift > 0: samples move in a common direction as age increases (per 1 unit)
      Samples: 12 used.   Issues: 2 warnings (see $issues)

---

    Code
      print(buildCacoaModel(mm, ~ group + site + batch, test = "group"))
    Output
      Model: ~group + site + batch      dispersion: ~group
      Test: group: B vs A  (reference 'A': most frequent level; to change: test = "group: A vs B")
             adjusted for site, batch; permutations: block within 2 strata (400 distinct)
             shift > 0: B samples differ from A samples in a common direction, beyond within-group variability
             note: term site dropped from the formula: constant in the samples used
      Samples: 12 used.   Issues: 1 warning, 1 note (see $issues)

