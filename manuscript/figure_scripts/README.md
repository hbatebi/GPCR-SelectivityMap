# Figure generation

Regenerates the four main manuscript figures from the exported source tables.

```bash
cd manuscript/figure_scripts
python fig1.py && python fig2.py && python fig3.py && python fig4.py
```

Each script writes PDF, SVG, and PNG to `../figures/`. `style.py` holds the shared
rcParams and colour palette; `receptor.py` holds the vector receptor primitives used
before the hand-drawn schematics replaced panels 1A and 2B.

Panels 1A and 2B embed raster schematics (`../figures/Figure_1A_schematic_composite.png`,
`../figures/Figure_2B_receptor_schematic.png`); all other panels are fully vector.
